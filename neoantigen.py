#!/opt/common/CentOS_6-dev/python/python-2.7.10/bin/python
import os, sys, subprocess, psutil
import argparse, re, errno
import time, datetime, dateutil.parser
import traceback
import pandas as pd
import logging, logging.handlers
import gzip
import copy
from ConfigParser import ConfigParser
from joblib import Parallel, delayed

#####
# Neoantigen prediction pipeline. Four main steps: 
#       (1) Genotype HLA using POLYSOLVER, 
#       (2) Constructed mutated peptide sequences from HGVSp/HGVSc 
#       (3) Run NetMHC-4.0 + NetMHCpan-4.0
#       (4) Gather results and generate: 
#               - <sample_id>.neoantigens.maf: original maf with neoantigen prediction columns added
#               - <sample_id>.all_neoantigen_predictions.txt: all the predictions made for all peptides by both the algorithms
#####

def main():
    prog_description = (
        '# Neoantigen prediction pipeline. Four main steps: \n'
        '\t\t(1) Genotype HLA using POLYSOLVER, \n'
        '\t\t(2) Constructed mutated peptide sequences from HGVSp/HGVSc \n'
        '\t\t(3) Run NetMHC-4.0 + NetMHCpan-4.0 \n'
        '\t\t(4) Gather results and generate: \n'
        '\t\t\t\t- <sample_id>.neoantigens.maf: original maf with neoantigen prediction columns added \n'
        '\t\t\t\t- <sample_id>.all_neoantigen_predictions.txt: all the predictions made for all peptides by both the algorithms \n'
    )
    prog_epilog = (
        '\n'
    )

    parser = argparse.ArgumentParser(description=prog_description,
                                     epilog=prog_epilog,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     add_help=True)
    required_arguments = parser.add_argument_group('Required arguments')
    required_arguments.add_argument('--config_file',
                        required=True,
                        help='See: neoantigen-luna.config in the repo')
    required_arguments.add_argument('--sample_id',
                        required=True,
                        help='sample_id used to limit neoantigen prediction to identify mutations '
                             'associated with the patient in the MAF (column 16). ')
    required_arguments.add_argument('--output_dir',
                        required=True,
                        help='output directory')
    required_arguments.add_argument('--maf_file',
                        required=True,
                        help='expects a CMO maf file (post vcf2maf.pl)')
    required_arguments.add_argument('--normal_bam',
                        required=False,
                        help='full path to normal bam file. Either --normal_bam or --hla_file are required.')

    optional_arguments = parser.add_argument_group('Optional arguments')                        
    optional_arguments.add_argument('--hla_file',
                        required=False,
                        help='POLYSOLVER output file (winners.hla.txt) for the sample. If not provided,'
                             'POLYSOLVER is run. Either --normal_bam or --hla_file are required.')
    optional_arguments.add_argument('--peptide_lengths', 
                        required=False,
                        help='comma-separated numbers indicating the lengths of peptides to generate. Default: 9,10')
    optional_arguments.add_argument('--keep_tmp_files',
                        required=False,
                        help='keeps POLYSOLVER\'s temporary files. for debugging purposes. ' +
                             'Note: TMP files can be more than 5GB. Default: true',
                        action='store_true')
    optional_arguments.add_argument('--force_rerun_polysolver',
                        required=False,
                        help='ignores any existing polysolver output and re-runs it. Default: false',
                        action='store_true')
    optional_arguments.add_argument('--force_rerun_netmhc',
                        required=False,
                        help='ignores any existing netMHCpan output and re-runs it. Default: false',
                        action='store_true')                        


    args = parser.parse_args()

    normal_bamfile = str(args.normal_bam)
    maf_file = str(args.maf_file)
    output_dir = str(args.output_dir)
    hla_file = str(args.hla_file)
    sample_id = str(args.sample_id)
    
    peptide_lengths = [9, 10]
    if args.peptide_lengths is not None:
        peptide_lengths = map(int, str(args.peptide_lengths).split(','))

    if args.config_file is None:
        config_file_path = os.path.dirname(os.path.realpath(__file__)) + '/neoantigen-luna.config'
    else:
        config_file_path = str(args.config_file)

    if not os.path.exists(config_file_path):
        print 'Error: could not open config file: ' + config_file_path + '. Exiting.'
        exit(1)

    keep_tmp_files = False
    force_polysolver = False
    force_netmhc = False

    if args.keep_tmp_files is not None and args.keep_tmp_files:
        keep_tmp_files = True

    if args.force_rerun_polysolver is not None and args.force_rerun_polysolver is True:
        force_polysolver = True

    if args.force_rerun_netmhc is not None and args.force_rerun_netmhc is True:
        force_netmhc = True

    if args.normal_bam is None and args.hla_file is None:
        print >> sys.stderr, 'Error: --normal_bam or --hla_file is required. Exiting.'
        exit(1)

    if args.normal_bam is not None and not os.path.exists(normal_bamfile):
        print >> sys.stderr, 'Error: --normal_bam '' + normal_bamfile + '' does not exist. Exiting.'
        exit(1)

    if args.hla_file is not None and not os.path.exists(hla_file):
        print >> sys.stderr, 'Error: --hla_file '' + hla_file + '' does not exist. Exiting.'
        exit(1)

    if not os.path.exists(maf_file):
        print >> sys.stderr, 'Error: --maf_file '' + maf_file + '' does not exist. Exiting.'
        exit(1)

    os.system('mkdir -p ' + output_dir)

    config = ConfigParser()
    config.read(config_file_path)

    #
    # initialize loggers
    #
    logger = logging.getLogger('neoantigen')
    logger.setLevel(logging.DEBUG)
    console_formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(message)s',
                                          datefmt='%m-%d-%Y %H:%M:%S')

    # logfile handler
    handler_file = logging.FileHandler(output_dir + '/neoantigen_run.log', mode='w')
    handler_file.setLevel(logging.DEBUG)
    handler_file.setFormatter(console_formatter)
    logger.addHandler(handler_file)

    # stdout handler
    handler_stdout = logging.StreamHandler(sys.stdout)
    handler_stdout.setFormatter(console_formatter)
    handler_stdout.setLevel(logging.INFO)
    logger.addHandler(handler_stdout)

    logger.info('Starting neoantigen prediction run')
    logger.info('\tLog file: ' + output_dir + '/neoantigen_run.log')
    logger.info('\t--config_file:' + config_file_path)
    logger.info('\t--normal_bam: ' + normal_bamfile)
    logger.info('\t--maf_file: ' + maf_file)
    logger.info('\t--output_dir: ' + output_dir)

    ### Load from config file
    reference_cdna_file = config.get('Reference Paths', 'GRCh37cdna')
    reference_cds_file = config.get('Reference Paths', 'GRCh37cds')
    netmhc4_bin = config.get('NetMHCpan', 'netmhc4_bin_path')
    netmhc4_alleleslist = config.get('NetMHCpan', 'netmhc4_alleleslist')
    netmhcpan4_bin = config.get('NetMHCpan', 'netmhcpan4_bin_path')

    if not os.path.exists(reference_cdna_file) or not os.path.exists(reference_cds_file):
        sys.exit('Could not find reference CDNA/CDS files: ' + '\n\t' + reference_cdna_file + '\n\t' + reference_cds_file)

    polysolver_bin = config.get('POLYSOLVER', 'polysolver_bin')

    #######################
    ### POLYSOLVER
    #######################
    hla_alleles = []
    if args.hla_file is None or not os.path.exists(hla_file):
        try:
            hla_file = output_dir + '/polysolver/winners.hla.txt'
            if not force_polysolver:
                logger.info('Checking for previous run of POLYSOLVER...')

            if os.path.isfile(os.path.abspath(hla_file)) and os.stat(
                    os.path.abspath(hla_file)) > 0 and not force_polysolver:
                logger.info('Previous run found! (using hla alleles in ' + hla_file + ')')
            else:
                if not force_polysolver:
                    logger.info('Running POLYSOLVER...')
                else:
                    logger.info('Running POLYSOLVER... (in --force_rerun_polysolver mode)')

                run_polysolver_cmd = 'export NUM_THREADS=8; ' + polysolver_bin + ' ' + normal_bamfile + \
                                     ' Unknown 1 hg19 STDFQ 0 ' + output_dir + '/polysolver/' + ' > ' + \
                                     output_dir + '/polysolver.log ' + ' 2> ' + output_dir + '/polysolver.err '

                execute_cmd(run_polysolver_cmd)

                if not keep_tmp_files:
                    logger.info('Cleaning up POLYSOLVER run directory')
                    cleanup_cmd = 'rm -rf ' + output_dir + '/polysolver/*fastq ' + \
                                  output_dir + '/polysolver/*bam ' + \
                                  output_dir + '/polysolver/*bai '
                    execute_cmd(cleanup_cmd)

        except Exception:
            logger.error('Could not run POLYSOLVER. Check polysolver.log and polysolver.err files in ' + output_dir)
            logger.error(traceback.format_exc())
            if not keep_tmp_files:
                logger.info('Cleaning up POLYSOLVER run directory')
                cleanup_cmd = 'rm -rf ' + output_dir + '/polysolver/*fastq ' + \
                              output_dir + '/polysolver/*bam ' + \
                              output_dir + '/polysolver/*bai '
                execute_cmd(cleanup_cmd)
            exit(1)

    ## parse hla-alleles into the format that is required by NetMHC
    if os.path.isfile(os.path.abspath(hla_file)):
        for allele in re.split('\n|\t', subprocess.check_output('cut -f 2-3 ' + hla_file, shell=True)):
            if allele == '':
                continue
            levels = allele.split('_')
            hla_alleles.append('HLA-' + levels[1].upper() + levels[2] + ':' + levels[3])

    #######################
    ### FASTA with mutated peptides
    #######################
    logger.info('Generating mutated peptides FASTA')
    sample_path_pfx = output_dir + '/' + sample_id
    mutated_sequences_fa = sample_path_pfx + '.mutated_sequences.fa'
    mutations = []
    out_fa = open(mutated_sequences_fa, 'w')

    ## generate .debug.fa for debugging purposes. 
    debug_out_fa = open(sample_path_pfx + '.mutated_sequences.debug.fa', 'w')

    try:
        logger.info('Loading reference CDS/cDNA sequences...')
        cds_seqs = load_transcript_fasta(reference_cds_file)
        cdna_seqs = load_transcript_fasta(reference_cdna_file)
        logger.info('Finished loading reference CDS/cDNA sequences...')
        
        logger.info('Reading MAF file and constructing mutated peptides...')
        maf_df = pd.read_table(maf_file, comment='#', low_memory=False, header=0)
        n_muts = n_non_syn_muts = n_missing_tx_id = 0
        for index, row in maf_df.iterrows():
            cds_seq = ''
            cdna_seq = ''
            
            n_muts += 1
            tx_id = row['Transcript_ID']
            if tx_id in cds_seqs:
                cds_seq = cds_seqs[tx_id]

            if tx_id in cdna_seqs:
                cdna_seq = cdna_seqs[tx_id]

            mut = mutation(row, cds_seq, cdna_seq)

            if mut.is_non_syn():
                n_non_syn_muts += 1

            if cds_seq == '':
                n_missing_tx_id += 1
                
            if cds_seq != '' and mut.is_non_syn():
                mut.generate_translated_sequences(max(peptide_lengths))         

            if len(mut.mt_altered_aa) > 5:
                out_fa.write('>' + mut.identifier_key + '_mut\n')
                out_fa.write(mut.mt_altered_aa + '\n')
                
                ### write out WT/MT CDS + AA for debugging purposes
                debug_out_fa.write('>' + mut.identifier_key + '_mut\n')
                debug_out_fa.write('mt_altered_aa: ' + mut.mt_altered_aa + '\n')
                debug_out_fa.write('wt_full_cds: ' + mut.wt_cds + '\n')
                debug_out_fa.write('wt_full_aa: ' + mut.wt_aa + '\n')
                debug_out_fa.write('mt_full_cds: ' + mut.mt_cds + '\n')
                debug_out_fa.write('mt_full_aa: ' + mut.mt_aa + '\n')
                

            mutations.append(mut)
        out_fa.close()
        debug_out_fa.close()

        logger.info('\tMAF mutations summary')
        logger.info('\t\t# mutations: ' + str(n_muts))
        logger.info('\t\t# non-syn: ' + str(n_non_syn_muts) + ' (# with missing CDS: ' + str(n_missing_tx_id) + ')')

        # create empty neoantigens.maf file if there are no mutations 
        if n_muts == 0:
            logger.info('No mutations in this tumor. Creating empty .neoantigens.maf file')
            logger.info('Exiting neoantigen pipeline')
            execute_cmd('touch ' + sample_path_pfx + '.neoantigens.maf')
            exit(1)

    except Exception:
        logger.error('Error while generating mutated peptides')
        logger.error(traceback.format_exc())
        exit(1)

    #######################
    ### Run NetMHCpan 4.0, netMHC 4.0
    #######################
    try:
        logger.info('')
        logger.info('Running MHC--peptide binding predictions using NetMHCpan 4.0...')
    
        netmhcpan_output_pfx = sample_path_pfx + '.netmhcpan4.output'
        if not force_netmhc and check_file_exists_and_not_empty(netmhcpan_output_pfx + '.txt'):
            logger.info('Previous run of NetMHCpan-4.0 found... Skipping!')
        else:
            #####
            logger.info('Starting NetMHCpan 4.0...')
            #####
            logger.info('Predicting on the following HLA-alleles: ' + ','.join(sorted(set(hla_alleles))))
            run_netmhcpan_cmd = 'rm -f ' + netmhcpan_output_pfx + '.txt; '+ \
                                netmhcpan4_bin + ' -s -BA ' + \
                                ' -a ' + ','.join(sorted(set(hla_alleles))) + \
                                ' -f ' + mutated_sequences_fa + \
                                ' -l ' + ','.join(map(str, peptide_lengths)) + \
                                ' -inptype 0 ' + \
                                ' -xls ' + \
                                ' -xlsfile ' + netmhcpan_output_pfx + '.xls ' + \
                                ' > ' + netmhcpan_output_pfx + '.txt'
            execute_cmd(run_netmhcpan_cmd)

        netmhc_output_pfx = sample_path_pfx + '.netmhc4.output'
        if not force_netmhc and check_file_exists_and_not_empty(netmhc_output_pfx + '.txt'):
            logger.info('Previous run of NetMHC-4.0 found... Skipping!')
        else:
            #####
            logger.info('Starting NetMHC 4.0...')
            #####
            # For netMHC-4 prediction, only predict on alleles for which data exists
            netmhc_alleles = list(pd.read_table(netmhc4_alleleslist, header=None, usecols=[0])[0])
            alleles_for_prediction = list(set(netmhc_alleles) & set([x.replace(':', '') for x in hla_alleles]))
            logger.info('Only predicting on the following HLA-alleles: ' + ','.join(sorted(set(alleles_for_prediction))))

            run_netmhc_cmd = 'rm -f ' + netmhc_output_pfx + '.txt; '+ \
                                netmhc4_bin + ' -s ' + \
                                ' -a ' + ','.join(sorted(set(alleles_for_prediction))) + \
                                ' -f ' + mutated_sequences_fa + \
                                ' -l ' + ','.join(map(str, peptide_lengths)) + \
                                ' -inptype 0 ' + \
                                ' -xls ' + \
                                ' -xlsfile ' + netmhc_output_pfx + '.xls ' + \
                                ' > ' + netmhc_output_pfx + '.txt'
            execute_cmd(run_netmhc_cmd)

        logger.info('Cleaning up NetMHCpan-4.0 run directory')
        cleanup_cmd = 'rm -rf ' + output_dir + '/scratch/*/*.pred '
        execute_cmd(cleanup_cmd)

    except Exception:
        logger.error('Error while running NetMHCpan and NetMHC')
        logger.error(traceback.format_exc())
        logger.info('Cleaning up NetMHCpan run directory')
        cleanup_cmd = 'rm -rf ' + output_dir + '/scratch/*/*.pred '
        execute_cmd(cleanup_cmd)
        exit(1)

    #######################
    ### Parse NetMHCpan 4.0, netMHC 4.0 output; and, generate output files
    #######################
    try:
        logger.info('')
        logger.info('Parse NetMHC-4.0 and NetMHCpan-4.0 output....')
        ###
        ### Parse NetMHCpan and NetMHC output into a single file with binding scores 
        ###
        parse_output_cmd = "grep -P \"^\\s*\\d+\\s*HLA\\-\"  | sed -r \'s/\\s+/\\t/g\' | sed -r \'s/^\\s*//g\' | cut -f 2-4,10,12-14 | "
        combined_output = sample_path_pfx + '.netmhcpan_netmhc_combined.output.txt'
        generate_output_cmd = 'echo -e \'algorithm\\thla_allele\\tpeptide\\tcore\\ticore\\tscore\\taffinity\\trank\'' + \
                                ' > ' + combined_output + '; ' + \
                                ' cat ' + netmhcpan_output_pfx + '.txt | ' + \
                                parse_output_cmd + \
                                ' awk \'{print \"NetMHCpan-4.0\\t\"$0}\' >> ' + combined_output + '; ' + \
                                ' cat ' + netmhc_output_pfx + '.txt | ' + \
                                parse_output_cmd + \
                                ' awk \'{print \"NetMHC-4.0\\t\"$0}\' >> ' + combined_output + '; '
        execute_cmd(generate_output_cmd)

        #logger.info('memory usage: ' + str((psutil.Process(os.getpid()).memory_info().rss/1e9)))

        # read combined_output file containing all neopeptides that have been evaluated by both prediction algorithms
        logger.info('Reading predictions from the two algorithms and evaluating binders')
        np_df = pd.read_table(combined_output).drop_duplicates()

        ## netMHC-4.0 requires and outputs alleles in a different format; just correct the name
        np_df['hla_allele'] = np_df['hla_allele'].map(lambda a: reformat_hla_allele(a))

        ##
        ## This determination of Strong/Weak binder is based on Swanton's PMID: 30894752. 
        ##
        np_df['binder_class'] = 'non binder' ## keep the casing for 'non' as is; for sorting purpose later
        np_df.loc[(np_df['affinity'] < 500) | (np_df['rank'] < 2), 'binder_class'] = 'Weak Binder'
        np_df.loc[(np_df['affinity'] < 50 ) | (np_df['rank'] < 0.5), 'binder_class'] = 'Strong Binder'

        ##
        ## for each 'peptide', multiple binding predictions will be generated for different HLA alleles and by
        ## different algorithms. 'best_binder_for_icore_group' flag (True/False) represents the best binding 
        ## prediction across the different alleles/algorithms for a given icore. Here, we are sorting by the 
        ## columns below to select the one with best binding prediction and only the top row is retained.
        ##
        np_by_peptide_df = np_df.sort_values(['binder_class', 'rank', 'affinity'], 
                                    ascending=[True, True, True]).groupby('icore').first().reset_index()
        np_by_peptide_df['best_binder_for_icore_group'] = True

        ## annotate np_df with the 'best_binder_for_icore_group'
        np_annotated_df = pd.merge(np_df, np_by_peptide_df, how='left', 
                        on=['algorithm', 'hla_allele', 'peptide', 'core', 
                            'icore', 'score', 'affinity', 'rank', 'binder_class'])

        np_annotated_df['best_binder_for_icore_group'] = np_annotated_df['best_binder_for_icore_group'].fillna(False)
        np_annotated_df.loc[(np_annotated_df['binder_class'] == 'non binder'), 'binder_class'] = 'Non Binder'

        # For each neopeptide, we want to check whether that peptide fragment could be generated by any other WT protein 
        # in the genome. Currently, optimal approach is to generate a string of the entire coding sequence in the genome
        # and search each neopeptide against it
        logger.info('Checking if the icore-peptide can be generated from WT sequence from the entire peptidome...')
        
        ref_aa_str = ''
        for cds in cds_seqs.values():
            if cds[0:3] == 'ATG':
                ref_aa_str += mutation.cds_to_aa(cds) + '|'

        # make a list of all unique peptides
        all_peptides = ({row['icore']:1 for index, row in np_df.iterrows()}).keys()

        # parallelize and search each icore peptide against the reference peptidome. Note: deliberately hard-coded 4 cores for now.
        results = Parallel(n_jobs=4)(delayed(find_in_reference_peptides)(all_peptides, i, 4, ref_aa_str) for i in range(1, 5))

        # construct a dataframe of the peptides that are found in other protein coding genes
        icore_in_reference = pd.DataFrame(data={item:1 for sublist in results for item in sublist}.keys(), columns=['icore'])
        icore_in_reference['is_in_wt_peptidome'] = True
        logger.info('...completed!')

        # annotate the neopeptide dataframe with 'is_in_wt_peptidome' -- that will be written to output file
        np_annotated_df = pd.merge(np_annotated_df, icore_in_reference, how='left', on=['icore'])
        np_annotated_df['is_in_wt_peptidome'] = np_annotated_df['is_in_wt_peptidome'].fillna(False)
        
        # make a neopeptide object of each row. 
        all_neopeptides = [neopeptide(row) for index, row in np_annotated_df.iterrows()]

        # parse the mutations (with neopeptides/binding predictions) and compile output files.
        maf_output = []
        predictions_output = []
        logger.info('Annotating MAF with neopeptides that have the stongest binding affinities')
        for mut in mutations:
            # find all neopeptides that are generated for the given mutation.
            mut.match_with_neopeptides(all_neopeptides)
            maf_output.append(mut.get_maf_row_to_print())
            predictions_output.extend(mut.get_predictions_rows_to_print())
            
        maf_output_df = pd.DataFrame.from_items([(s.name, s) for s in maf_output]).T
        maf_output_df.to_csv(sample_path_pfx + '.neoantigens.maf' , sep='\t', index=False)

        predictions_output_df = pd.DataFrame.from_items([(s.name, s) for s in predictions_output]).T
        predictions_output_df.to_csv(sample_path_pfx + '.all_neoantigen_predictions.txt', sep='\t', index=False)
        
    except Exception:
        logger.error('Error while processing NetMHCpan and NetMHC output')
        logger.error(traceback.format_exc())
        exit(1)
    logger.info('neoantigen-dev pipeline execution completed.\nExiting!')

# helper function to properly re-format hla_allele
def reformat_hla_allele(hla_allele):
    if re.match(r'HLA-\w\d\d\d\d$', hla_allele):
        gene, major, minor = re.match(r'(HLA-\w)(\d\d)(\d\d)$', hla_allele).groups()
        hla_allele = gene + '*' + major + ':' + minor
    return hla_allele

# helper function to split the list of peptides into batches and search against the reference peptidome.
def find_in_reference_peptides(peptides_list, batch_id, n_batches, ref_aa_str):
    batch_size =(len(peptides_list)/n_batches) + 1
    st = batch_size * (batch_id - 1)
    en = min(len(peptides_list), batch_size * batch_id)

    matches = []
    for i in range(st, en):
        if peptides_list[i] in ref_aa_str:
            matches.append(peptides_list[i])
    return matches

def check_file_exists_and_not_empty(filename):
    if os.path.isfile(filename) and os.stat(filename).st_size > 100:
        return True
    return False

def get_timestamp():
    return time.strftime('%Y-%m-%d  %H:%M:%S')

def execute_cmd(cmd):
    logger = logging.getLogger('neoantigen')
    logger.debug('Executing command: ' + cmd)
    output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    logger.debug('Done')

def load_transcript_fasta(fa_file):
    seqs = dict()
    if fa_file[-3:len(fa_file)] == '.gz':
        lines = gzip.open(fa_file, 'rb').readlines()
    else:
        lines = open(fa_file).readlines()
    idx = 0
    while idx < len(lines):
        line = lines[idx]
        m = re.search('^>(ENST\d+)\s', line)
        transcript_id = ''
        if not m:
            sys.exit('Error parsing transcript file ' + fa_file + ' at line: ' + line)
        else:
            transcript_id = m.group(1)

        idx = idx + 1
        seq_str = ''
        while idx < len(lines) and not re.match('^>ENST', lines[idx]):
            seq_str = seq_str + lines[idx].strip()
            idx = idx + 1
            seqs[transcript_id] = seq_str

    return seqs

#
# class to hold the binding prediction for each peptide/hla_allele and algorithm
#
class neopeptide(object):
    row = None
    algorithm = ''
    hla_allele = ''
    peptide = ''
    core = ''
    icore = ''
    score = 0
    binding_affinity = 10000
    rank = 100
    best_binder_for_icore_group = ''
    binder_class = ''
    is_in_wt_peptidome = False

    def __init__(self, row):
        self.row = row
        self.algorithm = row['algorithm']
        self.hla_allele = row['hla_allele']
        self.peptide = row['peptide']
        self.core = row['core']
        self.icore = row['icore']
        self.score = row['score']
        self.binding_affinity = row['affinity']
        self.rank = row['rank']
        self.binder_class = row['binder_class']
        self.best_binder_for_icore_group = row['best_binder_for_icore_group']
        self.is_in_wt_peptidome = row['is_in_wt_peptidome']
              
    def is_strong_binder(self):
        if self.binder_class == 'Strong Binder':
            return True
        return False

    def is_weak_binder(self):
        if self.binder_class == 'Weak Binder':
            return True
        return False

#
# class to hold the list of neopeptides and helper functions to identify strong/weak binders
#
class binding_predictions(object):
    neopeptides = None

    def __init__ (self, neopeptides):
        self.neopeptides = neopeptides

    def add_neopeptide(self, np):
        self.neopeptides.append(np)

    def get_strong_binders(self):
        return [x for x in self.neopeptides if x.is_strong_binder()]

    def get_weak_binders(self):
        return [x for x in self.neopeptides if x.is_weak_binder()]

    def get_all_binders(self):
        return [x for x in self.neopeptides if x.is_strong_binder() or x.is_weak_binder()]

    def get_best_binder(self):
        if (len(self.neopeptides) == 0):
            return None
        return sorted(self.neopeptides, key=lambda x: x.rank, reverse=False)[0]

#
# mutation class holds each row in the maf and has 
#
class mutation(object):
    maf_row = None
    cds_seq = ''
    cdna_seq = ''
    wt_cds = ''
    wt_aa = ''
    wt_altered_aa = ''
    mt_cds = ''
    mt_aa = ''
    mt_altered_aa = ''
    identifier_key = ''
    predicted_neopeptides = None

    def __init__(self, maf_row, cds_seq, cdna_seq):
        self.maf_row = maf_row
        self.cds_seq = cds_seq
        self.cdna_seq = cdna_seq
        self.predicted_neopeptides = binding_predictions([])
        self.identifier_key = self.maf_row['Tumor_Sample_Barcode'] + '_' + \
                              str(self.maf_row['Chromosome']) + '_' + \
                              str(self.maf_row['Start_Position']) + '-' + \
                              str(self.maf_row['End_Position']) + '_' + \
                              self.maf_row['Reference_Allele'] + '_' + \
                              self.maf_row['Tumor_Seq_Allele2']
        
    ### Check if the variant_classification is among those that can generate a neoantigen
    def is_non_syn(self):
        types = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 
                'In_Frame_Ins', 'Missense_Mutation', 'Nonstop_Mutation']

        return self.maf_row['Variant_Classification'] in types and not pd.isnull(self.maf_row['HGVSp_Short'])

    ### helper function #source: stackoverflow.
    def reverse_complement(self, dna):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return ''.join([complement[base] for base in dna[::-1]])

    ### helper function to translate cDNA sequence
    @staticmethod
    def cds_to_aa(cds):
        # https://www.geeksforgeeks.org/dna-protein-python-3/
        codon_table = {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_', 'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
        }
        protein = ''

        for i in range(0, len(cds), 3):
            codon = cds[i:i + 3]
            if len(codon) != 3:
                ## This is unusual; in some cases in Ensembl the CDS length is not a multiple of 3. Eg: ENST00000390464
                ## For this reason, decided not to throw an error and just stop translating if the CDS ends with a non-triplet
                #print 'CDS ends with non-triplet: ' + codon + ' ' + cds
                break
            if codon_table[codon] == '_': # stop codon reached
                break
            protein += codon_table[codon]
        return protein

    # function that parses the HGVSc and constructs the WT and mutated coding sequences for the given mutation.
    def generate_translated_sequences(self, pad_len=10):
        if not self.is_non_syn():
            return None
        
        ## append the 3'UTR to the CDS -- to account for non stop mutations and indels that shift the canonical stop
        if not self.cds_seq in self.cdna_seq:
            print 'Skipping because the CDS is not contained within cDNA. Note: only 2 transcripts/peptides are like this'
            return None

        hgvsc = self.maf_row['HGVSc']
        position, ref_allele, alt_allele, sequence, hgvsc_type = [-1, '', '', '', 'ONP']

        if re.match(r'^c\.(\d+).*([ATCG]+)>([ATCG]+)$', hgvsc):
            position, ref_allele, alt_allele = re.match(r'^c\.(\d+).*(\w+)>(\w+)', hgvsc).groups()

        elif re.match(r'^c\.(\d+).*del([ATCG]+)ins([ATCG]+)$', hgvsc):
            position, ref_allele, alt_allele = re.match(r'^c\.(\d+).*del([ATCG]+)ins([ATCG]+)$', hgvsc).groups()

        elif re.match(r'^c\.(\d+).*(dup|ins|del|inv)([ATCG]+)$', hgvsc):
            position, hgvsc_type, sequence = re.match(r'^c\.(\d+).*(dup|ins|del|inv)([ATCG]+)$', hgvsc).groups()

        else:
            sys.exit('Error: not one of the known HGVSc strings: ' + hgvsc)

        position = int(position) - 1
        if hgvsc_type in 'dup,ins':
            alt_allele = sequence
        elif hgvsc_type == 'del':
            ref_allele = sequence
        elif hgvsc_type == 'inv':
            ref_allele = sequence
            alt_allele = self.reverse_complement(sequence)

        ## start of mutated region in CDS
        cds = re.search(self.cds_seq + '.*', self.cdna_seq).group()

        seq_5p = cds[0:position]
        seq_3p = cds[position:len(cds)]

        #print self.hgvsp + '\t' + self.variant_class + '\t' + self.variant_type + '\t' + self.ref_allele + '\t' + self.alt_allele + \
        #      '\t' + self.cds_position + '\nFull CDS: ' + self.cds_seq + '\nSeq_5: ' + seq_5p + '\nSeq_3' + seq_3p + '\n>mut_1--' + mut_cds_1 + '\n>mut_2--' + mut_cds_2 + '\n>mut_3--' + mut_cds_3
        self.wt_cds = seq_5p + ref_allele + seq_3p[len(ref_allele):len(seq_3p)]
        self.mt_cds = seq_5p + alt_allele + seq_3p[len(ref_allele):len(seq_3p)]
        wt = mutation.cds_to_aa(self.wt_cds)
        mt = mutation.cds_to_aa(self.mt_cds)

        ### identify regions of mutation in WT and MT sequences. 
        ### logic is to match the wt and mt sequences first from the beginning until a mismatch is found; and, then,
        ### start from the end of both sequences until a mismatch is found. the intervening sequence represents the WT and MT sequences
        ### Note, aside from missenses, the interpretation of WT sequence is ambiguous.
        len_from_start = len_from_end = 0

        ## from start
        for i in range(0, min(len(wt), len(mt))):
            len_from_start = i
            if wt[i:i + 1] != mt[i:i + 1]:
                break

        ## from end
        wt_rev = wt[::-1]
        mt_rev = mt[::-1]
        for i in range(0, min(len(wt), len(mt))):
            len_from_end = i
            if len_from_end + len_from_start >= min(len(wt), len(mt)) or \
                wt_rev[i:i + 1] != mt_rev[i:i + 1]:
                break

        wt_start = len_from_start
        wt_end = len(wt) - len_from_end

        mt_start = len_from_start
        mt_end = len(mt) - len_from_end

        self.wt_aa = wt
        self.mt_aa = mt

        self.wt_altered_aa = wt[max(0, wt_start - pad_len + 1):min(len(wt), wt_end + pad_len-1)]
        self.mt_altered_aa = mt[max(0, mt_start - pad_len + 1):min(len(mt), mt_end + pad_len-1)]

    # function to iterate over all the the neopeptide predictions made in the entire MAF and identify
    # which neopeptides are generated by this mutation object
    def match_with_neopeptides(self, all_neopeptides):
        for np in all_neopeptides:
            # make sure the neopeptide is not a peptide fragment of the wild-type protein
            if np.icore in self.mt_altered_aa and np.icore not in self.wt_aa:
                self.predicted_neopeptides.add_neopeptide(copy.deepcopy(np))
    
    # simply prints the original row in the MAF file along with some neoantigen prediction specific 
    # appended at the end
    def get_maf_row_to_print(self):
        row = self.maf_row
        row['neo_maf_identifier_key'] = self.identifier_key

        if self.predicted_neopeptides.get_best_binder() is not None:
            best_binder = self.predicted_neopeptides.get_best_binder()

            strong_binders = self.predicted_neopeptides.get_strong_binders()
            weak_binders = self.predicted_neopeptides.get_weak_binders()
            row['neo_best_icore_peptide'] = best_binder.icore
            row['neo_best_rank'] = best_binder.rank
            row['neo_best_binding_affinity'] = best_binder.binding_affinity
            row['neo_best_binder_class'] = best_binder.binder_class
            row['neo_best_is_in_wt_peptidome'] = best_binder.is_in_wt_peptidome
            row['neo_best_algorithm'] = best_binder.algorithm
            row['neo_best_hla_allele'] = best_binder.hla_allele
            row['neo_n_all_binders'] = len(strong_binders) + len(weak_binders)
            row['neo_n_strong_binders'] = len(strong_binders)
            row['neo_n_weak_binders'] = len(weak_binders)
        else:
            row['neo_best_icore_peptide'] = ''
            row['neo_best_rank'] = ''
            row['neo_best_binding_affinity'] = ''
            row['neo_best_binder_class'] = ''
            row['neo_best_is_in_wt_peptidome'] = ''
            row['neo_best_algorithm'] = ''
            row['neo_best_hla_allele'] = ''
            row['neo_n_all_binders'] = 0
            row['neo_n_strong_binders'] = 0
            row['neo_n_weak_binders'] = 0
        return row

    # simply prints the original row of the 'combined_output' of neoantigen predictions along with additional columns
    def get_predictions_rows_to_print(self):
        rows = []
        for prediction in self.predicted_neopeptides.neopeptides:
            prediction.row['neo_maf_identifier_key'] = self.identifier_key
            rows.append(prediction.row)
        return rows           


if __name__ == '__main__':
    main()
