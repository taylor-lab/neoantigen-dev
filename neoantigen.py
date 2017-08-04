#!/opt/common/CentOS_6-dev/python/python-2.7.10/bin/python
import argparse, os, re, errno
import time, sys, datetime, dateutil.parser
import traceback
import pandas as pd
from ConfigParser import ConfigParser
import subprocess
import logging, logging.handlers
import requests
import commands

#####
#
#  cmo_neoantigen - This script is a wrapper to the Neoantigen pipeline written by
#  Claire Margolis in Van-Allen lab (https://github.com/vanallenlab/neoantigen_calling_pipeline).
#  Some paths are hard-coded in the Van-Allen code base. So, we forked a recent branch (05/15/2017)
#  and tweaked it to build a wrapper around it.  Any updates to the original branch should
#  be reviewed and integrated manually into ours'.  This will likely be not an issue in
#  the future as their pipeline matures.
#
#####

def main():
    prog_description = (
        'Wrapper to execute Van-Allen lab\'s neoantigen pipeline. In its present state, \n'
        'this wrapper is primarily for preliminary analysis.  '
    )
    prog_epilog = (
        'NOTES:  This pipeline is not optimized for run time efficiency.\n'
        '\n\n'
    )

    parser = argparse.ArgumentParser(description=prog_description,
                                     epilog=prog_epilog,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     add_help=True)
    parser.add_argument('--config_file',
                        required=True,
                        help='REQUIRED. See: neoantigen-luna.config in the repo')
    parser.add_argument('--sample_id',
                        required=True,
                        help='REQUIRED. sample_id used to limit neoantigen prediction to identify mutations '
                             'associated with the patient in the MAF (column 16). ')
    parser.add_argument('--output_dir',
                        required=True,
                        help='REQUIRED. output directory')
    parser.add_argument('--normal_bam',
                        required=False,
                        help='full path to normal bam file. Either --normal_bam or --hla_file are required.')
    parser.add_argument('--hla_file',
                        required=False,
                        help='POLYSOLVER output file (winners.hla.txt) for the sample. If not provided,'
                             'POLYSOLVER is run. Either --normal_bam or --hla_file are required.')
    parser.add_argument('--maf_file',
                        required=True,
                        help='REQUIRED. expects a CMO maf file (post vcf2maf.pl)')
    parser.add_argument('--keep_tmp_files',
                        required=False,
                        help='keeps POLYSOLVER\'s temporary files. for debugging purposes. ' +
                             'Note: TMP files can be more than 5GB. Default: true',
                        action='store_true')
    parser.add_argument('--force_rerun_polysolver',
                        required=False,
                        help='ignores any existing polysolver output and re-runs it. Default: false',
                        action='store_true')
    parser.add_argument('--force_rerun_netmhc',
                        required=False,
                        help='ignores any existing netMHCpan output and re-runs it. Default: false',
                        action='store_true')

    args = parser.parse_args()

    normal_bamfile = str(args.normal_bam)
    maf_file = str(args.maf_file)
    output_dir = str(args.output_dir)
    hla_file = str(args.hla_file)
    sample_id = str(args.sample_id)

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

    if args.force_rerun_polysolver is not None and args.force_rerun_polysolver is True:
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

    polysolver_bin = config.get('POLYSOLVER', 'polysolver_bin')
    maf2fasta_py = os.path.dirname(os.path.realpath(__file__)) + '/neoantigen_calling_pipeline/mafToFastaV2.py'
    runNetMHC_py = os.path.dirname(os.path.realpath(__file__)) + '/neoantigen_calling_pipeline/runNetMHCpan.py'
    mut_post_process_py = os.path.dirname(
        os.path.realpath(__file__)) + '/neoantigen_calling_pipeline/mutationPostProcess.py'

    #######################
    ### POLYSOLVER
    #######################
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

    #######################
    ### MAF2FASTA
    #######################
    sample_maf_file = output_dir + '/' + sample_id + '.maf'

    try:
        logger.info('Generating candidate neoantigen peptide sequences...')

        extract_sample_mutations_from_maf(maf_file, sample_id, sample_maf_file)

        oncotator_indel_maf_file = output_dir + '/sample.oncotator_formatted.INDEL.maf'
        oncotator_snp_maf_file = output_dir + '/sample.oncotator_formatted.SNP.maf'

        #
        # Convert MAF to oncotator format
        #
        maf_df = pd.read_table(sample_maf_file,
                               comment='#',
                               low_memory=False,
                               header=0)

        maf_df = maf_df[maf_df['Variant_Classification'].isin(['Missense_Mutation',
                                                               'Nonsense_Mutation',
                                                               'Nonstop_Mutation',
                                                               'Frame_Shift_Del',
                                                               'Frame_Shift_Ins',
                                                               'In_Frame_Del',
                                                               'In_Frame_Ins'])]

        expected_num_indels = len( maf_df[maf_df['Variant_Type'].isin(['INS', 'DEL'])] )
        expected_num_snps   = len( maf_df[maf_df['Variant_Type'] == 'SNP'] )

        wall_time = "0:59"
        if expected_num_indels + expected_num_snps > 10000:
            wall_time = "3:00"

        oncotator_file_pfx = output_dir + "/sample.oncotator_formatted"
        convert_script = os.path.dirname(os.path.realpath(__file__)) + "/convertToOncotator.py "
        bsub_cmd = "bsub -n 1 -K -R \"select[internet && mem=4]\" -We " + wall_time + \
                   " -o " + oncotator_file_pfx + ".bsub.output -e " + oncotator_file_pfx + ".bsub.err " + \
                   " python " + convert_script + \
                   " --input_maf " + sample_maf_file + \
                   " --output_maf_prefix " + oncotator_file_pfx

        execute_cmd(bsub_cmd)

        # check if the conversion succeeded.
        converter_log = open(output_dir + '/sample.oncotator_formatted.converter.log', 'r')
        lines = converter_log.readlines()
        if lines[ len(lines) - 1 ].find('Success') == -1:
            logger.error('Error occurred while converting input MAF to oncotator format:\n')
            logger.error('\n' + bsub_cmd)
            exit(1)

        maf2fasta_snv_cmd = 'export CMO_NEOANTIGEN_CONFIG=' + config_file_path + '; ' + \
                            'python ' + maf2fasta_py + ' ' + \
                            oncotator_snp_maf_file + ' 0 9,10 sample ' + output_dir
        maf2fasta_indel_cmd = 'export CMO_NEOANTIGEN_CONFIG=' + config_file_path + '; ' + \
                              'python ' + maf2fasta_py + ' ' + \
                              oncotator_indel_maf_file + ' 1 9,10 sample ' + output_dir

        # cleanup old peptide files (that is because mafToFastaV2.py attempts to 'append' to these files
        execute_cmd('rm -f ' + output_dir + '/len9pep*txt ' + output_dir + '/len10pep*txt')

        logger.info('Running mafToFasta for SNVs')
        execute_cmd(maf2fasta_snv_cmd)

        logger.info('Running mafToFasta for INDELs')
        execute_cmd(maf2fasta_indel_cmd)

    except Exception:
        logger.error('Error while running mafToFasta')
        logger.error(traceback.format_exc())
        logger.error('Exiting.')
        exit(1)

    #######################
    ### NetMHCpan
    #######################
    try:
        logger.info('Running MHC--peptide binding predictions using NetMHCpan...')
        if not force_netmhc and check_for_netmhc_completion(output_dir):
            logger.info('WARNING: Using output from previous run of NetHMCpan. '
                        'If you suspect the integrity of previous run, please re-run this with --force_rerun_mhc')
        else:
            logger.info('Starting NETMHCpan...')
            runNetMHCpan_cmd = 'export CMO_NEOANTIGEN_CONFIG=' + config_file_path + '; ' + \
                               'python ' + runNetMHC_py + ' ' + \
                               output_dir + '/len9pep_FASTA_indel.txt,' + \
                               output_dir + '/len9pep_FASTA_snv.txt,' + \
                               output_dir + '/len10pep_FASTA_indel.txt,' + \
                               output_dir + '/len10pep_FASTA_snv.txt' + \
                               ' ' + hla_file + ' 9,9,10,10 1 ' + output_dir

            execute_cmd(runNetMHCpan_cmd)
            logger.info('Cleaning up NetMHCpan run directory')
            cleanup_cmd = 'rm -rf ' + output_dir + '/scratch/*/*.pred '
            execute_cmd(cleanup_cmd)

    except Exception:
        logger.error('Error while running NetMHCpan')
        logger.error(traceback.format_exc())
        logger.info('Cleaning up NetMHCpan run directory')
        cleanup_cmd = 'rm -rf ' + output_dir + '/scratch/*/*.pred '
        execute_cmd(cleanup_cmd)
        exit(1)


    #######################
    ### Process NetMHCpan output
    #######################
    try:
        logger.info('Post-processing...')
        mutationPostProcess_cmd = 'export CMO_NEOANTIGEN_CONFIG=' + config_file_path + '; python ' + \
                                  mut_post_process_py + ' ' + \
                                  output_dir + '/NETMHCpan_out_9SNV.xls,' + \
                                  output_dir + '/NETMHCpan_out_9InDel.xls,' + \
                                  output_dir + '/NETMHCpan_out_10SNV.xls,' + \
                                  output_dir + '/NETMHCpan_out_10InDel.xls ' + \
                                  output_dir + '/len9pep_headermap_snv.txt,' + \
                                  output_dir + '/len9pep_headermap_indel.txt,' + \
                                  output_dir + '/len10pep_headermap_snv.txt,' + \
                                  output_dir + '/len10pep_headermap_indel.txt ' + \
                                  ' 9,9,10,10 sample 1 ' + output_dir
        execute_cmd(mutationPostProcess_cmd)

    except Exception:
        logger.error('Error: could not run post-processing')
        logger.error(traceback.format_exc())
        exit(1)

    #######################
    ### Annotate input maf with neoantigen binding affinities
    #######################
    logger.info('Starting annotation of input maf with binding affinities')

    netmhc_out_df = pd.read_table(output_dir + '/sample_processedcombinedNETMHCpan_out.txt',
                                  comment='#',
                                  low_memory=False,
                                  header=0)

    maf_df = pd.read_table(sample_maf_file,
                           comment='#',
                           low_memory=False,
                           header=0)

    ## construct a unique key in both dataframes to do a join
    maf_df['key'] = 'chr' + maf_df['Chromosome'].map(str) + ':' + maf_df['Start_Position'].map(str) + '-' + \
                    maf_df['End_Position'].map(str) + ',' + maf_df['Hugo_Symbol'].map(str) + ',' + \
                    maf_df['HGVSp_Short'].map(str)

    netmhc_out_df['key'] = netmhc_out_df['chrom_loc'].map(str) + ',' + netmhc_out_df['gene'].map(str) + ',' + \
                           netmhc_out_df['prot_change'].map(str)

    ## select specific columns
    netmhc_out_df = netmhc_out_df[['key', 'HLA', 'pep_length', 'pep_mut', 'aff_mut', 'rank_mut',
                                   'pep_wt', 'aff_wt', 'rank_wt']]

    #
    # for mutations with >1 nominated neoantigen peptides, select the one with the highest binding affinity
    #
    logger.info('For mutations with multiple neoantigens, selecting the MUT peptide with highest affinity ...')

    netmhc_out_summary = netmhc_out_df['aff_mut'].groupby(netmhc_out_df['key']).apply(get_stats).unstack()

    netmhc_out_summary['key'] = netmhc_out_summary.axes[0].tolist()

    netmhc_out_df = pd.merge(netmhc_out_summary,
                             netmhc_out_df,
                             how='left',
                             on=['key', 'aff_mut'])

    netmhc_out_df = netmhc_out_df.rename(columns={'count': 'total_num_neoantigens'})

    netmhc_out_df = netmhc_out_df[['key',
                                   'HLA',
                                   'pep_length',
                                   'pep_mut',
                                   'aff_mut',
                                   'rank_mut',
                                   'pep_wt',
                                   'aff_wt',
                                   'rank_wt',
                                   'total_num_neoantigens']]

    maf_df_netmhc = pd.merge(maf_df, netmhc_out_df, how='outer', on='key', indicator=True)

    unmatched_aff_file = output_dir + '/' + sample_id + '_unmatched_netMHCpan_affinities.txt'

    maf_df_netmhc.to_csv(unmatched_aff_file, sep='\t', index=False)

    # QC check to make sure there are no binding affinities in the netMHCpan output that cannot be paired
    # with the mutations in the MAF file.  This shouldn't happen.  Just to make sure that the 'key' we are
    # joining with, above, is correct.
    if len(maf_df_netmhc[maf_df_netmhc._merge == 'right_only']) > 0:
        logger.info('WARNING: There are binding affinities in the netHMCpan output that ' \
                    'cannot be properly paired with the mutations in the maf file. ' \
                    'See file: ' + unmatched_aff_file)
        maf_df_netmhc[maf_df_netmhc._merge == 'right_only'].to_csv(unmatched_aff_file,
                                                                   sep='\t',
                                                                   index=False)
    del maf_df_netmhc['_merge']

    maf_df_netmhc.to_csv(output_dir + '/' + sample_id + '.netMHCpan.neoantigens.maf',
                         sep='\t',
                         index=False)

    logger.info('neoantigen-dev pipeline execution completed.\nExiting!')

def check_for_netmhc_completion(output_dir):

    if not os.path.isfile(output_dir + '/NETMHCpan_out_9SNV.xls') or \
            not os.path.isfile(output_dir + '/NETMHCpan_out_9InDel.xls') or \
            not os.path.isfile(output_dir + '/NETMHCpan_out_10SNV.xls') or \
            not os.path.isfile(output_dir + '/NETMHCpan_out_10InDel.xls'):
        return False

    maf2fasta_9snv    = commands.getoutput('cat ' + output_dir +
                                         '/len9pep_FASTA_snv.txt    | grep seq_ | grep _mut | sort | uniq | wc -l')
    maf2fasta_10snv   = commands.getoutput('cat ' + output_dir +
                                         '/len10pep_FASTA_snv.txt   | grep seq_ | grep _mut | sort | uniq | wc -l')
    maf2fasta_9indel  = commands.getoutput('cat ' + output_dir +
                                         '/len9pep_FASTA_indel.txt  | grep seq_ | grep _mut | sort | uniq | wc -l')
    maf2fasta_10indel = commands.getoutput('cat ' + output_dir +
                                         '/len10pep_FASTA_indel.txt | grep seq_ | grep _mut | sort | uniq | wc -l')

    netmhc_9snv     = commands.getoutput('cat ' + output_dir +
                                         '/netMHCpanoutlen_9SNV.txt    | grep seq_ | grep _mut | cut -f 1 -d\'.\' | sort | uniq | wc -l')
    netmhc_10snv    = commands.getoutput('cat ' + output_dir +
                                         '/netMHCpanoutlen_10SNV.txt   | grep seq_ | grep _mut | cut -f 1 -d\'.\' | sort | uniq | wc -l')
    netmhc_9indel   = commands.getoutput('cat ' + output_dir +
                                         '/netMHCpanoutlen_9InDel.txt  | grep seq_ | grep _mut | cut -f 1 -d\'.\' | sort | uniq | wc -l')
    netmhc_10indel  = commands.getoutput('cat ' + output_dir +
                                         '/netMHCpanoutlen_10InDel.txt | grep seq_ | grep _mut | cut -f 1 -d\'.\' | sort | uniq | wc -l')

    incomplete = 0
    if maf2fasta_9snv != netmhc_9snv or \
                    maf2fasta_9indel != netmhc_9indel or \
                    maf2fasta_10snv != netmhc_10snv or \
                    maf2fasta_10indel != netmhc_10indel:
        incomplete = 1

    if incomplete == 1:
        #os.system('rm -f ' + output_dir + '/NETMHCpan_out*txt ' + output_dir + '/netMHCpanoutlen*txt ' )
        return False
    else:
        return True

def get_stats(group):
    return {'aff_mut': group.min(), 'count': group.count()}

def get_timestamp():
    return time.strftime('%Y-%m-%d  %H:%M:%S')

def execute_cmd(cmd):
    logger = logging.getLogger('neoantigen')
    logger.debug('Executing command: ' + cmd)
    subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    logger.debug('Done')

def extract_sample_mutations_from_maf(maf_file, sample_id, sample_maf_file):
    logger = logging.getLogger('neoantigen')
    logger.info('Parsing MAF file for mutations associated with sample: ' + sample_id)

    maf_df = pd.read_table(maf_file,
                         comment='#',
                         low_memory=False,
                         header=0)

    if maf_df.columns[15] != 'Tumor_Sample_Barcode':
        logger.error('Error: \'Tumor_Sample_Barcode\' not found. Use cmo_maf2maf before running this pipeline.')
        exit(1)

    maf_df = maf_df[maf_df.Tumor_Sample_Barcode == sample_id]
    logger.info('Found ' + str(len(maf_df)) + ' mutations')
    maf_df.to_csv(sample_maf_file, sep='\t', index=False)

#
# Function DEPRACATED.  Oncotator's maf annotation is different from every other annotator's output including vcf2maf.
# This function shuffles the columns to make them compatible with Van Allen group's pipeline
#
def convert_to_oncotator_format(inp_maf_file, oncotat_indel_maf_file, oncotat_snp_maf_file):
    logger = logging.getLogger('neoantigen')
    logger.info('Converting input maf to Oncotator format')
    inp_maf_file_obj = open(inp_maf_file, 'r')
    oncotat_indel_maf_file_obj = open(oncotat_indel_maf_file, 'w', 0)
    oncotat_snp_maf_file_obj = open(oncotat_snp_maf_file, 'w', 0)

    for line in inp_maf_file_obj.readlines():
        if re.match(r'^#', line):
            continue

        codon_change = ''
        columns = line.split('\t')

        cds_position = columns[52]  # col 53: "CDS_position"
        variant_type = columns[9]   # col 10: "Variant_Type"

        columns[39] = columns[34].replace('dup', 'ins')  # col 35: "HGVSc"
        columns[41] = columns[36]   # col 37: "HGVSp_Short"
        columns[35] = columns[37]   # col 38: "Transcript_ID"


        if re.match(r'^Hugo_Symbol', line):
            columns[40] = 'Codon_Change'      # Codon_Change in Oncotator is 41st column
            oncotat_indel_maf_file_obj.write('\t'.join(columns))
            oncotat_snp_maf_file_obj.write('\t'.join(columns))
        else:
            acds_start = 0
            acds_end = 0
            pat = re.compile('Missense|Nonsense|Nonstop|Frame_Shift|In_Frame')
            if not pat.match(columns[8]):
                continue

            if variant_type == 'SNP':
                position = int(cds_position.split('/')[0])

                if position % 3 == 0:
                    acds_start = position - 2
                    acds_end = position
                elif position % 3 == 1:
                    acds_start = position
                    acds_end = position + 2
                else:
                    acds_start = position - 1
                    acds_end = position + 1

                columns[55] = columns[55].replace('/', '>')     # col 56: "Codons"
                codon_change = 'c.(' + str(acds_start) + '-' + str(acds_end) + ')' + columns[55]
                columns[40] = codon_change
                oncotat_snp_maf_file_obj.write('\t'.join(columns))

            else:
                if cds_position.find('?') != -1:
                    continue  # ignoring cases where the translational consequence is ambiguous
                              # (eg: start of INDEL in exon and end in UTR/intron).

                positions = (cds_position.split('/')[0]).split('-')
                position_s = int(positions[0])
                if len(positions) > 1:
                    position_e = int(positions[1])
                else:
                    position_e = position_s

                if position_s % 3 == 0:
                    acds_start = position_s - 2
                elif position_s % 3 == 1:
                    acds_start = position_s
                else:
                    acds_start = position_s - 1

                if position_e % 3 == 0:
                    acds_end = position_e
                elif position_e % 3 == 1:
                    acds_end = position_e + 2
                else:
                    acds_end = position_e + 1

                columns[55] = columns[55].replace('/', '>')
                codon_change = 'c.(' + str(acds_start) + '-' + str(acds_end) + ')' + columns[55]
                columns[40] = codon_change
                oncotat_indel_maf_file_obj.write('\t'.join(columns))


if __name__ == '__main__':
    main()
