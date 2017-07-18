#!/opt/common/CentOS_6-dev/python/python-2.7.10/bin/python
import argparse, os, re, errno
import time, sys, datetime, dateutil.parser
import traceback
import pandas as pd
from ConfigParser import ConfigParser
import subprocess

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
        "Wrapper to execute Van-Allen lab's neoantigen pipeline. In its present state, \n"
        "this wrapper is primarily for preliminary analysis.  Many of paths to the scripts,\n "
        "binaries and data files are hard-coded (this will change later on).  "
    )
    prog_epilog = (
        "NOTES:  This pipeline is not optimized for run time efficiency.\n"
        "\n\n"
    )

    parser = argparse.ArgumentParser(description=prog_description,
                                     epilog=prog_epilog,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     add_help=True)
    parser.add_argument("--config_file",
                        required=True,
                        help="REQUIRED. See: /ifs/res/taylorlab/bandlamc/neoantigens/cmo_neoantigen/cmo_neoantigen.config")
    parser.add_argument("--normal_bam",
                        required=False,
                        help="full path to normal bam file. mutually exclusive with --hla_file.  One of the two is required. ")
    parser.add_argument("--maf_file",
                        required=True,
                        help="REQUIRED. expects a CMO maf file (post vcf2maf.pl)")
    parser.add_argument("--output_dir",
                        required=True,
                        help="REQUIRED. output directory")
    parser.add_argument("--hla_file",
                        required=False,
                        help="polysolver output (winners.hla.txt) for the patient. If not provided, \n"
                        "POLYSOLVER is run. Option mutually exclusive with --normal_bamfile. \n"
                        "Either option is required")
    parser.add_argument("--keep_tmp_files",
                        required=False,
                        help="keeps POLYSOLVER's temporary files. for debugging purposes. " +
                              "Note: TMP files can be more than 5GB. Default: true",
                        action='store_true')
    parser.add_argument("--force_rerun_polysolver",
                        required=False,
                        help="ignores any existing polysolver output and re-runs it. Default: false",
                        action='store_true')
    parser.add_argument("--force_rerun_netmhc",
                        required=False,
                        help="ignores any existing netMHCpan output and re-runs it. Default: false",
                        action='store_true')

    args = parser.parse_args()

    normal_bamfile = str(args.normal_bam)
    maf_file = str(args.maf_file)
    output_dir = str(args.output_dir)
    hla_file = str(args.hla_file)
    config_file_path = str(args.config_file)

    keep_tmp_files = False
    force_polysolver = False
    force_netmhc = False
    patient_prefix = "sample"  ## dummy prefix

    if args.keep_tmp_files is not None and args.keep_tmp_files:
        keep_tmp_files = True

    if args.force_rerun_polysolver is not None and args.force_rerun_polysolver is True:
        force_polysolver = True

    if args.force_rerun_polysolver is not None and args.force_rerun_polysolver is True:
        force_netmhc = True

    if args.normal_bam is None and args.hla_file is None:
        print >> sys.stderr, "Error: --normal_bam or --hla_file is required. Exiting."
        exit(1)

    if args.normal_bam is not None and not os.path.exists(normal_bamfile) :
        print >> sys.stderr, "Error: --normal_bam '" + normal_bamfile + "' does not exist. Exiting."
        exit(1)

    if args.hla_file is not None and not os.path.exists(hla_file) :
        print >> sys.stderr, "Error: --hla_file '" + hla_file + "' does not exist. Exiting."
        exit(1)

    if not os.path.exists(maf_file):
        print >> sys.stderr, "Error: --maf_file '" + maf_file + "' does not exist. Exiting."
        exit(1)

    os.system("mkdir -p " + output_dir)

    config = ConfigParser()
    config.read(os.path.dirname(os.path.realpath(__file__)) + "/cmo_neoantigen.config")
    polysolver_bin = config.get('POLYSOLVER', 'polysolver_bin')

    maf2fasta_py = os.path.dirname(os.path.realpath(__file__)) + "/neoantigen_calling_pipeline/mafToFastaV2.py"
    runNetMHC_py = os.path.dirname(os.path.realpath(__file__)) + "/neoantigen_calling_pipeline/runNetMHCpan.py"
    mut_post_process_py = os.path.dirname(os.path.realpath(__file__)) + "/neoantigen_calling_pipeline/mutationPostProcess.py"

    #######################
    ### POLYSOLVER
    #######################
    if args.hla_file is None or not os.path.exists(hla_file):
        try:
            hla_file = output_dir + "/polysolver/winners.hla.txt"
            if not force_polysolver:
                sys.stdout.write( time.strftime("%Y-%m-%d  %H:%M:%S") + ": Checking for previous run of POLYSOLVER..." )
                sys.stdout.flush()

            if os.path.isfile(os.path.abspath(hla_file)) and os.stat(os.path.abspath(hla_file)) > 0 and not force_polysolver:
                sys.stdout.write(" Found! (using hla alleles in " + hla_file + ")\n")
            else:
                if not force_polysolver:
                    sys.stdout.write(" Not Found!\n")
                    print time.strftime("%Y-%m-%d  %H:%M:%S") + ": Running POLYSOLVER..."
                else:
                    print time.strftime("%Y-%m-%d  %H:%M:%S") + ": Running POLYSOLVER... (in --force_rerun_polysolver mode)"

                run_polysolver_cmd = "export NUM_THREADS=8; " + polysolver_bin  + " " + normal_bamfile + \
                                     " Unknown 1 hg19 STDFQ 0 " + output_dir + "/polysolver/" + " > " + \
                                     output_dir + "/polysolver.log " + " 2> " + output_dir + "/polysolver.err "

                print time.strftime("%Y-%m-%d  %H:%M:%S") + ": Executing command: \n\n" + run_polysolver_cmd
                sys.stdout.flush()
                subprocess.check_output(run_polysolver_cmd, shell=True, stderr=subprocess.STDOUT)
                print time.strftime("%Y-%m-%d  %H:%M:%S") + ": Done\n"

                if not keep_tmp_files:
                    print time.strftime("%Y-%m-%d  %H:%M:%S") + ": Cleaning up POLYSOLVER run directory\n"
                    os.system("rm -rf " + output_dir + "/polysolver/*fastq " +
                              output_dir + "/polysolver/*bam " +
                              output_dir + "/polysolver/*bai " )

        except Exception, e:
            print >> sys.stderr, time.strftime("%Y-%m-%d  %H:%M:%S") + \
                                 ": ERROR: could not run POLYSOLVER. " \
                                 "Check polysolver.log and polysolver.err files in " + output_dir
            print >> sys.stderr, traceback.format_exc()
            if not keep_tmp_files:
                print >> sys.stderr, time.strftime("%Y-%m-%d  %H:%M:%S") + ": Cleaning up directory"
                os.system("rm -rf " + output_dir + "/polysolver/*fastq " +
                          output_dir + "/polysolver/*bam " +
                          output_dir + "/polysolver/*bai ")
            exit(1)

    #######################
    ### MAF2FASTA
    #######################
    try:
        print time.strftime("%Y-%m-%d  %H:%M:%S") + ": Generating candidate neoantigen peptide sequences..."

        oncotator_indel_maf_file = output_dir + "/" + patient_prefix + ".oncotator_formatted.INDEL.maf"
        oncotator_snp_maf_file = output_dir + "/" + patient_prefix + ".oncotator_formatted.SNP.maf"

        convert_to_oncotator_format(maf_file, oncotator_indel_maf_file, oncotator_snp_maf_file)

        maf2fasta_snv_cmd   = "export CMO_NEOANTIGEN_CONFIG=" + config_file_path + "; " + \
                              "python " + maf2fasta_py + " " + \
                              oncotator_snp_maf_file   + " 0 9,10 " + patient_prefix + " " + \
                              output_dir
        maf2fasta_indel_cmd = "export CMO_NEOANTIGEN_CONFIG=" + config_file_path + "; " + \
                              "python " + maf2fasta_py + " " + \
                              oncotator_indel_maf_file + " 1 9,10 " + \
                              patient_prefix + " " + output_dir

        #cleanup old peptide files (that is because mafToFastaV2.py attempts to "append" to these files
        os.system("rm -f " + output_dir + "/len9pep*txt " + output_dir + "/len10pep*txt")

        print time.strftime("%Y-%m-%d  %H:%M:%S") + ": Executing command: \n\n" + maf2fasta_snv_cmd
        sys.stdout.flush()
        subprocess.check_output(maf2fasta_snv_cmd, shell=True, stderr=subprocess.STDOUT)
        print time.strftime("%Y-%m-%d  %H:%M:%S") + ": Done \n"

        print time.strftime("%Y-%m-%d  %H:%M:%S") + ": Executing command: \n\n" + maf2fasta_indel_cmd
        sys.stdout.flush()
        subprocess.check_output(maf2fasta_indel_cmd, shell=True, stderr=subprocess.STDOUT)
        print time.strftime("%Y-%m-%d  %H:%M:%S") + ": Done \n"

    except Exception, e:
        print >> sys.stderr, time.strftime("%Y-%m-%d  %H:%M:%S") + ": Error: could not run mafToFastaV2.py."
        print >> sys.stderr, time.strftime("%Y-%m-%d  %H:%M:%S") + ": Error message: "
        print >> sys.stderr, traceback.format_exc()
        exit(1)

    #######################
    ### NetMHCpan
    #######################
    try:
        sys.stdout.write(time.strftime("%Y-%m-%d  %H:%M:%S") + ": Checking for previous run of NetMHCpan...")
        sys.stdout.flush()
        if (
                not force_netmhc and
                os.path.isfile(output_dir + "/NETMHCpan_out_9SNV.xls")    and os.stat(output_dir + "/NETMHCpan_out_9SNV.xls").st_size > 0 and
                os.path.isfile(output_dir + "/NETMHCpan_out_9InDel.xls")  and os.stat(output_dir + "/NETMHCpan_out_9InDel.xls").st_size > 0 and
                os.path.isfile(output_dir + "/NETMHCpan_out_10SNV.xls")   and
                os.path.isfile(output_dir + "/NETMHCpan_out_10InDel.xls")
            ):  ## this is an imperfect check for completeness;
            sys.stdout.write(" found\n")
            sys.stdout.write("\nWARNING!!! Using output from previous run of NetMHCpan. Please re-run this with --force_rerun_netmhc "
                             "if you previous NetMHC run failed.\n\n")
        else:
            print time.strftime("%Y-%m-%d  %H:%M:%S") + ": Starting NETMHCpan..."

            runNetMHCpan_cmd = "export CMO_NEOANTIGEN_CONFIG=" + config_file_path + "; " + \
                               "python " + runNetMHC_py + " " + \
                               output_dir + "/len9pep_FASTA_indel.txt,"  + \
                               output_dir + "/len9pep_FASTA_snv.txt," + \
                               output_dir + "/len10pep_FASTA_indel.txt," + \
                               output_dir + "/len10pep_FASTA_snv.txt" + \
                               " " + hla_file + " 9,9,10,10 1 " + output_dir

            print time.strftime("%Y-%m-%d  %H:%M:%S") + ": Executing command: \n\n" + runNetMHCpan_cmd
            sys.stdout.flush()
            subprocess.check_output(runNetMHCpan_cmd, shell=True, stderr=subprocess.STDOUT)
            print time.strftime("%Y-%m-%d  %H:%M:%S") + ": Done\n"

    except Exception, e:
        print >> sys.stderr, time.strftime("%Y-%m-%d  %H:%M:%S") + ": Error: could not run NetMHCpan"
        print >> sys.stderr, time.strftime("%Y-%m-%d  %H:%M:%S") + ": Error message: "
        print >> sys.stderr, traceback.format_exc()
        exit(1)

    #######################
    ### Process NetMHCpan output
    #######################
    try:
        print time.strftime("%Y-%m-%d  %H:%M:%S") + ": Post-processing..."
        mutationPostProcess_cmd = "export CMO_NEOANTIGEN_CONFIG=" + config_file_path + "; python " + \
                                  mut_post_process_py + " " + \
                                  output_dir + "/NETMHCpan_out_9SNV.xls," + \
                                  output_dir + "/NETMHCpan_out_9InDel.xls," + \
                                  output_dir + "/NETMHCpan_out_10SNV.xls," + \
                                  output_dir + "/NETMHCpan_out_10InDel.xls " + \
                                  output_dir + "/len9pep_headermap_snv.txt," + \
                                  output_dir + "/len9pep_headermap_indel.txt," + \
                                  output_dir + "/len10pep_headermap_snv.txt," + \
                                  output_dir + "/len10pep_headermap_indel.txt " + \
                                  " 9,9,10,10 " + patient_prefix + " 1 " + output_dir

        print time.strftime("%Y-%m-%d  %H:%M:%S") + ": Executing command: \n\n" + mutationPostProcess_cmd
        subprocess.check_output(mutationPostProcess_cmd, shell=True, stderr=subprocess.STDOUT)
        print time.strftime("%Y-%m-%d  %H:%M:%S") + ": Done \n"

    except Exception, e:
        print >> sys.stderr, time.strftime("%Y-%m-%d  %H:%M:%S") + ": Error: could not run post-processing"
        print >> sys.stderr, time.strftime("%Y-%m-%d  %H:%M:%S") + ": Error message: "
        print >> sys.stderr, traceback.format_exc()
        exit(1)

    #######################
    ### Annotate input maf with neoantigen binding affinities
    #######################
    print time.strftime("%Y-%m-%d  %H:%M:%S") + ": Starting annotation of input maf with binding affinities"

    processed_combined_netmhcpan_outputfile = output_dir + "/" + patient_prefix + "_processedcombinedNETMHCpan_out.txt"
    netmhcpan_results_file = open(processed_combined_netmhcpan_outputfile, "r")
    for line in netmhcpan_results_file.readlines():
        columns = line.strip().split("\t")
        if columns[0] == "patient":
            continue

    netmhc_out_df = pd.read_table(processed_combined_netmhcpan_outputfile )

    maf_df = pd.read_table(maf_file)
    maf_df['key'] = "chr" + maf_df['Chromosome'].map(str) + ":" + maf_df['Start_Position'].map(str) + "-" + \
                    maf_df['End_Position'].map(str) + "," + maf_df['Hugo_Symbol'].map(str) + "," + \
                    maf_df['HGVSp_Short'].map(str)
    netmhc_out_df['key'] = netmhc_out_df['chrom_loc'].map(str) + "," + netmhc_out_df['gene'].map(str) + "," + \
                           netmhc_out_df['prot_change'].map(str)

    netmhc_out_df = netmhc_out_df[['key', 'HLA', 'pep_length', 'pep_mut', 'aff_mut', 'rank_mut',
                                   'pep_wt', 'aff_wt', 'rank_wt']]

    ### for mutations with >1 nominated neoantigen peptides, select the one with the highest binding affinity
    print time.strftime("%Y-%m-%d  %H:%M:%S") + \
          ": For mutations with multiple neoantigens, selecting the MUT peptide with highest affinity ..."
    netmhc_out_summary = netmhc_out_df['aff_mut'].groupby(netmhc_out_df['key']).apply(get_stats).unstack()
    netmhc_out_summary['key'] = netmhc_out_summary.axes[0].tolist()
    netmhc_out_df = pd.merge(netmhc_out_summary, netmhc_out_df, how='left', on=['key', 'aff_mut'])
    netmhc_out_df = netmhc_out_df.rename(columns={'count': 'total_num_neoantigens'})
    netmhc_out_df = netmhc_out_df[['key', 'HLA', 'pep_length', 'pep_mut', 'aff_mut', 'rank_mut',
                                   'pep_wt', 'aff_wt', 'rank_wt', 'total_num_neoantigens']]

    maf_df_netmhc = pd.merge(maf_df, netmhc_out_df, how='outer', on='key', indicator=True)

    unmatched_aff_file = output_dir + "/" + patient_prefix + "_unmached_netMHCpan_affinities.txt"
    maf_df_netmhc.to_csv(unmatched_aff_file, sep="\t", index=False)

    # QC check to make sure there are no binding affinities in the netMHCpan output that cannot be paired
    # with the mutations in the MAF file.  This shouldn't happen.  Just to make sure that the "key" we are
    # joining with, above, is correct.
    if len(maf_df_netmhc[maf_df_netmhc._merge == 'right_only']) > 0:
        unmatched_aff_file = output_dir + "/" + patient_prefix + "_unmached_netMHCpan_affinities.txt"
        print time.strftime("%Y-%m-%d  %H:%M:%S") + "WARNING: There are binding affinities in the netHMCpan output that " \
              "cannot be properly paired with the mutations in the maf file. " \
              "See file: " + unmatched_aff_file
        maf_df_netmhc[maf_df_netmhc._merge == 'right_only'].to_csv(unmatched_aff_file, sep="\t", index=False)
    del maf_df_netmhc['_merge']
    maf_df_netmhc.to_csv(output_dir + "/" + os.path.basename(maf_file) + ".netMHCpan.neoantigens.maf",  sep="\t", index=False)

    print time.strftime("%Y-%m-%d  %H:%M:%S") + ": Done. Exiting!"


def get_stats(group):
    return{'aff_mut': group.min(), 'count': group.count()}

def convert_to_oncotator_format(inp_maf_file, oncotat_indel_maf_file, oncotat_snp_maf_file):
    inp_maf_file_obj = open(inp_maf_file, 'r')
    oncotat_indel_maf_file_obj = open(oncotat_indel_maf_file, "w", 0)
    oncotat_snp_maf_file_obj = open(oncotat_snp_maf_file, "w", 0)

    for line in inp_maf_file_obj.readlines():
        if re.match(r'^#', line):
            continue

        codon_change = ""
        columns = line.split("\t")

        pat = re.compile('Missense|Nonsense|Nonstop|Frame_Shift|In_Frame')
        if not pat.match(columns[8]):
            continue

        cds_position = columns[52]
        variant_type = columns[9]

        columns[39] = columns[34].replace("dup", "ins")
        columns[41] = columns[36]
        columns[35] = columns[37]
        if re.match(r'^Hugo_Symbol', line):
            columns[40] = "Codon_Change"
            oncotat_indel_maf_file_obj.write( "\t".join(columns))
            oncotat_snp_maf_file_obj.write( "\t".join(columns))
        else:
            acds_start = 0
            acds_end = 0
            if variant_type == "SNP":
                position = int(cds_position.split("/")[0])

                if position % 3 == 0:
                    acds_start = position - 2
                    acds_end = position
                elif position % 3 == 1:
                    acds_start = position
                    acds_end = position + 2
                else:
                    acds_start = position - 1
                    acds_end = position + 1
                columns[55] = columns[55].replace("/", ">")
                codon_change = "c.(" + str(acds_start) + "-" + str(acds_end) + ")" + columns[55]
                columns[40] = codon_change
                oncotat_snp_maf_file_obj.write( "\t".join(columns))

            else:
                if cds_position.find('?') != -1:
                    continue  # ignoring cases where the translational consequence is ambiguous
                              # (eg: start of INDEL in exon and end in UTR/intron).

                positions = (cds_position.split("/")[0]).split("-")
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

                columns[55] = columns[55].replace("/", ">")
                codon_change = "c.(" + str(acds_start) + "-" + str(acds_end) + ")" + columns[55]
                columns[40] = codon_change
                oncotat_indel_maf_file_obj.write( "\t".join(columns))

if __name__ == '__main__':
    main()
