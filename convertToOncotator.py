#!/opt/common/CentOS_6-dev/python/python-2.7.10/bin/python
import argparse, os, re, errno
import time, sys, datetime, dateutil.parser
import traceback
import pandas as pd
from ConfigParser import ConfigParser
import subprocess
import logging, logging.handlers
import requests

#####
#
# Takes an MAF and annotates HGVSc, codon_change columns with Oncotator. To be used
# only by the neoantigen-dev pipeline.
#
#####

def main():
    prog_description = (
        'Convert MAF to oncotator format.  Only for neoantigen-dev.  Only coding \n'
        'variants (Missense|Nonsense|Nonstop|Frame_Shift|In_Frame) are written to output maf.  '
    )

    parser = argparse.ArgumentParser(description=prog_description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     add_help=True)
    parser.add_argument('--input_maf',
                        required=True,
                        help='REQUIRED. Input MAF that is formatted with maf2maf')
    parser.add_argument('--output_maf_prefix',
                        required=True,
                        help='Prefix of SNP and INDEL mafs annoated by Oncotator. '
                             '<output_maf_prefix>.INDEL.maf and <output_maf_prefix>.SNP.maf are output')

    args = parser.parse_args()

    input_maf = str(args.input_maf)
    output_maf_prefix = str(args.output_maf_prefix)

    log_out = open(output_maf_prefix + '.converter.log', 'w')
    try:
        maf_df = pd.read_table(input_maf,
                      comment='#',
                      low_memory=False,
                      header=0,
                      dtype=str,
                      index_col=False)

        maf_df = maf_df[maf_df['Variant_Classification'].isin(['Missense_Mutation',
                                                               'Nonsense_Mutation',
                                                               'Nonstop_Mutation',
                                                               'Frame_Shift_Del',
                                                               'Frame_Shift_Ins',
                                                               'In_Frame_Del',
                                                               'In_Frame_Ins'])]

        # reindex the dataframe
        maf_df.index = range(len(maf_df))

        # first rename columns in old MAF to avoid collisions
        maf_df.rename(columns={
            'HGVSc': 'HGVSc_',
            'HGVSp_Short' : 'HGVSp_Short_',
            'Transcript_ID' : 'Transcript_ID_'},
            inplace=True)

        # update column names of the new columns
        maf_df.rename(columns={
            'HGVSp': 'Transcript_ID',
            't_depth': 'HGVSc',
            't_ref_count': 'codon_change',
            't_alt_count': 'HGVSp_Short'},
            inplace=True)

        rows_to_rem = []
        for index, row in maf_df.iterrows():
            if row['Chromosome'] == "MT":
                rows_to_rem.append(index)
                continue

            tag = str(row['Chromosome']) + '_' + \
                  str(row['Start_Position']) + '_' + \
                  str(row['End_Position']) + '_' + \
                  row['Reference_Allele'] + '_' + \
                  row['Tumor_Seq_Allele2']

            response_json = requests.post("http://www.broadinstitute.org/oncotator/mutation/" + tag + "/").json()

            #print str(index) + " > " + row['Chromosome'] + "_" + row['Start_Position'] + '_' + row['End_Position'] + '_' + row['Reference_Allele'] + '_' + row['Tumor_Seq_Allele2']
            #log_out.write(row['Chromosome'] + "_" + row['Start_Position'] + '_' + row['End_Position'] + '_' + row['Reference_Allele'] + '_' + row['Tumor_Seq_Allele2'] + '\n')
            if response_json['transcript_change'] == "" or \
                            response_json['transcript_id'] == "" or \
                            response_json['codon_change'] == "":
                rows_to_rem.append(index)
                continue

            hgvsc = response_json['transcript_change'].replace('dup', 'ins')
            maf_df.set_value(index, 'Transcript_ID', response_json['transcript_id'])
            maf_df.set_value(index, 'HGVSc', hgvsc)
            maf_df.set_value(index, 'codon_change', response_json['codon_change'])
            maf_df.set_value(index, 'HGVSp_Short', row['HGVSp_Short_'])

        maf_df[maf_df.index.isin(rows_to_rem)].to_csv(output_maf_prefix + '.conversion_failed.maf', sep='\t', index=False)

        maf_df = maf_df.drop(maf_df.index[[rows_to_rem]])

        maf_df[maf_df['Variant_Type'] == 'SNP'].to_csv( output_maf_prefix + '.SNP.maf',
                                                                    sep='\t',
                                                                    index=False)

        maf_df[maf_df['Variant_Type'].isin(['INS', 'DEL'])].to_csv( output_maf_prefix + '.INDEL.maf',
                                                                    sep='\t',
                                                                    index=False)

        log_out.write('' + str(len(rows_to_rem)) + ' coding mutations failed Oncotator annotation.\n')
        log_out.write('See: ' + output_maf_prefix + '.conversion_failed.maf.\n')
        log_out.write('Success')
    except Exception:
        log_out.write('Error occurred while converting MAF to oncotator format.')
        log_out.write(traceback.format_exc())
        log_out.write('Failed\n')
        exit(1)

if __name__ == '__main__':
    main()

def get_timestamp():
    return time.strftime('%H:%M:%S %s')
