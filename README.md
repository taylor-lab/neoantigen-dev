
# neoantigen-dev
neoantigen prediction from WES/WGS - a wrapper to Van Allen lab's neoantigen prediction pipeline (https://github.com/vanallenlab/neoantigen_calling_pipeline). 

## Install
Clone the repo and you are ready to go:
```bash
git clone https://github.com/taylor-lab/neoantigen-dev.git
```

## Usage
```
usage: neoantigen.py [-h] --config_file CONFIG_FILE [--normal_bam NORMAL_BAM]
                     --maf_file MAF_FILE --output_dir OUTPUT_DIR
                     [--hla_file HLA_FILE] [--keep_tmp_files]
                     [--force_rerun_polysolver] [--force_rerun_netmhc]

Wrapper to execute Van-Allen lab's neoantigen pipeline. In its present state,
this wrapper is primarily for preliminary analysis.  Many of paths to the scripts,
 binaries and data files are hard-coded (this will change later on).

optional arguments:
  -h, --help            show this help message and exit
  --config_file CONFIG_FILE
                        REQUIRED. 
  --normal_bam NORMAL_BAM
                        full path to normal bam file. mutually exclusive with
                        --hla_file. One of the two is required.
  --maf_file MAF_FILE   REQUIRED. expects a CMO maf file (post vcf2maf.pl)
  --output_dir OUTPUT_DIR
                        REQUIRED. output directory
  --hla_file HLA_FILE   polysolver output (winners.hla.txt) for the patient.
                        If not provided, POLYSOLVER is run. Option mutually
                        exclusive with --normal_bamfile. Either option is
                        required
  --keep_tmp_files      keeps POLYSOLVER's temporary files. for debugging
                        purposes. Note: TMP files can be more than 5GB.
                        Default: true
  --force_rerun_polysolver
                        ignores any existing polysolver output and re-runs it.
                        Default: false
  --force_rerun_netmhc  ignores any existing netMHCpan output and re-runs it.
                        Default: false

```
## Output

### HLA genotypes
```
<output_dir>/polysolver/winners.hla.txt
```

### Neoantigen binding affinties annotated MAF
```
<output_dir>/<maf_file>.netMHCpan.neoantigens.maf  (only the peptide with the highest binding affinity is selected for each mutation) 
<output_dir>/sample_processedcombinedNETMHCpan_out.txt  (all predicted peptides for each mutation)
```
The following columns are appended to the input maf.

| Column Name        | Description           |
| ------------- |:-------------|
| key      | a unique key that can be used to find other peptides predicted for the same mutation (in `sample_processedcombinedNETMHCpan_out.txt`)  |
| HLA | HLA allele to which this peptide is predicted to bind |
| pep_length | length of the peptide |
| pep_mut | mutated peptide |
| aff_mut | binding affinity of mutated peptide (nM) |
| rank_mut | Rank of the predicted affinity compared to a set of 400.000 random natural peptides (see http://www.cbs.dtu.dk/services/NetMHC/output.php) |
| pep_wt | WT peptide |
| aff_wt | binding affinity of WT peptide (nM) |
| rank_wt | Rank of the predicted affinity compared to a set of 400.000 random natural peptides (see http://www.cbs.dtu.dk/services/NetMHC/output.php) |
| total_num_neoantigens | # of other neoantigen peptides predicted for this mutation  |

## Example

```
python neoantigen.py --config_file neoantigen-luna.config \
                     --normal_bam <normal.bam> \
                     --output_dir <output_dir> \
                     --maf_file <cmo_maf_file>
```

