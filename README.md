

# neoantigen-dev
MHC Class I neoantigen prediction pipeline from IM5/IM6/WES/WGS. Takes Normal `.bam` and somatic `.maf` and generates neoantigen predictions for HLA-A/B/C.

The pipeline has four main steps:

1. **Genotype HLA**. genotyping performed using POLYSOLVER. 
2. **Construct mutated peptides**.  For non-synonymous mutations, generates mutated peptide sequences based on `HGVSc`.  _NOTE_: `.maf` file should be VEP annotated using `cmo_maf2maf  --version 1.6.14 --vep-release 88` **using this EXACT VERSION**. TODO: Generate mutated sequences for fusions.
3. **Run NetMHCpan-4.0 and NetMHC-4.0**. using default parameters for each algorithm. 
4. **Post-processing**. compiles predictions from both algorithms and finds strongest binder for each non-synonymous mutation. Also, each predicted neopeptide is searched against the entire reference peptidome to make sure it is a true neopeptide. `is_in_wt_peptidome` column reflects that. TODO: Incorporate neoantigen quality from [Lukzsa et al., Nature 2017](https://www.nature.com/articles/nature24473)


## Install

Clone the repo and install any necessary Python libraries from `requirements.txt`. Note that this repo is currently only compatible with Python 3, not Python 2.x:

```bash
git clone https://github.com/taylor-lab/neoantigen-dev.git
cd neoantigen-dev
pip install -r requirements.txt
```


## Usage
NOTE: For POLYSOLVER step, the pipeline requires 8 cores. 
```
# Neoantigen prediction pipeline. Four main steps:
		(1) Genotype HLA using POLYSOLVER,
		(2) Constructed mutated peptide sequences from HGVSp/HGVSc
		(3) Run NetMHC-4.0 + NetMHCpan-4.0
		(4) Gather results and generate:
				- <sample_id>.neoantigens.maf: original maf with neoantigen prediction columns added
				- <sample_id>.all_neoantigen_predictions.txt: all the predictions made for all peptides by both the algorithms

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  --config_file CONFIG_FILE
                        See: neoantigen-luna.config in the repo
  --sample_id SAMPLE_ID
                        sample_id used to limit neoantigen prediction to
                        identify mutations associated with the patient in the
                        MAF (column 16).
  --output_dir OUTPUT_DIR
                        output directory
  --maf_file MAF_FILE   expects a CMO maf file (post vcf2maf.pl)
  --normal_bam NORMAL_BAM
                        full path to normal bam file. Either --normal_bam or
                        --hla_file are required.

Optional arguments:
  --hla_file HLA_FILE   POLYSOLVER output file (winners.hla.txt) for the
                        sample. If not provided,POLYSOLVER is run. Either
                        --normal_bam or --hla_file are required.
  --peptide_lengths PEPTIDE_LENGTHS
                        comma-separated numbers indicating the lengths of
                        peptides to generate. Default: 9,10
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
<sample_id>.neoantigens.maf. (peptide with the highest binding affinity is incorporated into the original .maf for each non-syn mutation) 
<sample_id>.all_neoantigen_predictions.txt: all the predictions made for all peptides by both the algorithms
```
The following columns are appended to the input `.maf`.

| Column Name        | Description           |
| ------------- |:-------------|
| neo_maf_identifier_key      | a unique key that can be used to find other peptides predicted for the same mutation (in `.all_neoantigen_predictions.txt`)  |
| neo_best_icore_peptide | neopeptide sequence for the strongest binder | 
| neo_best_rank | binding rank for the strongest binder | 
| neo_best_binding_affinity | binding affinity for the strongest binder | 
| neo_best_binder_classification | binding classification for the strongest binder (`Non Binder`, `Strong Binder`, `Weak Binder`) | 
| neo_best_is_in_wt_peptidome |  `TRUE`/`FALSE` indicating whether the strongest binder peptide is in the reference peptidome | 
| neo_best_algorithm | algorithm predicting the strongest binder | 
| neo_best_hla_allele | hla allele for the strongest binder | 
| neo_n_peptides_evaluated | total # of all peptides evaluated (unique icore peptides) | 
| neo_n_strong_binders |  total # of strong binders | 
| neo_n_weak_binders | total # of weak binders |

The column description for `.all_neoantigen_predictions.txt` can be found in: [http://www.cbs.dtu.dk/services/NetMHC/output.php](http://www.cbs.dtu.dk/services/NetMHC/output.php). Additional columns are:

The following columns are appended to the input `.maf`.

| Column Name        | Description           |
| ------------- |:-------------|
| binder_class      | `Non Binder`, `Strong Binder` (`rank < 0.5 or affinity < 50`), `Weak Binder` (`rank < 2 or affinity < 500`) |
| best_binder_for_icore_group | `TRUE`/`FALSE` indicating if the binding prediction is the strongest among all the HLA-alleles/algorithms for the given icore peptide. | 
| is_in_wt_peptidome | if the peptide is present in any other protein in the entire peptidome | 
| neo_maf_identifier_key | a unique key that can be used to find other peptides predicted for the same mutation (in `.neoantigens.maf`)  | 

## Example

```
python neoantigen.py --config_file neoantigen-luna.config \
                     --sample_id <sample_id> \
                     --normal_bam <normal.bam> \
                     --output_dir <output_dir> \
                     --maf_file <cmo_vep_annotated_maf_file>
```

