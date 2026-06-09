# Elasticnet training inputs

This repository contains R scripts for training elastic-net gene expression prediction models.
The main wrapper documented here is `code/updated_gtex_tiss_chrom_training.R`, which calls
`code/updated_gtex_v7_nested_cv_elnet.R`.

## What the script does

The wrapper is run with two command-line arguments:

- `pop`: population folder name such as `aapop`, `europop`, `chnpop`, or `hispop`
- `chrom`: chromosome number such as `1`, `2`, or `22`

Example:

```bash
Rscript code/updated_gtex_tiss_chrom_training.R aapop 1
```

The script reads SNP annotation, gene annotation, genotype, and expression files, then writes model summaries, SNP weights, and covariance estimates under the population-specific output folder, that are then used to create Elastic net db models.

## Required packages

The underlying R code uses:

- `dplyr`
- `glmnet`
- `reshape2`
- `data.table`
- `methods`
- `tictoc`

## Required inputs

The wrapper builds these example paths from `pop` and `chrom`:

| Input | Path pattern |
| --- | --- |
| SNP annotation | `~/Elasticnet/<pop>/data/split_snp/snp.annot.chr<chrom>.txt.gz` |
| Gene annotation | `~/Elasticnet/input/geneannotation.txt` |
| Genotype matrix | `~/Elasticnet/<pop>/data/split_genotype/imputing/genotype.impute.chr<chrom>.txt.gz` |
| Expression matrix | `~/Elasticnet/<pop>/data/geneexpression.txt` |

The script also expects these output directories to exist before running:

- `~/Elasticnet/<pop>/cis/summary`
- `~/Elasticnet/<pop>/cis/weights`
- `~/Elasticnet/<pop>/cis/covariances`

## Expected input shapes

### 1) SNP annotation

The SNP annotation table is used to find variants inside the cis-window around each gene.
At minimum, it should include `varID`, `pos`, `ref_vcf`, and `alt_vcf`.

Example:

| chr | varID | pos | ref_vcf | alt_vcf |
| --- | --- | ---: | --- | --- |
| 1 | chr1_100012_A_G | 100012 | A | G |
| 1 | chr1_100250_C_T | 100250 | C | T |
| 1 | chr1_101004_G_A | 101004 | G | A |

### 2) Gene annotation

The gene annotation table is filtered by chromosome and gene type.
The code expects at least these columns: `gene_id`, `gene_name`, `gene_type`, `chr`, `start`, and `end`.

Example:

| gene_id | gene_name | gene_type | chr | start | end |
| --- | --- | --- | ---: | ---: | ---: |
| ENSG00000111111.1 | GENE1 | protein_coding | 1 | 99000 | 103000 |
| ENSG00000122222.1 | GENE2 | protein_coding | 1 | 150000 | 152000 |
| ENSG00000133333.1 | GENE3 | protein_coding | 1 | 210000 | 214000 |

### 3) Genotype matrix

The genotype file is read as a matrix with sample IDs in the first column and variant dosages in the remaining columns.
Column names must match the `varID` values in the SNP annotation file.

Example:

| sample_id | chr1_100012_A_G | chr1_100250_C_T | chr1_101004_G_A |
| --- | ---: | ---: | ---: |
| S001 | 0 | 1 | 2 |
| S002 | 1 | 0 | 1 |
| S003 | 2 | 1 | 0 |

### 4) Expression matrix

The expression file is read with sample IDs in the first column (`sidno`) and gene IDs in the remaining columns.
Only genes present in the gene annotation file are used.

Simulated example:

| sidno | ENSG00000111111.1 | ENSG00000122222.1 | ENSG00000133333.1 |
| --- | ---: | ---: | ---: |
| S001 | 4.2 | 1.8 | 0.3 |
| S002 | 3.9 | 2.1 | 0.5 |
| S003 | 5.1 | 1.4 | 0.2 |

## Output files

The script writes files using the prefix `~/Elasticnet/<pop>/cis`.

| Output | Path pattern | Content |
| --- | --- | --- |
| Model summaries | `~/Elasticnet/<pop>/cis/summary/model_chr<chrom>_model_summaries.txt` | Per-gene model performance and cross-validation metrics |
| SNP weights | `~/Elasticnet/<pop>/cis/weights/model_chr<chrom>_weights.txt` | Non-zero SNP weights for fitted genes |
| Covariances | `~/Elasticnet/<pop>/cis/covariances/model_chr<chrom>_covariances.txt` | Covariance values among selected SNPs |
| Chromosome summary | `~/Elasticnet/<pop>/cis/summary/model_chr<chrom>_tiss_chr_summary.txt` | Sample count, chromosome, seed, and gene count |

## Fine-mapped data

If you want to train models using fine-mapped data, use the workflow in `code/fm_gtex_tiss_chrom_training.R` or `code/fm_cistrans_gtex_tiss_chrom_training.R`.
That wrapper follows the same overall pattern as the standard script, but it also requires a fine-mapping file.

Example command:

```bash
Rscript code/fm_gtex_tiss_chrom_training.R aapop 1
```

Additional fine-mapped input example:

| Input | Path pattern |
| --- | --- |
| Fine-mapped gene/SNP file | `~/Elasticnet/fine_map/fm_data/<pop>_fm.txt` |

Example:

| ensg | cs_label | variant_id | pip |
| --- | ---: | ---: | ---: |
| ENSG00000111111.1 | L1 | chr1_100012_A_G | 0.9 |
| ENSG00000122222.1 | L1 | chr1_100250_C_T | 0.8 |
| ENSG00000133333.1 | L1 | chr1_101004_G_A | 0.85 |

Where cs_label is each credible set identified within the region of the gene. 

Subnote: if you want to use cis+trans fine-mapped data, use `code/fm_cistrans_gtex_tiss_chrom_training.R`.
That workflow needs the cis fine-mapping file above plus trans-specific inputs such as:

| Input | Path pattern |
| --- | --- |
| Cis fine-mapped file | `~/Elasticnet/fine_map/fm_data/<pop>_fm.txt` |
| Trans fine-mapped file | `~/Elasticnet/fine_map/fm_data/trans/<pop>_fm.txt` |
| Trans genotype file | `~/Elasticnet/fine_map/fm_data/trans/<pop>_trans_genotypes.txt.gz` |

Example of trans genotype file:

| sample_id | chr1_100012_A_G | chr1_100250_C_T | chr1_101004_G_A |
| --- | ---: | ---: | ---: |
| S001 | 0 | 1 | 2 |
| S002 | 1 | 0 | 1 |
| S003 | 2 | 1 | 0 |

## Notes

- The cis-window used by the R code is 1 Mb by default.
- The script filters SNPs by minor allele frequency using the `maf` argument in `main()`; the wrapper uses the default value of `0.01`.

