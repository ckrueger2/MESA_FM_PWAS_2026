# Elastic-net training inputs and outputs

This repository contains R scripts for training elastic-net gene expression
prediction models. Each wrapper takes a population name and chromosome number,
then calls the corresponding training implementation in this folder. It is recommended to run a for loop to run all chromosomes sequentially, or run batches in parallel. Adjust the paths as needed in the wrappers. 

## Workflows

| Workflow | Wrapper | Training script | Output prefix |
| --- | --- | --- | --- |
| Standard cis | `en_training_wrapper.R` | `en_nested_cv_elnet.R` | `~/Elasticnet/<pop>/cis-keepam-harmpost` |
| Cis fine-mapped | `en_cis_fm_training_wrapper.R` | `en_cis_fm_nested_cv_elnet.R` | `~/Elasticnet/<pop>/cis_fm` |
| Cis+trans fine-mapped | `en_cistrans_fm_training_wrapper.R` | `en_cistrans_fm_nested_cv_elnet.R` | `~/Elasticnet/<pop>/cistrans_fm_redo` |

The `en_*_nested_cv_elnet.R` scripts contain the model-fitting functions. The
`en_*_training_wrapper.R` scripts define the input paths and call `main()`.

## What the scripts do

The wrappers are run with two command-line arguments:

- `pop`: population folder name such as `aapop`, `europop`, `chnpop`, or `hispop`
- `chrom`: chromosome number such as `1`, `2`, or `22`

Example:

```bash
Rscript en_training_wrapper.R aapop 1
```

The scripts read SNP annotation, gene annotation, genotype, and protein expression
files. They fit nested cross-validated elastic-net models and write model
summaries, SNP weights, and genotype covariance estimates.

## Required packages

The R code uses:

- `dplyr`
- `glmnet`
- `reshape2`
- `data.table`
- `methods`
- `tictoc`

Fine-mapped workflows also use the packages required by the corresponding
fine-mapping and trans-genotype functions.

## Required inputs

The standard wrapper builds these paths from `pop` and `chrom`:

| Input | Path pattern |
| --- | --- |
| SNP annotation | `~/Elasticnet/<pop>/data/split_snp/snp.annot.chr<chrom>.txt.gz` |
| Gene annotation | `~/Elasticnet/input/geneannotation.txt` |
| Genotype matrix | `~/Elasticnet/<pop>/data/split_genotype/imputing/genotype.impute.chr<chrom>.txt.gz` |
| Expression matrix | `~/Elasticnet/<pop>/data/geneexpression.txt` |

The fine-mapped wrappers use harmonized genotype inputs:

| Input | Path pattern |
| --- | --- |
| SNP annotation | `~/Elasticnet/<pop>/data/harmonized_genotypes/genotype.chr<chrom>.snp_annot.txt.gz` |
| Genotype matrix | `~/Elasticnet/<pop>/data/harmonized_genotypes/imputing/genotype.impute.chr<chrom>.txt.gz` |
| Gene annotation | `~/Elasticnet/input/geneannotation.txt` |
| Expression matrix | `~/Elasticnet/<pop>/data/geneexpression.txt` |

The output directories under each prefix must exist before running:

- `<prefix>/summary`
- `<prefix>/weights`
- `<prefix>/covariances`

## Expected input shapes

### 1) SNP annotation

The SNP annotation table is used to identify variants in the cis-window around
each gene. It should include `varID`, `pos`, `ref_vcf`, and `alt_vcf`.

Example:

| chr | varID | pos | ref_vcf | alt_vcf |
| --- | --- | ---: | --- | --- |
| 1 | chr1_100012_A_G | 100012 | A | G |
| 1 | chr1_100250_C_T | 100250 | C | T |
| 1 | chr1_101004_G_A | 101004 | G | A |

### 2) Gene annotation

The gene annotation table is filtered by chromosome and gene type. It should
include `gene_id`, `gene_name`, `gene_type`, `chr`, `start`, and `end`.

Example:

| gene_id | gene_name | gene_type | chr | start | end |
| --- | --- | --- | ---: | ---: | ---: |
| ENSG00000111111.1 | GENE1 | protein_coding | 1 | 99000 | 103000 |
| ENSG00000122222.1 | GENE2 | protein_coding | 1 | 150000 | 152000 |

### 3) Genotype matrix

The genotype file contains sample IDs in the first column and variant dosages in
the remaining columns. Column names must match the `varID` values in the SNP
annotation file.

Example:

| sample_id | chr1_100012_A_G | chr1_100250_C_T |
| --- | ---: | ---: |
| S001 | 0 | 1 |
| S002 | 1 | 0 |
| S003 | 2 | 1 |

### 4) Expression matrix

The expression file contains sample IDs in the first column (`sidno`) and gene
IDs in the remaining columns. Only genes present in the gene annotation file are
used.

Example:

| sidno | ENSG00000111111.1 | ENSG00000122222.1 |
| --- | ---: | ---: |
| S001 | 4.2 | 1.8 |
| S002 | 3.9 | 2.1 |
| S003 | 5.1 | 1.4 |

## Fine-mapped inputs

The cis fine-mapped wrapper reads:

```text
~/Elasticnet/fine_map/fm_data/<pop>_fm.txt
```

The file should include `ensg`, `variant_id`, `cs_label`, and `pip`. The
`variant_id` values are used to calculate fine-mapping penalty factors.

The cis+trans fine-mapped wrapper additionally reads:

| Input | Path pattern |
| --- | --- |
| Cis fine-mapping | `~/Elasticnet/fine_map/fm_data/<pop>_fm.txt` |
| Trans fine-mapping | `~/Elasticnet/fine_map/fm_data/trans/<pop>_fm.txt` |
| Trans genotype matrix | `~/Elasticnet/fine_map/fm_data/trans/<pop>_trans_genotypes.txt.gz` |

Example fine-mapping rows:

| ensg | cs_label | variant_id | pip |
| --- | --- | --- | ---: |
| ENSG00000111111.1 | L1 | chr1_100012_A_G | 0.90 |
| ENSG00000111111.1 | L2 | chr1_100250_C_T | 0.80 |

## Output files

Each workflow writes the following files below its workflow-specific prefix:

| Output | Path pattern | Content |
| --- | --- | --- |
| Model summaries | `<prefix>/summary/model_chr<chrom>_model_summaries.txt` | Gene-level model parameters and nested cross-validation performance |
| SNP weights | `<prefix>/weights/model_chr<chrom>_weights.txt` | Fitted SNP weights with gene, variant, reference, alternate, and beta columns |
| Covariances | `<prefix>/covariances/model_chr<chrom>_covariances.txt` | Pairwise genotype covariance values for model SNPs |
| Chromosome summary | `<prefix>/summary/model_chr<chrom>_tiss_chr_summary.txt` | Sample count, chromosome, random seed, and gene count |

For the standard, cis fine-mapped, and cis+trans fine-mapped workflows, replace
`<prefix>` with `cis-keepam-harmpost`, `cis_fm`, or `cistrans_fm_redo`, respectively.

## Fine-mapped models

Run the cis fine-mapped workflow with:

```bash
Rscript en_cis_fm_training_wrapper.R aapop 1
```

Run the cis+trans fine-mapped workflow with:

```bash
Rscript en_cistrans_fm_training_wrapper.R aapop 1
```

The fine-mapped workflows use the fine-mapping probabilities to construct
penalty factors. The cis+trans workflow combines cis and trans variants and
reports separate cis and trans SNP counts in the model summaries.

## Notes

- The default cis-window is 1 Mb.
- The default minor allele frequency filter is `maf=0.01`.
- The wrappers currently use `null_testing=FALSE`.
- The wrapper source paths and data roots are configured for `~/Elasticnet`.
