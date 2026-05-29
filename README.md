# MESA_FM_PWAS_2026

This directory contains analysis scripts and resources used to build and evaluate fine-mapping (FM) and PWAS models for the MESA proteomics project (2026). Below is a concise summary of the current folder structure and the main scripts found in each subdirectory.

## Structure

- **Proteomics-QC/**: Proteomics quality-control and population-specific processing
	- `europop_lm.R` — population linear-model adjustment for covariates

- **WGS-QC/**: Whole-genome sequencing QC scripts and notebooks
	- `PCAIR_europop.R` — PCA / related population QC
	- `final_qc_plinkfiles.rmd` — RMarkdown documenting / creating final PLINK files
	- `x_chromosome_qc.rmd` — chromosome X QC workflow

- **elastic-net/**: Elastic net model training and nested cross-validation for GTEx-derived features
	- `fm_cistrans_gtex_tiss_chrom_training.R` - create cis+trans FM EN models
	- `fm_cistrans_gtex_v7_nested_cv_elnet.R` - create cis+trans FM EN models
	- `fm_gtex_tiss_chrom_training.R` - create cis FM EN models
	- `fm_gtex_v7_nested_cv_elnet.R` - create cis FM EN models
	- `updated_gtex_tiss_chrom_training.R` - create cis-baseline EN models
	- `updated_gtex_v7_nested_cv_elnet.R` - create cis-baseline EN models

- **mashr/**: MASHR model generation utilities
	- `make_MASHR_db_allsnps.R` — create MASHR DB using all SNPs
	- `make_MASHR_db_finemapped.R` — create MASHR DB from fine-mapped variants
	- `make_MASHR_db_top250_snps.R` — create MASHR DB using top 250 SNPs
	- `run_MASHR.R` — script to run MASHR

- **udr/**: UDR model generation utilities (alternative/complimentary model set)
	- `make_UDR_db_allsnps.R` — create UDR DB using all SNPs
	- `make_UDR_db_finemapped.R` — create MASHR DB from fine-mapped variants
	- `make_UDR_db_top250_snps.R` — create MASHR DB using top 250 SNPs
	- `run_udr.R` — script to run UDR

## Usage

- Most scripts are written for R; run with `Rscript <script.R>` or open the `.rmd` files in RStudio.
- Check the top of each script for required package lists and input/output path expectations.

## Dependencies

- R (>= 4.0 recommended) and common genomics / modeling packages (e.g., `glmnet`, `mashr`, `tidyverse`, `data.table`, `plink` for system workflows). Refer to individual scripts for exact requirements.

## Notes

- This README reflects the repository state as of the current directory listing. For details on inputs, outputs, and parameters, inspect the headers and comments inside each script.
