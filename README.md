# MESA_FM_PWAS_2026

This directory contains analysis scripts and resources used to build and evaluate fine-mapping (FM) and PWAS models for the MESA proteomics project (2026). Below is a concise summary of the current folder structure and the main scripts found in each subdirectory.

Our developed multi-ancestry TOPMed MESA proteome prediction models, fine-mapping output, and gene coding region boundaries used for analysis are available at https://doi.org/10.5281/zenodo.15483844

## Structure

- **Proteomics-QC/**: Proteomics quality-control and population-specific processing
	- `europop_lm.R` — population linear-model adjustment for covariates

- **WGS-QC/**: Whole-genome sequencing QC scripts and notebooks
	- `PCAIR_europop.R` — PCA / related population QC
	- `final_qc_plinkfiles.rmd` — RMarkdown documenting / creating final PLINK files
	- `x_chromosome_qc.rmd` — chromosome X QC workflow

- **qtl_mapping/**: Protein quantitative trait loci mapping with TensorQTL
	- `01tensorqtl_formatting.Rmd` - formatting genotype files and phenotype .bed file for TensorQTL input
	- `02tensorqtl_mapping.Rmd` - cis- and trans- pQTL mapping (cis_nominal, cis_empirical, trans)

- **fine_mapping/**: Cross-ancestry and cross-model fine-mapping in TOPMed MESA and UKB
	- `01tensorqtl_susie.Rmd` - TOPMed MESA TensorQTL Sum of Single Effects (SuSiE) implementation across ALL, EUR, AFR, HIS, and CHN - cis and trans
	- `02susieR.Rmd` - TOPMed MESA SusieR implementation across ALL, EUR, AFR, HIS, and CHN - cis and trans
	- `03sushie.Rmd` - TOPMed MESA Sum of Shared Effects (SuShiE) - cis META output
	- `04multi.Rmd` - TOPMed MESA MultiSuSiE - cis META output
	- `05susieX.Rmd` - TOPMed MESA SuSiEx - cis META output
	- `06ukb_sum_stats_processing.Rmd` - Download and formatting of UKB pQTL summary statistics
	- `07ukb_mesa_pqtl_replication.Rmd` - pSNP replication and pi1 calculation of TOPMed MESA in UKB
	- `08ukb_mesa_finemap_replication.Rmd` - UKB SusieR, SuShiE, MultiSuSiE, and SuSiEx fine-mapping followed by precision, recall, and F1 calculations with TOPMed MESA ALL fine-mapping output

- **elastic-net/**: Elastic net model training and nested cross-validation for GTEx-derived features
	- [`README.md`](elastic-net/README.md) - inputs, outputs, workflows, and naming conventions
	- `en_cistrans_fm_training_wrapper.R` - wrapper to run cis+trans fine-mapped EN models
	- `en_cistrans_fm_nested_cv_elnet.R` - cis+trans fine-mapped EN model implementation
	- `en_cis_fm_training_wrapper.R` - wrapper to run cis fine-mapped EN models
	- `en_cis_fm_nested_cv_elnet.R` - cis fine-mapped EN model implementation
	- `en_training_wrapper.R` - wrapper to run standard cis EN models
	- `en_nested_cv_elnet.R` - standard cis EN model implementation

- **mashr/**: MASHR model generation utilities
	- [`README.md`](mashr/README.md) - inputs, posterior outputs, and database workflows
	- `make_MASHR_db_allsnps.R` — create MASHR PrediXcan db files using all SNPs
	- `make_MASHR_db_finemapped.R` — create MASHR PrediXcan db files from fine-mapped variants
	- `make_MASHR_db_top250_snps.R` — create MASHR PrediXcan db files using top 250 SNPs
	- `run_MASHR.R` — script to model data using MASHR

- **udr/**: UDR model generation utilities (alternative/complimentary model set)
	- [`README.md`](udr/README.md) - inputs, posterior outputs, and database workflows
	- `make_UDR_db_allsnps.R` — create UDR PrediXcan db files using all SNPs
	- `make_UDR_db_finemapped.R` — create UDR PrediXcan db files from fine-mapped variants
	- `make_UDR_db_top250_snps.R` — create UDR PrediXcan db files using top 250 SNPs
	- `run_udr.R` — script to model data using UDR

- **pwas/**: PWAS with TOPMed MESA trained EN, MASHR, and UDR models in AoU and Pan-UKB
	- `01pwas_in_aou.Rmd` - pre-processing of data and running PWAS analysis in All of Us (AoU) Researcher Workbench
 	- `02pwas_in_ukb.Rmd` - pre-processing of Pan-UKB phenotype summary statistics and PWAS analysis
  	- `03pwas_replication.Rmd` - Bonferroni-significant protein-trait associations in AoU that replicated in Pan-UKB
  	- `04mesa_aou_coloc.Rmd` - Colocalization analysis with TOPMed MESA cis-pQTLs and AoU protein-coding regions

## Usage

- Scripts were written in Python, R, and bash; run with `Rscript <script.R>` or open the `.Rmd` files in RStudio and paste scripts in appropriate environments, changing paths and ancestry names where necessary.
- Check the top of each script for required package lists and input/output path expectations.

## Dependencies

- R (>= 4.0 recommended) and common genomics / modeling packages (e.g., `glmnet`, `mashr`, `tidyverse`, `data.table`, `plink` for system workflows). Refer to individual scripts for exact requirements.

## Notes

- This README reflects the repository state as of the current directory listing. For details on inputs, outputs, and parameters, inspect the headers and comments inside each script.
