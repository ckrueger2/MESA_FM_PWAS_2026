# MASHR model training inputs and outputs

This repository contains R scripts for fitting MASHR models across conditions and
creating condition-specific proteomic prediction model files. The workflow has
two stages:

1. `run_MASHR.R` fits one MASHR model for each gene/protein.
2. One of the `make_MASHR_db_*.R` scripts creates PrediXcan ready SQL databases from the MASHR model fitting.

## What the scripts do

### 1) Run MASHR

`run_MASHR.R` reads one input file per gene, estimates MASHR posterior means,
posterior standard deviations, and local false sign rates (LFSRs), and writes
three compressed output files per gene.

Example:

```bash
Rscript run_MASHR.R \
  --input /path/to/mashr_input \
  --geneannotation /path/to/geneannotation.txt \
  --output /path/to/mashr_output
```

Required arguments:

- `--input` (`-i`): directory containing the per-protein MASHR input .txt.gz files
- `--geneannotation` (`-g`): gene annotation file
- `--output` (`-o`): directory for the per-protein MASHR posterior files

### 2) Create prediction model files

The database-making scripts use the same command-line arguments:

- `--filesdirectory` (`-f`): directory containing MASHR posterior files
- `--geneannotation` (`-g`): gene annotation file
- `--finemapping` (`-m`): fine-mapping input file; required by `make_MASHR_db_finemapped.R`
- `--codes` (`-c`): population codes separated by hyphens, such as `aapop-europop`
- `--outpath` (`-o`): directory for the prediction model files

Example:

```bash
Rscript make_MASHR_db_allsnps.R \
  --filesdirectory /path/to/mashr_output \
  --geneannotation /path/to/geneannotation.txt \
  --codes aapop-europop \
  --outpath /path/to/mashr_models
```

Use `make_MASHR_db_finemapped.R` or `make_MASHR_db_top250_snps.R` in place of
`make_MASHR_db_allsnps.R` for the corresponding restricted SNP models.

## Required packages

The R code uses:

- `data.table`
- `dplyr`
- `argparse`
- `mashr`
- `tictoc`
- `stringr`
- `RSQLite`

## Required inputs

### 1) Gene annotation

The gene annotation file is read to obtain the list of genes and gene names. It
must include at least `gene_id` and `gene_name`.

Example:

| gene_id | gene_name |
| --- | --- |
| ENSG00000111111.1 | GENE1 |
| ENSG00000122222.1 | GENE2 |

### 2) Per-gene MASHR input files

The MASHR input directory must contain one tab-delimited, gzip-compressed file
for each `gene_id`:

```text
<input>/<gene_id>_data.txt.gz
```

Each file must include the annotation columns `gene`, `snp_ID`, and `snps`, plus
one `*_beta` and one `*_se` column for each condition/population. The beta and standard-error
columns are used to calculate Z scores as `beta / se`.

Example:

| gene | snp_ID | snps | aapop_beta | aapop_se | europop_beta | europop_se |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| ENSG00000111111.1 | chr1:100012:A:G | 1:100012:A:G | 0.12 | 0.04 | 0.10 | 0.03 |
| ENSG00000111111.1 | chr1:100250:C:T | 1:100250:C:T | -0.08 | 0.05 | -0.06 | 0.04 |

Genes with only one SNP are skipped because MASHR requires more than one SNP.
Missing Z scores are imputed with the mean Z score for that SNP across conditions (the level of missingness in your data should be determined by you and may or may not be appropriate depending on the study).

### 3) Fine-mapping input

`make_MASHR_db_finemapped.R` takes the fine-mapping file through the required
`--finemapping` (`-m`) argument. The file can be stored anywhere accessible to R.
For example:

```text
/path/to/finemapping/top_pip_snp_finemapped_data.txt.gz
```

The file must be tab-delimited and include `ensg` and `variant_id`. Each unique
`ensg`/`variant_id` pair identifies a fine-mapped SNP retained in the model.

Example:

| ensg | variant_id |
| --- | --- |
| ENSG00000111111.1 | chr1:100012:A:G |
| ENSG00000122222.1 | chr1:100250:C:T |

The script joins the fine-mapping `variant_id` values directly to the MASHR
`snp_ID` values. The identifiers must therefore use the same allele orientation.

## MASHR outputs

For every gene with at least two SNPs, `run_MASHR.R` writes files under
`--output`:

| Output | Path pattern | Content |
| --- | --- | --- |
| Posterior betas | `<output>/<gene_id>_MASHR_beta.txt.gz` | MASHR posterior mean beta for every SNP and condition |
| Posterior standard deviations | `<output>/<gene_id>_MASHR_SD.txt.gz` | MASHR posterior standard deviation for every SNP and condition |
| Local false sign rates | `<output>/<gene_id>_MASHR_lfsr.txt.gz` | Local false sign rate for every SNP and condition |

The posterior files retain the `gene`, `snps`, and `snp_ID` columns. Condition
columns are named with the condition code and suffix, such as `aapop_beta`,
`aapop_SD`, and `aapop_lfsr`.

## Prediction model outputs

Each database-making script writes one set of files per condition code in
`--codes`:

| Output | Path pattern | Content |
| --- | --- | --- |
| Model summaries | `<outpath>/<code>_MASHR_summaries.txt` | Gene names, number of SNPs, and placeholder prediction-performance fields |
| SNP weights | `<outpath>/<code>_MASHR_weights.txt` | Non-zero MASHR posterior beta weights in PrediXcan format |
| SQLite model database | `<outpath>/<code>_MASHR.db` | `extra` summary table and `weights` table with indexes |

The weights file and the `weights` SQLite table contain `gene`, `rsid`, `varID`,
`ref_allele`, `eff_allele`, and `weight`. SNP identifiers are converted from
colon-separated form to underscore-separated form and receive the `_b38` suffix (i.e. chr1_57908200_G_C_b38).

## All-SNP models

Use `make_MASHR_db_allsnps.R` to include every SNP with a non-zero posterior beta
from the MASHR beta files.

```bash
Rscript make_MASHR_db_allsnps.R \
  --filesdirectory /path/to/mashr_output \
  --geneannotation /path/to/geneannotation.txt \
  --codes aapop-europop \
  --outpath /path/to/mashr_models/allsnps
```

## Fine-mapped models

Use `make_MASHR_db_finemapped.R` to retain only gene/SNP pairs listed in
`finemapping/top_pip_snp_finemapped_data.txt.gz`. The script uses the MASHR beta
posterior for the retained SNPs and writes the standard summary, weights, and
SQLite outputs described above.

```bash
Rscript make_MASHR_db_finemapped.R \
  --filesdirectory /path/to/mashr_output \
  --geneannotation /path/to/geneannotation.txt \
  --finemapping /path/to/finemapping/top_pip_snp_finemapped_data.txt.gz \
  --codes aapop-europop \
  --outpath /path/to/mashr_models/finemapped
```

## Top 250 SNP models

Use `make_MASHR_db_top250_snps.R` to select the 250 SNPs with the lowest LFSR for
each gene and condition. The script reads the condition-specific LFSR column,
such as `<code>_lfsr`, joins those SNPs to the posterior beta file, and then
removes zero-weight SNPs.

```bash
Rscript make_MASHR_db_top250_snps.R \
  --filesdirectory /path/to/mashr_output \
  --geneannotation /path/to/geneannotation.txt \
  --codes aapop-europop \
  --outpath /path/to/mashr_models/top250
```

## Notes

- `--codes` must use the same condition names as the prefixes in the beta, SE,
  and LFSR columns.
- Output model files are generated separately for each condition code.
- Model performance fields in the summaries are currently populated with `NA`.
