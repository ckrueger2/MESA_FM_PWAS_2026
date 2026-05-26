# source("~/training-elasticnet/PredictDB-Tutorial/code/rsidgtex_v7_nested_cv_elnet.R")
# "%&%" <- function(a,b) paste(a,b, sep='')
# 
# argv <- commandArgs(trailingOnly = TRUE)
# chrom <- argv[1]
# 
# #tiss <- argv[1]
# #chrom <- argv[2]
# 
# snp_annot_file <- "../output/snp_annot.chr" %&% chrom %&% ".txt"
# gene_annot_file <- "../output/gene_annot.parsed.txt"
# genotype_file <- "../output/genotype.chr" %&% chrom %&% ".txt"
# expression_file <- "../output/transformed_expression.txt"
# #covariates_file <- "../output/covariates.txt"
# prefix <- "worsid-model_training_testing"
# 
# main(snp_annot_file, gene_annot_file, genotype_file, expression_file, as.numeric(chrom), prefix, null_testing=FALSE)

source("~/Elasticnet/code/fm_gtex_v7_nested_cv_elnet.R")
suppressMessages(library(tictoc))

"%&%" <- function(a,b) paste(a,b, sep='')

argv <- commandArgs(trailingOnly = TRUE)
pop <- argv[1]
chrom <- argv[2]

#tiss <- argv[1]
#chrom <- argv[2]

snp_annot_file <- "~/Elasticnet/" %&% pop %&% "/data/harmonized_genotypes/genotype.chr" %&% chrom %&% ".snp_annot.txt.gz"
gene_annot_file <- "~/Elasticnet/input/geneannotation.txt"
genotype_file <- "~/Elasticnet/" %&% pop %&% "/data/harmonized_genotypes/imputing/genotype.impute.chr" %&% chrom %&% ".txt.gz"
expression_file <- "~/Elasticnet/" %&% pop %&% "/data/geneexpression.txt"
#covariates_file <- "../output/covariates.txt"
fm_file <- "~/Elasticnet/fine_map/fm_data/" %&% pop %&% "_fm.txt"
prefix <- "~/Elasticnet/" %&% pop %&% "/cis_fm"


main(snp_annot_file, gene_annot_file, genotype_file, expression_file, fm_file, as.numeric(chrom), prefix, null_testing=FALSE)


warnings()