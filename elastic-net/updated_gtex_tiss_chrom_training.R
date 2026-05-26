source("/home/matt/Elasticnet/code/updated_gtex_v7_nested_cv_elnet.R")
suppressMessages(library(tictoc))

"%&%" <- function(a,b) paste(a,b, sep='')

argv <- commandArgs(trailingOnly = TRUE)
pop <- argv[1]
chrom <- argv[2]

#tiss <- argv[1]
#chrom <- argv[2]

snp_annot_file <- "~/Elasticnet/" %&% pop %&% "/data/split_snp/snp.annot.chr" %&% chrom %&% ".txt.gz"
gene_annot_file <- "~/Elasticnet/input/geneannotation.txt"
genotype_file <- "~/Elasticnet/" %&% pop %&% "/data/split_genotype/imputing/genotype.impute.chr" %&% chrom %&% ".txt.gz"
expression_file <- "~/Elasticnet/" %&% pop %&% "/data/geneexpression.txt"
#covariates_file <- "../output/covariates.txt"
prefix <- "~/Elasticnet/" %&% pop %&% "/cis-keepam-harmpost"


main(snp_annot_file, gene_annot_file, genotype_file, expression_file, as.numeric(chrom), prefix, null_testing=FALSE)


warnings()