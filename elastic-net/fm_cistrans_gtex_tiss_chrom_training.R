source("/home/matt/Elasticnet/code/fm_cistrans_gtex_v7_nested_cv_elnet.R")
suppressMessages(library(tictoc))

"%&%" <- function(a,b) paste(a,b, sep='')

argv <- commandArgs(trailingOnly = TRUE)
pop <- argv[1]
chrom <- argv[2]

snp_annot_file <- "~/Elasticnet/" %&% pop %&% "/data/harmonized_genotypes/genotype.chr" %&% chrom %&% ".snp_annot.txt.gz"
gene_annot_file <- "~/Elasticnet/input/geneannotation.txt"
genotype_file <- "~/Elasticnet/" %&% pop %&% "/data/harmonized_genotypes/imputing/genotype.impute.chr" %&% chrom %&% ".txt.gz"
expression_file <- "~/Elasticnet/" %&% pop %&% "/data/geneexpression.txt"
cis_fm_file = "~/Elasticnet/fine_map/fm_data/" %&% pop %&% "_fm.txt"
trans_fm_file = "~/Elasticnet/fine_map/fm_data/trans/" %&% pop %&% "_fm.txt"
trans_genotype_file = "~/Elasticnet/fine_map/fm_data/trans/" %&% pop %&% "_trans_genotypes.txt.gz"
prefix <- "~/Elasticnet/" %&% pop %&% "/cistrans_fm_redo"

main(snp_annot_file, gene_annot_file, genotype_file, expression_file, cis_fm_file, trans_fm_file, trans_genotype_file, as.numeric(chrom), prefix, null_testing=FALSE)


warnings()