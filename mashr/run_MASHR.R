# Loading libraries and defining arguments
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(argparse))
suppressMessages(library(mashr))
suppressMessages(library(tictoc))
'%&%' = function(a,b) paste (a,b,sep='')
parser <- ArgumentParser()
parser$add_argument('-i', '--input', help='path of the directory with input files')
parser$add_argument('-g', '--geneannotation', help='file path of the gene annotation file')
parser$add_argument('-o', '--output', help='path of the output directory')
args <- parser$parse_args()

# Working directory to where input files are 
mashr_dir = args$input

# Get list of genes
gene_list <- fread(args$geneannotation) %>% pull(gene_id) %>% unique()  

# Run MASHR for each gene at a time
for (working_gene in gene_list){
  print('INFO: Running MASHR with gene ' %&% working_gene)
  
  # Load beta and SE dfs as matrices
  initialdata <- read.table(mashr_dir %&% '/' %&% working_gene %&% '_data.txt.gz', header=T, stringsAsFactors=F, sep = '\t')
  beta <- initialdata %>% select(contains('beta')) %>% as.matrix()
  se <- initialdata %>% select(contains('se')) %>% as.matrix() %>% abs() 

  # Obtain Z scores and scale them for MASHR input
  Z <- beta/se
  Z <- scale(Z)
  
  # Impute missing Z scores with the mean of the non-missing Z scores for that SNP across populations, using data that allows for one population to have missing data. This allows us to increase SNP coverage. 
  Z_imputed <- as.matrix(t(apply(Z, 1, function(row) {
    row[is.na(row)] <- mean(row, na.rm = TRUE)
    return(row)
  })))
  
  # For each protein, run MASHR 
  df_anno <- initialdata %>%  mutate(snps = snp_ID) %>% select(gene, snps, snp_ID)
  
  # Set up if statement to skip genes with only 1 SNP in the input data frame, as MASHR cannot be run with only 1 SNP.
  if (nrow(df_anno)==1){
    print('INFO: Not enough SNPs in the input data frame to run MASHR for ' %&% working_gene %&%'. Skipping gene')
    next
  } else {

    # Set up main MASHR data object
    data = mash_set_data(Z_imputed)
  
    # Get covariance matrices
    data.c = cov_canonical(data) # canonical
    data.pca = cov_pca(data, min(ncol(beta),nrow(beta))) # pca
    data.ed = cov_ed(data, data.pca) # data-driven
  
    # Fit model
    print('INFO: Fitting model')
    m = mash(data, Ulist=c(data.ed,data.c))
  
    # Get posterior summaries
    posterior_lfsr <- get_lfsr(m) # local false sign rate
    posterior_lfsr <- cbind(df_anno, posterior_lfsr)
    colnames(posterior_lfsr) <- gsub('_beta', '_lfsr', colnames(posterior_lfsr))
    posterior_mean <- get_pm(m) # new betas
    posterior_mean <- cbind(df_anno, posterior_mean)
    posterior_sd <- get_psd(m) # standard deviantion
    posterior_sd <- cbind(df_anno, posterior_sd)
    colnames(posterior_sd) <- gsub('_beta', '_SD', colnames(posterior_sd))

    # Write output
    fwrite(posterior_lfsr, file=args$output %&% '/' %&% working_gene %&% '_MASHR_lfsr.txt.gz', quote=F, sep=' ')
    fwrite(posterior_mean, file=args$output %&% '/' %&% working_gene %&% '_MASHR_beta.txt.gz', quote=F, sep=' ')
    fwrite(posterior_sd, file=args$output %&% '/' %&% working_gene %&% '_MASHR_SD.txt.gz', quote=F, sep=' ')
    print('INFO: Successfully ran MASHR for ' %&% working_gene)
  }
}
