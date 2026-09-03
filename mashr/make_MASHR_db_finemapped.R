# Loading libraries and defining arguments
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(RSQLite))
suppressMessages(library(data.table))
suppressMessages(library(argparse))

"%&%" <- function(a,b) paste(a,b, sep='')
driver <- dbDriver('SQLite')
parser <- ArgumentParser()
parser$add_argument('-f', '--filesdirectory', help='path of the directory with files containing MASHR outputs')
parser$add_argument('-g', '--geneannotation', help='file path of the gene annotation file')
parser$add_argument('-m', '--finemapping', help='file path of the fine-mapping input file')
parser$add_argument('-c', '--codes', help='conditions code used, separated by a hyphen ("-")')
parser$add_argument('-o', '--outpath', help='output directory path')
args <- parser$parse_args()

# working directory to where mashr results files are 
mashr_dir = args$filesdirectory

# Get conditions codes
codes <- args$codes %>% str_split(pattern='-') %>% unlist()

fm_file <- args$finemapping

# Read in fine mapping data
top_SNPs_df <- fread(fm_file, sep = '\t', header = T) %>% select(gene = ensg, snp_ID = variant_id) %>% distinct(gene, snp_ID)

# Get gene names in the MASHR output files directory
gene_list <- unique(top_SNPs_df$gene)

process_gene_data <- function(working_gene, c, top_snps_df, mashr_dir) {
  
  # A. Read mashr data (assumed to be the reference)
  mashr_file <- mashr_dir %&% "/" %&% working_gene %&% '_MASHR_beta.txt.gz'
  if (!file.exists(mashr_file)) {
    warning("Mashr file not found for gene: ", working_gene)
    return(NULL)
  }
  
  mashr_in <- fread(mashr_file) %>% 
    select(gene, snps, snp_ID, contains(c))

  final_joined_data <- mashr_in %>%
    inner_join(
      top_snps_df %>% filter(gene == working_gene),
      by = c("gene", "snp_ID")
    )
  
  return(final_joined_data)
}

# Get MASHR-adjusted betas for the top SNPs
print('INFO: Making condition-specific transcriptome models')
for (c in codes){
  print('INFO: Current condition code is ' %&% c)
  
  # Get betas for each gene
  for (working_gene in gene_list){
    # Call the new function to read, check, correct, and join the data
    mashr_in <- process_gene_data(working_gene, c, top_SNPs_df, mashr_dir)
    
    if (exists('weights_df')){
      weights_df <- rbind(weights_df, mashr_in)
    } else {weights_df <- mashr_in}
  }
  
  
  # Make column with REF and ALT alleles
  weights_df <- weights_df %>% 
    mutate(refAllele=substr(snp_ID, nchar(snp_ID)-2, nchar(snp_ID)-2), effectAllele=substr(snp_ID, nchar(snp_ID), nchar(snp_ID))) %>% 
    rename(beta=contains(c)) %>% filter(beta!=0)
  
  # Get weight data frame into PrediXcan format
  weights_df <- weights_df %>% select('gene','snps','snp_ID','refAllele','effectAllele','beta') %>%
    rename(gene=gene, rsid=snps, varID=snp_ID, ref_allele=refAllele, eff_allele=effectAllele, weight=beta) 
  
  weights_df$rsid <- paste0(gsub(":", "_", weights_df$rsid), "_b38")
  weights_df$varID <- paste0(gsub(":", "_", weights_df$varID), "_b38")
  
  #Change X:38577675:A:G to chrX:38577675:A:G
  weights_df <- weights_df %>%
    mutate(
      # Check if the snp_id starts with "X_"
      is_x_snp = str_detect(rsid, pattern = "^X_"),
      
      # Conditionally prepend "chr"
      rsid = if_else(
        is_x_snp,
        paste0("chr", rsid), # Value if TRUE (add "chr")
        rsid               # Value if FALSE (keep original)
      ),
      
      varID = if_else(
        is_x_snp,
        paste0("chr", varID), # Value if TRUE (add "chr")
        varID               # Value if FALSE (keep original)
      )
      
    ) %>%
    select(-is_x_snp) # Remove the temporary checking column
  
  # Get number of SNPs per model 
  genes_table <- table(weights_df$gene) %>% as.data.frame() %>% rename(gene=Var1, n_snps=Freq)
  
  # Make summary table for df file
  model_summaries <- fread(args$geneannotation, header=T, stringsAsFactors=F) %>% 
    select(gene_id, gene_name) %>% unique() %>% filter(gene_id %in% gene_list) %>% 
    inner_join(genes_table, by=c('gene_id'='gene'))
  model_summaries$rho_avg_squared <- rep(NA, nrow(model_summaries))
  model_summaries$zscore_pval <- rep(NA, nrow(model_summaries))
  model_summaries$zscore_qval <- rep(NA, nrow(model_summaries))
  model_summaries <- model_summaries %>% rename(gene=gene_id, genename=gene_name, n.snps.in.model=n_snps, pred.perf.R2=rho_avg_squared,
                            pred.perf.pval=zscore_pval, pred.perf.qval=zscore_qval)
  
  # Make final files
  fwrite(model_summaries, args$outpath %&% '/' %&% c %&%'_MASHR_summaries.txt', col.names=T, quote=F, sep=' ')
  fwrite(weights_df, args$outpath %&% '/' %&% c %&%'_MASHR_weights.txt', col.names=T, quote=F, sep=' ')
  conn <- dbConnect(drv = driver, args$outpath %&% '/' %&% c %&%'_MASHR.db')
  dbWriteTable(conn, 'extra', model_summaries, overwrite = TRUE)
  dbExecute(conn, "CREATE INDEX gene_model_summary ON extra (gene)")
  dbWriteTable(conn, 'weights', weights_df, overwrite = TRUE)
  dbExecute(conn, "CREATE INDEX weights_rsid ON weights (rsid)")
  dbExecute(conn, "CREATE INDEX weights_gene ON weights (gene)")
  dbExecute(conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
  dbDisconnect(conn)
  rm(weights_df)
  print('INFO: Successfully made transcriptome prediction model files for '%&% c)
}
