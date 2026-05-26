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
parser$add_argument('-c', '--codes', help='conditions code used, separated by a hyphen ("-")')
parser$add_argument('-o', '--outpath', help='output directory path')
args <- parser$parse_args()

# working directory to where mashr results files are 
mashr_dir = args$filesdirectory

# Get conditions codes
codes <- args$codes %>% str_split(pattern='-') %>% unlist()

# Get gene names in the MASHR output files directory
gene_list <- fread(args$geneannotation, header = T,stringsAsFactors=F)
gene_list <- gene_list$gene_id

# Figure out what is the most significant SNP per gene, per pop
print('INFO: Assessing top SNPs per conditions for each gene')

# Get MASHR-adjusted betas for the top SNPs
print('INFO: Making condition-specific transcriptome models')
# Loop through conditions (e.g., populations)
for (c in codes){
  print('INFO: Current condition code is ' %&% c)
  
  weights_df <- data.frame()
  
  # Loop through each gene
  for (working_gene in gene_list){
    
    #Read in weights
    mashr_in <- read.table(
      mashr_dir %&% "/" %&% working_gene %&% '_MASHR_beta.txt.gz', 
      header = T
    ) %>% select(gene, snp_ID, snps, contains(c)) 
    
    #Read in lfsr data and filter by top 250 snp by lowest lfsr 
    lfsr_file <- file.path(mashr_dir %&% "/" %&% working_gene %&% "_MASHR_lfsr.txt.gz")
    
    if (file.exists(lfsr_file)) {
      lfsr_data <- fread(lfsr_file, header = TRUE)
      lfsr_col <- paste0(c, "_lfsr")
      
      if (!(lfsr_col %in% colnames(lfsr_data))) {
        warning("Missing LFSR column: ", lfsr_col, " in gene: ", working_gene)
        next
      }
      
      # Filter LFSR data: Select, rename, sort, and keep top 250 SNP_IDs
      top_250_snps <- lfsr_data %>%
        select(snp_ID, !!sym(lfsr_col)) %>%
        rename(lfsr = !!sym(lfsr_col)) %>%
        arrange(lfsr) %>%
        slice_head(n = 250) %>%
        select(snp_ID) # Only need the IDs for joining
      
      # --- 3. FILTER WEIGHTS AGAINST THE GENE-SPECIFIC TOP 250 LIST ---
      filtered_gene_weights <- mashr_in %>%
        inner_join(top_250_snps, by = "snp_ID")
      
    }
    
    # Accumulate the strictly filtered weights
    weights_df <- bind_rows(weights_df, filtered_gene_weights)
    
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
