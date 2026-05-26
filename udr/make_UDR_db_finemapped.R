# Loading libraries and defining arguments
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(RSQLite))
suppressMessages(library(data.table))
suppressMessages(library(argparse))
suppressMessages(library(tidyr))

"%&%" <- function(a,b) paste(a,b, sep='')
driver <- dbDriver('SQLite')
parser <- ArgumentParser()
parser$add_argument('-f', '--filesdirectory', help='path of the directory with files containing UDR outputs')
parser$add_argument('-g', '--geneannotation', help='file path of the gene annotation file')
parser$add_argument('-c', '--codes', help='conditions code used, separated by a hyphen ("-")')
parser$add_argument('-o', '--outpath', help='output directory path')
args <- parser$parse_args()

# working directory to where UDR results files are 
udr_dir = args$filesdirectory

# Get conditions codes
codes <- args$codes %>% str_split(pattern='-') %>% unlist()

# Figure out what is the most significant SNP per gene, per pop
print('INFO: Assessing top SNPs per conditions for each gene')

# Get the fine mapped file. Will change depending on cis or cis+trans SNP input strategies.
fm_file <- udr_dir %&% "/" %&% "finemapping" %&% "/" %&% "top_pip_snp_finemapped_data.txt.gz"

# Read in fine mapping data
top_SNPs_df <- fread(fm_file, sep = '\t', header = T) %>% select(gene = ensg, snp_ID = variant_id) %>% distinct(gene, snp_ID)

# Get gene names in the MASHR output files directory
gene_list <- unique(top_SNPs_df$gene)

#Function to correct strand flips from top_SNPs_df with UDR as reference (we determined that the fine-mapped SNPs had some strand flips (due to using plink 1.9) that needed to be corrected for the join with UDR data, which is used for the weights in the final model files, the UDR data is correctly harmonized with the reference genome and thus is used as the reference for allele information in the strand flip correction) 
process_gene_data <- function(working_gene, c, top_snps_df, udr_dir) {
  
  # A. Read UDR data (assumed to be the reference)
  udr_file <- udr_dir %&% "/" %&% working_gene %&% '_udr_beta.txt.gz'
  if (!file.exists(udr_file)) {
    warning("UDR file not found for gene: ", working_gene)
    return(NULL)
  }
  
  udr_in <- fread(udr_file) %>% 
    select(gene, snps, snp_ID, contains(c))
  
  # B. Extract Alleles from udr_in (Reference)
  # snp_ID format: chr1:196855643:T:C --> Ref:T, Eff:C
  udr_ref <- udr_in %>%
    # Use separate() to split by the last two colons.
    separate(snp_ID, 
             into = c("chr", "pos", "udr_ref", "udr_eff"), 
             sep = ":", 
             extra = "merge", # Merges extra fields into the first one
             fill = "right"   # Fills from the right
    ) %>%
    # Now, udr_ref and udr_eff are correctly isolated.
    select(gene, chr, pos, udr_ref, udr_eff) 

  # C. Extract Alleles from top_snps_df (Target to be corrected)
  # Assuming top_snps_df has 'variant_id' column in the same format: chr1:196855643:T:C --> Ref:T, Eff:C
  top_snps_alleles <- top_SNPs_df %>% 
    # Create temporary columns for the alleles using separation from the right
    separate(snp_ID, 
             into = c("chr", "pos", "top_ref", "top_eff"), 
             sep = ":", 
             extra = "merge", 
             fill = "right",
             remove = FALSE # Keep the original snp_ID column
    ) %>%
    select(gene, chr, pos, top_ref, top_eff)
  
  # D. Merge and Check for Strand Flips
  
  # Merge udr_ref alleles with top_snps_alleles based on common identifiers (gene, snp_ID)
  # Note: The 'gene' column is used here assuming top_snps_df is pre-filtered by gene.
  # If top_snps_df contains all genes, you might join only by 'snp_ID' and drop 'gene' from the selection.
  
  flip_check_df <- inner_join(
    udr_ref, 
    top_snps_alleles, 
    by = c("gene", "chr", "pos")
  ) %>%
    mutate(
      # Check 1: Perfect match (No flip needed)
      is_aligned = (udr_ref == top_ref) & (udr_eff == top_eff),
      
      # Check 2: Strand Flip required (udr is complement of top_snps)
      is_flip_required = (udr_ref == top_eff) & (udr_eff == top_ref),
      
      # Determine the corrected snp_ID
      corrected_snp_ID = case_when(
        is_aligned ~ paste0(chr,":", pos,":", udr_ref, ":", udr_eff),
        # If flip is required, swap the alleles in the ID string
        is_flip_required ~ paste0(chr,":", pos,":", udr_ref, ":", udr_eff), 
        # For unusable SNPs, use a placeholder or remove them later
        TRUE ~ NA_character_
      ),
      snp_ID = paste0(chr,":", pos,":", top_ref, ":", top_eff))
  
  # E. Apply Correction to top_snps_df for the join
  
  # Create a map of corrected snp_IDs
  correction_map <- flip_check_df %>%
    filter(!is.na(corrected_snp_ID)) %>%
    select(snp_ID, corrected_snp_ID)
  
  # Apply the correction to the specific gene subset of top_snps_df
  # This finds the rows in top_snps_df that match the current gene's SNPS and replaces the ID
  
  # You need to filter top_snps_df for the current gene first, or apply a complex join
  # The simplest approach is to use the corrected IDs for the join only:
  
  # F. Perform the corrected Inner Join
  corrected_top_snps <- top_SNPs_df %>%
    inner_join(correction_map, by = "snp_ID") %>%
    select(-snp_ID) %>% # Remove the old ID
    rename(snp_ID = corrected_snp_ID) # Use the corrected ID for the join
  
  # Perform the inner join with the udr input using the corrected snp_ID
  final_joined_data <- udr_in %>% 
    inner_join(corrected_top_snps, by=c('gene', 'snp_ID'))
  
  return(final_joined_data)
}

# Get UDR-adjusted betas for the top SNPs
print('INFO: Making condition-specific transcriptome models')
for (c in codes){
  print('INFO: Current condition code is ' %&% c)
  
  # Get betas for each gene
  for (working_gene in gene_list){
    # Call the new function to read, check, correct, and join the data
    udr_in <- process_gene_data(working_gene, c, top_SNPs_df, udr_dir)
    
    if (exists('weights_df')){
      weights_df <- rbind(weights_df, udr_in)
    } else {weights_df <- udr_in}
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
  fwrite(model_summaries, args$outpath %&% '/' %&% c %&%'_UDR_summaries.txt', col.names=T, quote=F, sep=' ')
  fwrite(weights_df, args$outpath %&% '/' %&% c %&%'_UDR_weights.txt', col.names=T, quote=F, sep=' ')
  conn <- dbConnect(drv = driver, args$outpath %&% '/' %&% c %&%'_UDR.db')
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
