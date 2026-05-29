suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(Matrix))
suppressMessages(library(mvtnorm))
suppressMessages(library(ashr))
suppressMessages(library(udr))
suppressMessages(library(mashr))
suppressMessages(library(flashr))
suppressMessages(library(LaplacesDemon))
suppressMessages(library(flashier))
suppressMessages(library(argparse))

'%&%' = function(a,b) paste (a,b,sep='')
parser <- ArgumentParser()
parser$add_argument('-i', '--input', help='path of the directory with input files')
parser$add_argument('-g', '--geneannotation', help='file path of the gene annotation file')
parser$add_argument('-o', '--output', help='path of the output directory')
args <- parser$parse_args()

# working directory to where input files are 
mashr_dir = args$input

# Get list of genes
gene_list <- fread(args$geneannotation) %>% pull(gene_id) %>% unique()  

# Run UDR for each gene at a time
for (working_gene in gene_list){
  print('INFO: Running UDR with gene ' %&% working_gene)
  
  # Load beta and SE dfs as matrices
  initialdata <- read.table(mashr_dir %&% '/' %&% working_gene %&% '_data.txt.gz', header=T, stringsAsFactors=F, sep = '\t')
  beta <- initialdata %>% select(contains('beta')) %>% as.matrix()
  se <- initialdata %>% select(contains('se')) %>% as.matrix() %>% abs() 

  # Obtain Z-scores and scale them
  Z <- beta/se
  Z <- scale(Z)
  
  # Impute missing Z-scores with the mean of the non-missing Z-scores for that SNP, using data that allowed for one population to have missing data. 
  Z_imputed <- as.matrix(t(apply(Z, 1, function(row) {
    row[is.na(row)] <- mean(row, na.rm = TRUE)
    return(row)
  })))

  # For each protein, run UDR
  df_anno <- initialdata %>%  mutate(snps = snp_ID) %>% select(gene, snps, snp_ID)

  if (nrow(df_anno)==1){
    print('INFO: Not enough SNPs in the input data frame to run UDR for ' %&% working_gene %&%'. Skipping gene')
    next
  } else {
    print("Starting UDR for " %&% working_gene)
    # Set up main MASHR data object
    mash_data = mash_set_data(Z_imputed)
    
    # Get covariance matrices
    # Smart Initialization 
    smart_initialization <- function(mash_data){
      U.f <- tryCatch(cov_flash(mash_data, factors = "nonneg"),
                      error = function(e) { message("FLASH failed."); NULL })
      
      U.pca <- cov_pca(mash_data, 4)
      U.canonical <- cov_canonical(mash_data)
      
      U.init <- c(U.pca, U.canonical)  # Skip U.f if NULL
      if (!is.null(U.f)) {
        U.init <- c(U.f, U.init)
      }
      
      return(U.init)
    }
    
    #Run smart initialization function
    U.smart <- smart_initialization(mash_data)

    #Run initialization of udr object
    f0 = ud_init(mash_data, U_scaled = NULL, U_unconstrained = U.smart, n_rank1 = 0)
  
    #Run UDR - TED algorithm with IW penalty. Set maxiter to 20 and tol to 1e-2 to speed up computation, testing has shown that these parameters still allow for convergence and similar results to the default parameters, but with much faster runtime.
    ted.iw_smart = ud_fit(f0, control = list(unconstrained.update = "ted", resid.update = 'none',
                                            tol = 1e-02, tol.lik = 1e-2, lambda = ncol(beta), penalty.type = "iw", maxiter = 20), verbose=FALSE)

    # Fit model
    print('INFO: Fitting model')
    U_smartlist_ed <- get_Ulist(ted.iw_smart)
    
    # Finishing up with UDR, now running mash() with the UDR-derived covariance matrices. We are using tryCatch to catch any errors that may occur during the mash() function, and if an error occurs, we will retry once before giving up and moving on to the next gene. This is because sometimes mash() can fail due to convergence issues, and a retry can often resolve this issue.
    # First attempt
    m <- tryCatch({
      mash(mash_data, U_smartlist_ed)
    }, error = function(e1) {
      cat("First mash() attempt failed for", working_gene, "- retrying...\n")
      
      # Retry once
      tryCatch({
        mash(mash_data, U_smartlist_ed)
      }, error = function(e2) {
        cat("Second mash() attempt also failed for", working_gene, "\n")
        stop(e2)
      })
    })

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
    fwrite(posterior_lfsr, file=args$output %&% '/' %&% working_gene %&% '_udr_lfsr.txt.gz', quote=F, sep=' ')
    fwrite(posterior_mean, file=args$output %&% '/' %&% working_gene %&% '_udr_beta.txt.gz', quote=F, sep=' ')
    fwrite(posterior_sd, file=args$output %&% '/' %&% working_gene %&% '_udr_SD.txt.gz', quote=F, sep=' ')
    print('INFO: Successfully ran UDR for ' %&% working_gene)
  }
}
