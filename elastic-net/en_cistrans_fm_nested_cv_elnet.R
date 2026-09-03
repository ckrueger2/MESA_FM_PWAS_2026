#! /usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(glmnet))
suppressMessages(library(reshape2))
suppressMessages(library(methods))
suppressMessages(library(data.table))
suppressMessages(library(tictoc))
"%&%" <- function(a,b) paste(a,b, sep = "")


get_filtered_snp_annot <- function(snp_annot_file_name) {
  snp_annot <- fread(snp_annot_file_name, header = TRUE, stringsAsFactors = FALSE)
  setDF(snp_annot) 
  snp_annot <- snp_annot %>%
    distinct(varID, .keep_all = TRUE)
  snp_annot
}

get_maf_filtered_genotype <- function(genotype_file_name,  maf, samples) {
  gt_df <- fread(genotype_file_name, header = T, stringsAsFactors = F)
  setDF(gt_df)  # Convert to data.frame without copying
  rownames(gt_df) <- gt_df$V1
  gt_df$V1 <- NULL
  effect_allele_freqs <- colMeans(gt_df) / 2
  gt_df <- gt_df[,which((effect_allele_freqs >= maf) & (effect_allele_freqs <= 1 - maf))]
  gt_df
}

get_gene_annotation <- function(gene_annot_file_name, chrom, gene_types=c('protein_coding', 'pseudogene', 'lincRNA')){
  gene_df <- fread(gene_annot_file_name, header = TRUE, stringsAsFactors = FALSE) 
  setDF(gene_df)
  gene_df <- gene_df %>% filter((chr == chrom) & gene_type %in% gene_types)
  gene_df
}

get_gene_type <- function(gene_annot, gene) {
  filter(gene_annot, gene_id == gene)$gene_type
}

get_gene_expression <- function(gene_expression_file_name, gene_annot) {
  expr_df <- fread(gene_expression_file_name, header = T, stringsAsFactors = F, sep = '\t', check.names = F)
  setDF(expr_df)
  rownames(expr_df) <- expr_df$sidno
  expr_df$sidno <- NULL
  colnames(expr_df) <- make.unique(colnames(expr_df))
  expr_df <- expr_df %>% select(one_of(intersect(gene_annot$gene_id, colnames(expr_df))))
  expr_df
}

get_gene_coords <- function(gene_annot, gene) {
  row <- gene_annot[which(gene_annot$gene_id == gene),]
  c(row$start, row$end)
}

get_cis_genotype <- function(gt_df, snp_annot, coords, cis_window) {
  snp_info <- snp_annot %>% filter((pos >= (coords[1] - cis_window)) & (pos <= (coords[2] + cis_window)))
  if (nrow(snp_info) == 0)
    return(NA)
  if (TRUE %in% (snp_info$varID %in% names(gt_df))) {
    cis_gt <- gt_df %>% select(one_of(intersect(snp_info$varID, colnames(gt_df))))
  } else {
    return(NA)
  }
  column_labels <- colnames(cis_gt)
  row_labels <- rownames(cis_gt)
  # Convert cis_gt to a matrix for glmnet
  cis_gt <- matrix(as.matrix(cis_gt), ncol=ncol(cis_gt))
  colnames(cis_gt) <- column_labels
  rownames(cis_gt) <- row_labels
  cis_gt
}

# Read trans genotype file (already transposed: samples as rows, variants as columns)
read_trans_genotype <- function(trans_genotype_file, maf, samples) {
  trans_gt <- fread(trans_genotype_file, header = T, stringsAsFactors = F)
  setDF(trans_gt)
  # First column is sample IDs (varID column from transposed file)
  rownames(trans_gt) <- trans_gt$varID
  trans_gt$varID <- NULL
  # No need to transpose - file is already in correct orientation
  # Filter to samples present in expression data
  trans_gt <- trans_gt[rownames(trans_gt) %in% samples, , drop = FALSE]
  # Apply MAF filter
  effect_allele_freqs <- colMeans(trans_gt) / 2
  trans_gt <- trans_gt[, which((effect_allele_freqs >= maf) & (effect_allele_freqs <= 1 - maf)), drop = FALSE]
  trans_gt
}

# Create annotation for trans SNPs from variant IDs
create_trans_snp_annot <- function(trans_varIDs) {
  # Parse variant IDs in format chr:pos:ref:alt
  annot_list <- strsplit(trans_varIDs, ":")
  trans_annot <- data.frame(
    chr = as.character(sapply(annot_list, function(x) x[1])),
    pos = as.integer(sapply(annot_list, function(x) x[2])),
    varID = trans_varIDs,
    ref_vcf = sapply(annot_list, function(x) x[3]),
    alt_vcf = sapply(annot_list, function(x) x[4]),
    stringsAsFactors = FALSE
  )
  trans_annot
}

# Get trans SNPs for a specific gene
get_trans_genotype <- function(trans_gt_full, pip_df, gene) {
  # Get trans variants for this gene from fine-mapping
  trans_variants <- pip_df %>% 
    filter(ensg == gene) %>% 
    pull(variant_id)
  
  if (length(trans_variants) == 0)
    return(NULL)
  
  # Select trans variants that are in the genotype data
  available_variants <- intersect(trans_variants, colnames(trans_gt_full))
  
  if (length(available_variants) == 0)
    return(NULL)
  
  trans_gt <- trans_gt_full[, available_variants, drop = FALSE]
  # Ensure numeric matrix for downstream glmnet
  column_labels <- colnames(trans_gt)
  row_labels <- rownames(trans_gt)
  trans_gt <- matrix(as.matrix(trans_gt), ncol=ncol(trans_gt))
  colnames(trans_gt) <- column_labels
  rownames(trans_gt) <- row_labels
  trans_gt
}

read_fm <- function(fm_file){
  pip <- fread(fm_file, header = T, sep = '\t')
  return(pip)
}

get_fm_penalty <- function(pip_df, gt, gene, pip_threshold=0.00001, filt_cluster=F){
  pip_gene <- pip_df %>% filter(ensg==gene) 
  
  if (nrow(pip_gene) == 0)
    return(NULL)
  
  pip_scores <- pip_gene$pip
  
  if(filt_cluster==T){
    pip_scores <- pip_scores %>% group_by(cluster) %>% top_n(n=1,wt=pip)
  }
  
  names(pip_scores) <- pip_gene$variant_id
  matched_pips <- pip_scores[colnames(gt)]
  matched_pips[is.na(matched_pips)] <- 0
  
  penalty_factors <- 1 - matched_pips  
  
  return(penalty_factors)
}

generate_fold_ids <- function(n_samples, n_folds=10) {
  n <- ceiling(n_samples / n_folds)
  fold_ids <- rep(1:n_folds, n)
  sample(fold_ids[1:n_samples])
}

calc_R2 <- function(y, y_pred) {
  tss <- sum(y**2)
  rss <- sum((y - y_pred)**2)
  1 - rss/tss
}

calc_corr <- function(y, y_pred) {
  sum(y*y_pred) / (sqrt(sum(y**2)) * sqrt(sum(y_pred**2)))
}

nested_cv_elastic_net_perf <- function(x, y, n_samples, n_train_test_folds, n_k_folds, alpha, samples, penalty) {
  # Gets performance estimates for k-fold cross-validated elastic-net models.
  # Splits data into n_train_test_folds disjoint folds, roughly equal in size,
  # and for each fold, calculates a n_k_folds cross-validated elastic net model. Lambda parameter is
  # cross validated. Then get performance measures for how the model predicts on the hold-out
  # fold. Get the coefficient of determination, R^2, and a p-value, where the null hypothesis
  # is there is no correlation between prediction and observed.
  #
  # The mean and standard deviation of R^2 over all folds is then reported, and the p-values
  # are combined using Fisher's method.
  R2_folds <- rep(0, n_train_test_folds)
  corr_folds <- rep(0, n_train_test_folds)
  zscore_folds <- rep(0, n_train_test_folds)
  pval_folds <- rep(0, n_train_test_folds)
  # Outer-loop split into training and test set.
  train_test_fold_ids <- generate_fold_ids(n_samples, n_folds=n_train_test_folds)
  for (test_fold in 1:n_train_test_folds) {
    train_idxs <- which(train_test_fold_ids != test_fold)
    test_idxs <- which(train_test_fold_ids == test_fold)
    x_train <- x[(rownames(x) %in% samples[train_idxs]), ]
    y_train <- y[(rownames(y) %in% rownames(x_train))]
    x_test <- x[(rownames(x) %in% samples[test_idxs]), ]
    y_test <- y[(rownames(y) %in% rownames(x_test))]
    # Inner-loop - split up training set for cross-validation to choose lambda.
    cv_fold_ids <- generate_fold_ids(length(y_train), n_k_folds)
    y_pred <- tryCatch({
      # Fit model with training data.
      fit <- cv.glmnet(x_train, y_train, nfolds = n_k_folds, alpha = alpha, type.measure='mse', foldid = cv_fold_ids,penalty.factor=penalty)
      # Predict test data using model that had minimal mean-squared error in cross validation.
      predict(fit, x_test, s = 'lambda.min')},
      # if the elastic-net model did not converge, predict the mean of the y_train (same as all non-intercept coef=0)
      error = function(cond) rep(mean(y_train), length(y_test)))
    R2_folds[test_fold] <- calc_R2(y_test, y_pred)
    # Get p-value for correlation test between predicted y and actual y.
    # If there was no model, y_pred will have var=0, so cor.test will yield NA.
    # In that case, give a random number from uniform distribution, which is what would
    # usually happen under the null.
    corr_folds[test_fold] <- ifelse(sd(y_pred) != 0, cor(y_pred, y_test), 0)
    zscore_folds[test_fold] <- atanh(corr_folds[test_fold])*sqrt(length(y_test) - 3) # Fisher transformation
    pval_folds[test_fold] <- ifelse(sd(y_pred) != 0, cor.test(y_pred, y_test)$p.value, runif(1))
  }
  R2_avg <- mean(R2_folds)
  R2_sd <- sd(R2_folds)
  rho_avg <- mean(corr_folds)
  rho_se <- sd(corr_folds)
  rho_avg_squared <- rho_avg**2
  # Stouffer's method for combining z scores.
  zscore_est <- sum(zscore_folds) / sqrt(n_train_test_folds)
  zscore_pval <- 2*pnorm(abs(zscore_est), lower.tail = FALSE)
  # Fisher's method for combining p-values: https://en.wikipedia.org/wiki/Fisher%27s_method
  pval_est <- pchisq(-2 * sum(log(pval_folds)), 2*n_train_test_folds, lower.tail = F)
  list(R2_avg=R2_avg, R2_sd=R2_sd, pval_est=pval_est, rho_avg=rho_avg, rho_se=rho_se, rho_zscore=zscore_est, rho_avg_squared=rho_avg_squared, zscore_pval=zscore_pval)
}

do_covariance <- function(gene_id, gt, varIDs) {
  model_gt <- gt[,varIDs, drop=FALSE]
  geno_cov <- cov(model_gt)
  geno_cov[lower.tri(geno_cov)] <- NA
  cov_df <- reshape2::melt(geno_cov, varnames = c("varID1", "varID2"), na.rm = TRUE) %>%
    mutate(gene=gene_id) %>%
    select(GENE=gene, VARID1=varID1, VARID2=varID2, VALUE=value) %>%
    arrange(GENE, VARID1, VARID2)
  cov_df
}

main <- function(snp_annot_file, gene_annot_file, genotype_file, expression_file,
                 cis_fm_file, trans_fm_file, trans_genotype_file, chrom, prefix, 
                 maf=0.01, n_folds=10, n_train_test_folds=5,
                 seed=NA, cis_window=1e6, alpha=0.5, null_testing=FALSE) {
  tic()
  gene_annot <- get_gene_annotation(gene_annot_file, chrom)
  expr_df <- get_gene_expression(expression_file, gene_annot)
  
  # Read fine-mapping data
  cis_pip_df <- read_fm(cis_fm_file)
  trans_pip_df <- read_fm(trans_fm_file)
  
  samples <- rownames(expr_df)
  n_samples <- length(samples)
  genes <- colnames(expr_df)
  n_genes <- length(expr_df)
  snp_annot <- get_filtered_snp_annot(snp_annot_file)
  snp_annot$chr <- as.character(snp_annot$chr)
  gt_df <- get_maf_filtered_genotype(genotype_file, maf, samples)
  
  # Read trans genotypes once for all genes
  cat("Reading trans genotype file...\n")
  trans_gt_full <- read_trans_genotype(trans_genotype_file, maf, samples)
  cat("Trans genotype dimensions:", dim(trans_gt_full), "\n")
  
  # Create trans SNP annotation from variant IDs
  cat("Creating trans SNP annotation...\n")
  trans_snp_annot <- create_trans_snp_annot(colnames(trans_gt_full))
  
  #print(head(snp_annot))
  #print(head(trans_snp_annot))
  # Combine cis and trans annotations
  combined_snp_annot <- bind_rows(snp_annot, trans_snp_annot) %>%
    distinct(varID, .keep_all = TRUE)
  cat("Combined annotation has", nrow(combined_snp_annot), "SNPs (", nrow(snp_annot), "cis +", nrow(trans_snp_annot), "trans )\n")
  
  # Set seed----
  seed <- ifelse(is.na(seed), sample(1:1000000, 1), seed)
  set.seed(seed)
  
  # Prepare output data----
  model_summary_file <- prefix %&% '/summary/model_chr' %&% chrom %&% '_model_summaries.txt'
  model_summary_cols <- c('gene_id', 'gene_name', 'gene_type', 'alpha', 'n_cis_snps', 'n_trans_snps', 
                          'n_snps_in_model', 'lambda_min_mse',
                          'test_R2_avg', 'test_R2_sd', 'cv_R2_avg', 'cv_R2_sd', 'in_sample_R2',
                          'nested_cv_fisher_pval', 'rho_avg', 'rho_se', 'rho_zscore', 'rho_avg_squared', 'zscore_pval',
                          'cv_rho_avg', 'cv_rho_se', 'cv_rho_avg_squared', 'cv_zscore_est', 'cv_zscore_pval', 'cv_pval_est')
  write(model_summary_cols, file = model_summary_file, ncol = 25, sep = '\t')
  
  weights_file <- prefix %&% '/weights/model_chr' %&% chrom %&% '_weights.txt'
  weights_col <- c('gene_id', 'varID', 'ref', 'alt', 'beta', 'snp_type')
  write(weights_col, file = weights_file, ncol = 6, sep = '\t')
  
  tiss_chr_summ_f <- prefix %&% '/summary/model_chr' %&% chrom %&% '_tiss_chr_summary.txt'
  tiss_chr_summ_col <- c('n_samples', 'chrom', 'cv_seed', 'n_genes')
  tiss_chr_summ <- data.frame(n_samples, chrom, seed, n_genes)
  colnames(tiss_chr_summ) <- tiss_chr_summ_col
  write.table(tiss_chr_summ, file = tiss_chr_summ_f, quote = FALSE, row.names = FALSE, sep = '\t')
  
  covariance_file <- prefix %&% '/covariances/model_chr' %&% chrom %&% '_covariances.txt'
  covariance_col <- c('GENE', 'VARID1', 'VARID2', 'VALUE')
  write(covariance_col, file = covariance_file, ncol = 4, sep = ' ')
  
  # Attempt to build model for each gene----
  for (i in 1:n_genes) {
    cat(i, "/", n_genes, "\n")
    gene <- genes[i]
    gene_name <- gene_annot$gene_name[gene_annot$gene_id == gene]
    gene_type <- get_gene_type(gene_annot, gene)
    coords <- get_gene_coords(gene_annot, gene)
    
    # Get cis genotype
    cis_gt <- get_cis_genotype(gt_df, snp_annot, coords, cis_window)
    
    # Get trans genotype
    trans_gt <- get_trans_genotype(trans_gt_full, trans_pip_df, gene)
    
    # Combine cis and trans genotypes
    combined_gt <- NULL
    n_cis_snps <- 0
    n_trans_snps <- 0
    
    if (!all(is.na(cis_gt))) {
      combined_gt <- cis_gt
      n_cis_snps <- ncol(cis_gt)
    }
    
    if (!is.null(trans_gt)) {
      n_trans_snps <- ncol(trans_gt)
      if (is.null(combined_gt)) {
        combined_gt <- trans_gt
      } else {
        # Ensure samples align
        common_samples <- intersect(rownames(cis_gt), rownames(trans_gt))
        
        # Prefer cis: drop overlapping trans variants
        trans_vars_to_keep <- setdiff(colnames(trans_gt), colnames(cis_gt))
        
        # Combine with no intersection (cis first, then trans unique)
        combined_gt <- cbind(
          cis_gt[common_samples, , drop = FALSE],
          trans_gt[common_samples, trans_vars_to_keep, drop = FALSE]
        )
        
        # Sanity: enforce unique columns
        combined_gt <- combined_gt[, !duplicated(colnames(combined_gt)), drop = FALSE]
      }
    }
    
    # Check if we have any SNPs
    if (is.null(combined_gt)) {
      model_summary <- c(gene, gene_name, gene_type, alpha, 0, 0, 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      write(model_summary, file = model_summary_file, append = TRUE, ncol = 25, sep = '\t')
      next
    }
    
    # Get combined penalty factors from both cis and trans fine-mapping
    combined_pip_df <- rbind(
      cis_pip_df %>% filter(ensg == gene),
      trans_pip_df %>% filter(ensg == gene)
    ) %>% distinct(variant_id, .keep_all = TRUE)
    penalty <- get_fm_penalty(combined_pip_df, combined_gt, gene)
    #print(paste0("Length of combined_gt: ", ncol(combined_gt), " and Length of penalty vector: ", length(penalty)))
    print(gene)
    print(paste0("Is Penalty Null (no fm data)?", is.null(penalty)))
    if (is.null(penalty))
    #print("Creating Penalty Vector with 1s")
    if (is.null(penalty) || length(penalty) != ncol(combined_gt)) {
      penalty <- rep(1, ncol(combined_gt))
    }
    #print(paste0("Length of combined_gt: ", ncol(combined_gt), " and New Length of penalty vector: ", length(penalty)))
    
    model_summary <- c(gene, gene_name, gene_type, alpha, n_cis_snps, n_trans_snps, 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
    
    if (ncol(combined_gt) >= 2) {
      expression_vec <- expr_df[[i]]
      #print(head(expression_vec))
      adj_expression <- as.matrix(expression_vec)
      rownames(adj_expression) <- rownames(expr_df)
      
      adj_expression <- adj_expression[rownames(adj_expression) %in% rownames(combined_gt), , drop = FALSE]
      #cat("Expression data dimensions (after filtering):", dim(adj_expression), "\n")
      #cat("Combined genotype data dimensions:", dim(combined_gt), "\n")
      cat("  - Cis SNPs:", n_cis_snps, ", Trans SNPs:", n_trans_snps, "\n")
      
      if (null_testing)
        adj_expression <- sample(adj_expression)
      
      perf_measures <- nested_cv_elastic_net_perf(combined_gt, adj_expression, n_samples, n_train_test_folds, n_folds, alpha, samples, penalty)
      R2_avg <- perf_measures$R2_avg
      R2_sd <- perf_measures$R2_sd
      pval_est <- perf_measures$pval_est
      rho_avg <- perf_measures$rho_avg
      rho_se <- perf_measures$rho_se
      rho_zscore <- perf_measures$rho_zscore
      rho_avg_squared <- perf_measures$rho_avg_squared
      zscore_pval <- perf_measures$zscore_pval
      
      # Fit on all data
        cv_fold_ids <- generate_fold_ids(length(adj_expression), n_folds)
        fit <- tryCatch(cv.glmnet(combined_gt, adj_expression, nfolds = n_folds, alpha = 0.5, type.measure='mse', foldid = cv_fold_ids, keep = TRUE,penalty.factor=penalty),
                      error = function(cond) {message('Error'); message(geterrmessage()); list()})
      
      if (length(fit) > 0) {
        cv_R2_folds <- rep(0, n_folds)
        cv_corr_folds <- rep(0, n_folds)
        cv_zscore_folds <- rep(0, n_folds)
        cv_pval_folds <- rep(0, n_folds)
        best_lam_ind <- which.min(fit$cvm)
        for (j in 1:n_folds) {
          fold_idxs <- which(cv_fold_ids == j)
          adj_expr_fold_pred <- fit$fit.preval[fold_idxs, best_lam_ind]
          cv_R2_folds[j] <- calc_R2(adj_expression[fold_idxs], adj_expr_fold_pred)
          cv_corr_folds[j] <- ifelse(sd(adj_expr_fold_pred) != 0, cor(adj_expr_fold_pred, adj_expression[fold_idxs]), 0)
          cv_zscore_folds[j] <- atanh(cv_corr_folds[j])*sqrt(length(adj_expression[fold_idxs]) - 3)
          cv_pval_folds[j] <- ifelse(sd(adj_expr_fold_pred) != 0, cor.test(adj_expr_fold_pred, adj_expression[fold_idxs])$p.value, runif(1))
        }
        cv_R2_avg <- mean(cv_R2_folds)
        cv_R2_sd <- sd(cv_R2_folds)
        adj_expr_pred <- predict(fit, as.matrix(combined_gt), s = 'lambda.min')
        training_R2 <- calc_R2(adj_expression, adj_expr_pred)
        
        cv_rho_avg <- mean(cv_corr_folds)
        cv_rho_se <- sd(cv_corr_folds)
        cv_rho_avg_squared <- cv_rho_avg**2
        cv_zscore_est <- sum(cv_zscore_folds) / sqrt(n_folds)
        cv_zscore_pval <- 2*pnorm(abs(cv_zscore_est), lower.tail = FALSE)
        cv_pval_est <- pchisq(-2 * sum(log(cv_pval_folds)), 2*n_folds, lower.tail = F)
        
        if (fit$nzero[best_lam_ind] > 0) {
          weights <- fit$glmnet.fit$beta[which(fit$glmnet.fit$beta[,best_lam_ind] != 0), best_lam_ind]
          weighted_snps <- names(fit$glmnet.fit$beta[,best_lam_ind])[which(fit$glmnet.fit$beta[,best_lam_ind] != 0)]
          
          # Determine which SNPs are cis vs trans
          cis_varids <- if (!all(is.na(cis_gt))) colnames(cis_gt) else character(0)
          snp_types <- ifelse(weighted_snps %in% cis_varids, "cis", "trans")
          
          weighted_snps_info <- combined_snp_annot %>% filter(varID %in% weighted_snps) %>% select(varID, ref_vcf, alt_vcf)
          weighted_snps_info$gene <- gene
          weighted_snps_info <- weighted_snps_info %>%
            merge(data.frame(weights = weights, varID=weighted_snps, snp_type=snp_types), by = 'varID') %>%
            select(gene, varID, ref_vcf, alt_vcf, weights, snp_type)
          write.table(weighted_snps_info, file = weights_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')
          
          covariance_df <- do_covariance(gene, combined_gt, weighted_snps_info$varID)
          write.table(covariance_df, file = covariance_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ")
          
          model_summary <- c(gene, gene_name, gene_type, alpha, n_cis_snps, n_trans_snps, fit$nzero[best_lam_ind], fit$lambda[best_lam_ind], R2_avg, R2_sd, cv_R2_avg, cv_R2_sd, training_R2, pval_est,
                             rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval, cv_rho_avg, cv_rho_se, cv_rho_avg_squared, cv_zscore_est, cv_zscore_pval, cv_pval_est)
        } else {
          model_summary <- c(gene, gene_name, gene_type, alpha, n_cis_snps, n_trans_snps, 0, fit$lambda[best_lam_ind], R2_avg, R2_sd, cv_R2_avg, cv_R2_sd, training_R2, pval_est, rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval,
                             cv_rho_avg, cv_rho_se, cv_rho_avg_squared, cv_zscore_est, cv_zscore_pval, cv_pval_est)
        }
      } else {
        model_summary <- c(gene, gene_name, gene_type, alpha, n_cis_snps, n_trans_snps, 0, NA, R2_avg, R2_sd, NA, NA, NA, pval_est, rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval,
                           NA, NA, NA, NA, NA, NA)
      }
    }
    toc()
    write(model_summary, file = model_summary_file, append = TRUE, ncol = 25, sep = '\t')
  }
}
