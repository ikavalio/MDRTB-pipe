plugin_do <- function(base_dir, out_dir) {
  res <- prepare_data(
    base_dir, 
    phe = phe_as_fam, 
    no_dups = remove_dups, 
    maf_thr = maf_thresh,
    ignore_cols = ignore_cols
  )
  summary_t <- NULL
  drugs <- colnames(res$p)
  
  for (i in 1:ncol(res$p)) {
    drug <- drugs[i]
    y <- res$p[, drug]
    
    gemma <- .get_gemma_data(base_dir, drug, colnames(res$X))
    X.e <- res$X[!is.na(y), gemma$rs]
    y <- as.character(y[!is.na(y)])
    b <- gemma$beta
    pred.raw <- -1 * X.e %*% b
    pred <- ifelse(pred.raw > 0, pheno_pos, pheno_neg)[,1]
    summary_t <- update_classification_stats(summary_t, drug, pred, y)
    p.vals.s <- sort(gemma$p_lrt)
    p.trs.pos <- min(length(p.vals.s), p_maxsel)
    p.trs.val <- min(p.vals.s[p.trs.pos], p_threshold)
    
    result <- gemma[,c("ps", "allele1", "allele0", "af", "p_lrt")]
    for (m in p_adj_m)
      result[m] <- p.adjust(gemma$p_lrt, m)
    
    write.csv(
      result, 
      file = file.path(out_dir, sprintf("%s.ps.gemma.csv", drug))
    )
    write.csv(
      result[result$p_lrt <= p.trs.val,], 
      file = file.path(out_dir, sprintf("%s.signif.gemma.csv", drug))
    )
  }
  write.csv(
    summary_t,
    file = file.path(out_dir, "summary.gemma.csv")
  )
}

.get_gemma_data <- function(base_dir, drug, allowed_cols = NULL) {
  path <- file.path(base_dir, sprintf("%s-lmm.c.assoc.txt", drug))
  test_data <- read.csv(
    path,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    na.strings = c("NA"),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(allowed_cols))
    test_data[test_data$rs %in% allowed_cols,]
  else 
    test_data
}
