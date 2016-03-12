includer("genMOSSplus")

plugin_do <- function(base_dir, out_dir) {
  res <- prepare_data(
    base_dir, 
    phe = phe_as_fam, 
    no_dups = remove_dups, 
    maf_thr = maf_thresh,
    ignore_cols = ignore_cols
  )
  X <- res$d
  
  for (i in 1:ncol(res$p)) {
    drug <- colnames(res$p)[i]
    drug_test <- res$p[,i]
    moss_data <- data.frame(cbind(X, drug_test))
    moss_data <- moss_data[!is.na(drug_test),]
    
    moss_res <- MOSS.GWAS(
      alpha = alpha,
      c = c,
      cPrime = cPrime,
      q = q,
      replicates = replicas,
      maxVars = vars,
      moss_data,
      rep(2, ncol(moss_data)),
      k = cv_k
    )
    
    nm <- function(t) file.path(out_dir, sprintf("%s.moss.%s.csv", drug, t))
    
    write.csv(moss_res$topRegressions, nm("regr"))
    write.csv(moss_res$postIncProbs, nm("probs"))
    write.csv(moss_res$interactionModels, nm("mod"))
    write.csv(moss_res$cvMatrix, nm("cv"))
  }
}
