includer("genMOSSplus")

plugin_do <- function(base_dir, out_dir) {
  res <- .prepare_data(base_dir)
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

.prepare_data <- function(base_dir) {
  pheno_desc <- read_pheno_ordering(base_dir)
  read_snps_plink_binary(
    base_dir,
    sub(".bed", "", list.files(base_dir, pattern = "*.bed")[1]), 
    pheno_desc, 
    use.phe = phe_as_fam, 
    remove.dups = remove_dups,
    maf.threshold = maf_thresh
  )
}