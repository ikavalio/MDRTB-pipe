includer("randomForest")

plugin_do <- function(base_dir, out_dir) {
  res <- .prepare_data(base_dir)
  X <- res$d
  vars <- min(vars, ncol(X))
  drugs <- colnames(res$p)
  
  for (i in 1:ncol(res$p)) {
    drug <- drugs[i]
    drug.test <- res$p[,drug]
    rf.data <- X[!is.na(drug.test),]
    drug.test <- factor(drug.test[!is.na(drug.test)])
    
    if (length(levels(drug.test)) == 2) {
      res.rf <- randomForest(
        x = rf.data,
        y = drug.test,
        ntree = trees,
        mtry = vars,
        importance = TRUE
      )
      
      predicted <- res.rf$predicted
      pred.summary <- table(predicted, drug.test, dnn = list("prediciton", "actual"))
      pred.summary <- as.data.frame(pred.summary)
      
      nm <- function(t) file.path(out_dir, sprintf("%s.rf.%s.csv", drug, t))
      imp <- data.frame(res.rf$importance, check.names = FALSE)
      
      write.csv(imp, nm("imp"))
      write.csv(pred.summary, nm("cv"))
      write.csv(imp[order(-abs(imp$MeanDecreaseGini)),], nm("imp.ordered"))
      write.csv(res.rf$confusion, nm("confusion"))
    } else {
      echo(sprintf("Bad phenotypes for drug %s\n", drug))
    }
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