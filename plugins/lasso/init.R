includer("glmnet")

plugin_do <- function(base_dir, out_dir) {
  res <- .prepare_data(base_dir)
  summary_t <- NULL
  summary_t_red <- NULL
  drugs <- colnames(res$p)
  
  for (i in 1:ncol(res$p)) {
    drug <- drugs[i]
    tM <- res$p[, i, drop = FALSE]
    tM <- tM[!is.na(tM[,1]), 1, drop = FALSE]
    
    # correct available factor values
    tM[,1] <- as.factor(tM[, 1])
    tM.summ <- summary(tM[,1])
    tM.noNA <- length(levels(tM[,1]))
    
    if (all(tM.summ[1:tM.noNA] > 2) && tM.noNA > 1) {
      Xt <- merge(tM, res$d, all = FALSE, by = "row.names")
      X <- as.matrix(Xt[,3:ncol(Xt)])
      y <- Xt[,2]
      ares <- .analyze(
        X, y, drug, summary_t, alpha = alpha, signif = signif_level
      )
      summary_t <- ares$summary
      snps_full <- ares$stats
      lasso_stats <- ares$lasso.stats
      if (nrow(snps_full) > 0) {
        ares_red <- .analyze(
          X[, snps_full$position], y, drug, summary_t_red, alpha = alpha,
          signif = signif_level
        )
        summary_t_red <- ares_red$summary
        snps_reduced <- ares_red$stats
      } else {
        y.r <- as.integer(y)
        summary_t_red <- update_classification_stats(
          summary_t_red, drug, rep(1, length(y.r)), y.r
        )
        snps_reduced <- data.frame(
          position = character(0), 
          bp = character(0), 
          coeficient = character(0)
        )
      }
      write.csv(
        snps_full, 
        file = file.path(out_dir, paste0(drug, ".lm.csv"))
      )
      write.csv(
        snps_reduced, 
        file = file.path(out_dir, paste0(drug, ".lm2.csv"))
      )
      write.csv(
        lasso_stats, 
        file = file.path(out_dir, paste0(drug, ".lm.stats.csv"))
      )
    } else {
      echo(sprintf("Ignore drug '%s': too few samples in one of classes", drug))
    }
  }
  write.csv(
    summary_t,
    file = file.path(out_dir, "summary.lm.csv")
  )
  write.csv(
    summary_t_red, 
    file = file.path(out_dir, "summary.lm2.csv")
  )
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

.analyze <- function(X, y, drug, summary, alpha = 0.7, signif = 0.01) {
  y.r <- as.integer(y)
  y.r.min <- min(y.r)
  y.r.max <- max(y.r)
  
  fit.regularized <- cv.glmnet(
    X, y, family = "binomial", alpha = alpha, lambda = seq(1, 0, by = -0.02), 
    nfolds = length(y)
  )
  lambda.best <- fit.regularized$lambda.min
  best.pos <- which(fit.regularized$lambda == lambda.best)
  best.cv <- fit.regularized$cvm[best.pos]
  best.nzero <- fit.regularized$nzero[best.pos]
  
  p.raw <- predict(fit.regularized, X, s = lambda.best)[,1]
  p <- ifelse(p.raw < 0, y.r.min, y.r.max)
  
  model.coef <- coef(fit.regularized, s = lambda.best)
  model.coef.nointerc <- model.coef[2:length(model.coef)]
  act.idx <- which(abs(model.coef.nointerc) > signif)
  
  act.snps <- as.integer(colnames(X)[act.idx])
  act.coef <- model.coef.nointerc[act.idx]
  
  list(
    stats = data.frame(position = act.idx, bp = act.snps, coeficient = act.coef),
    summary = update_classification_stats(summary, drug, p, y.r),
    roc = roc_create_data(y.r, p.raw),
    lasso.stats = data.frame(
      lambda = lambda.best,
      cv.err = best.cv,
      nzero  = best.nzero
    )
  )
}