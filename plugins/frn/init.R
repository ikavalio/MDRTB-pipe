includer(c("glmnet", "data.table"))

plugin_do <- function(base_dir, out_dir) {
  res <- prepare_data(
    base_dir, 
    phe = phe_as_fam, 
    no_dups = remove_dups, 
    maf_thr = maf_thresh,
    ignore_cols = ignore_cols
  )
  
  X <- data.table(res$d)
  data.p <- data.table(res$p)
  
  snps_cor <- if (is.null(corr_matrix)) {
    cor(X)
  } else {
    d <- read.table(corr_matrix, header = TRUE, 
                    sep = ",", na.strings = c("NA", "nan"))
    d.m <- as.matrix(d)
    ns <- sub("X", "", colnames(d.m))
    rownames(d.m) <- ns
    colnames(d.m) <- ns
    d.m[!is.finite(d.m)] <- 0
    d.m
  }
  
  type <- .detect_type(base_dir)
  if (all(!is.na(type$files))) {
    for (file in type$files) {
      proc <- type$func(file, X)
      common.ps <- intersect(colnames(X), proc$stats$mutations)
      Xt <- X[, as.character(common.ps), with = FALSE]
      
      if (!is.null(proc)) {
        ms <- c()
        y <- data.p[, proc$drug, with = FALSE]
        rows.valid <- !is.na(y[[1]])
        y <- factor(y[rows.valid][[1]], levels = c(s_neg, s_pos))
        y.table <- table(y)
        lasso.possible <- length(y.table) == 2 && all(y.table > 1)
        mutations <- intersect(proc$stats$mutations, colnames(snps_cor))
        snps_cor_t <- snps_cor[mutations, mutations]
        if (lasso.possible) {
          signif <- frn_find_significant(
            proc$stats, snps_cor_t, lambda = lambda,
            corr.threshold = cor_threshold
          )
          ms <- length(signif)
          write.table(
            signif, 
            sprintf("%s/%s.frn.imp.txt", out_dir, proc$drug), 
            quote = FALSE, 
            sep = ",", 
            row.names = FALSE
          )
          
          r <- if (length(signif) > 1) {
            Xt <- Xt[, signif, with = FALSE]
            Xt <- as.matrix(Xt[rows.valid])
            
            fit.regularized <- cv.glmnet(
              Xt, y, family = "binomial", alpha = lasso_alpha, 
              lambda = seq(1, 0, by = -0.02)
            )
            lasso.lambda.best <- fit.regularized$lambda.min
            p <- factor(
              ifelse(
                predict(fit.regularized, Xt, s = lasso.lambda.best)[,1] < 0, 
                s_neg, s_pos
              ), 
              levels = c(s_neg, s_pos)
            )
            
            result <- table(p, y)
            tp <- result[s_pos, s_pos]
            fp <- result[s_pos, s_neg]
            fn <- result[s_neg, s_pos]
            m1 <- tp/(tp + fp)
            m2 <- tp/(tp + fn)
            f1 <- 2 * m1 * m2 / (m1 + m2)
            c(f1, ms)
          } else {
            c(0, ms)
          }
          
          stats <- sprintf(
            "Lambda, %f\nF1, %.3f\nTotal, %.3f\n", lambda, r[1], r[2]
          )
          writeLines(
            stats, con = sprintf("%s/%s.frn.stat.txt", out_dir, proc$drug)
          )
        } else {
          echo("Error: unable to apply lasso for set quality estimation")
        }
      }
    }
  } else {
    echo("Error: cannot detect results type...")
  }
}

.lasso <- function(file, X) {
  allowed_ms <- colnames(X)
  drug <- sub(".*?(\\w+)[.]lm[.]csv", "\\1", file)
  s <- data.table(read.table(file, header = TRUE, stringsAsFactors = FALSE, sep = ","))
  if (length(setdiff(c("bp", "coeficient"), colnames(s))) == 0) {
    s <- s[as.character(bp) %in% allowed_ms]
    s <- s[,c("bp", "coeficient"), with = FALSE]
    
    cols.missing <- setdiff(allowed_ms, s$bp)
    missing.tbl <- data.table(bp = cols.missing, coeficient = rep(0, length(cols.missing)))
    s <- rbind(s, missing.tbl)[order(as.numeric(bp))]
    
    stat.observed <- abs(s$coeficient) > 0
    list(
      drug = drug, 
      stats = data.table(mutations = s$bp, relevance = stat.observed)
    )
  }
}

.gemma <- function(file, X) {
  allowed_ms <- colnames(X)
  drug <- sub(".*?(\\w+)-lmm[.]c[.]assoc[.]txt", "\\1", file)
  test.data <- read.csv(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    na.strings = c("NA"),
    stringsAsFactors = FALSE
  )
  s <- test.data
  
  if (length(setdiff(c("ps", "p_lrt"), colnames(s))) == 0) {
    common.ps <- intersect(allowed_ms, s$ps)
    s <- s[s$ps %in% common.ps,]
    stat.observed <- s$p_lrt < p_threshold
    list(
      drug = drug, 
      stats = data.table(mutations = s$ps, relevance = stat.observed)
    )
  }
}

.moss  <- function(file, X) {
  allowed_ms <- colnames(X)
  drug <- sub(".*?(\\w+)[.]moss[.]probs[.]csv", "\\1", file)
  s <- data.table(read.table(file, sep = ",", header = TRUE))
  if (length(setdiff(c("variable", "postIncProb"), colnames(s))) == 0) {
    s[, bp := sub("X(\\d+)", "\\1", variable)]
    s <- s[bp %in% allowed_ms]
    s <- s[,c("bp", "postIncProb"), with = FALSE]
    
    cols.missing <- setdiff(allowed_ms, s$bp)
    missing.tbl <- data.table(bp = cols.missing, postIncProb = rep(0, length(cols.missing)))
    s <- rbind(s, missing.tbl)[order(as.numeric(bp))]
    
    stat.observed <- s$postIncProb > 0
    list(
      drug = drug, 
      stats = data.table(mutations = s$bp, relevance = stat.observed)
    )
  }
}

.rf <- function(file, X) {
  allowed_ms <- colnames(X)
  drug <- sub(".*?(\\w+)[.]rf[.]imp[.]csv", "\\1", file)
  s <- data.table(read.table(file, header = TRUE, stringsAsFactors = FALSE, sep = ","))
  if (length(setdiff(c("X", "MeanDecreaseGini"), colnames(s))) == 0) {
    s <- s[as.character(X) %in% allowed_ms]
    s <- s[,c("X", "MeanDecreaseGini"), with = FALSE]
    
    cols.missing <- setdiff(allowed_ms, s$X)
    missing.tbl <- data.table(X = cols.missing, MeanDecreaseGini = rep(0, length(cols.missing)))
    s <- rbind(s, missing.tbl)[order(as.numeric(X))]
    
    stat.observed <- abs(s$MeanDecreaseGini) > 1
    list(
      drug = drug, 
      stats = data.table(mutations = s$X, relevance = stat.observed)
    )
  }
}

.detect_type <- function(base_dir) {
  lookup <- list(
    lasso = list(
      files = Sys.glob(sprintf("%s/*.lm.csv", base_dir)),
      func = .lasso
    ),
    gemma = list(
      files = Sys.glob(sprintf("%s/*-lmm.c.assoc.txt", base_dir)),
      func = .gemma
    ),
    moss = list(
      files = Sys.glob(sprintf("%s/*.moss.probs.csv", base_dir)),
      func = .moss
    ),
    rf = list(
      files = Sys.glob(sprintf("%s/*.rf.imp.csv", base_dir)),
      func = .rf
    ),
    dummy = list(
      files = NA,
      func = NULL
    )
  )
  lookup[[head(which(lapply(lookup, function(e) length(e$files)) != 0), 1)]]
}