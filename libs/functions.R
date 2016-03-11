.cran <- "http://cran.rstudio.com/"

echo <- function(...) do.call(cat, c("####", list(...), "\n"))

#' Dependency resolution and library inclusion function.
#' 
#' @param pkgs character vector of package names. Packages will be checked,
#' installed if required and included into global environment.
#' @param repositories character vector of repository CRAN URLs.
includer <- function(pkgs, repositories = .cran, from_bioc = FALSE) {
  lib <- .libPaths()[1]
  .isPkgInstalled <- function(pkg) 
    pkg %in% installed.packages(lib.loc = lib)[, "Package"]
  if (from_bioc)
    source("http://bioconductor.org/biocLite.R")
  func <- if (from_bioc) biocLite else install.packages
  
  lapply(pkgs[!.isPkgInstalled(pkgs)], function(dep) {
    do.call(func, list(dep, dependencies = TRUE, repos = repositories, lib = lib))
    if (!.isPkgInstalled(dep))
      stop(sprintf("Unable to install the package '%s'.", dep))
  })
  
  res <- sapply(pkgs, function(x) do.call(require, list(x)))
  
  if (!all(res))
    stop(sprintf(
      "Unable to load dependencies: %s",
      do.call(paste, c(pkgs[res], list(sep = ", ")))
    ))
}

inject_args <- function(l, env = globalenv()) 
  mapply(function(nm, val) assign(nm, val, envir = env), names(l), l)

#### READERS ####

includer("snpStats", from_bioc = TRUE)

read_snps_plink_binary <- function(base_dir, 
                                   name_pattern, 
                                   phenotypes,
                                   use.phe = FALSE, 
                                   remove.dups = FALSE, 
                                   maf.threshold = 0.01,
                                   na.strings = "-9",
                                   verbose = TRUE,
                                   use.cols = get("use_cols"),
                                   ignore.cols = get("use_cols")) {
  sample.data <- read.plink(
    file.path(base_dir, sprintf("%s.bed", name_pattern)),
    file.path(base_dir, sprintf("%s.bim", name_pattern)),
    file.path(base_dir, sprintf("%s.fam", name_pattern)),
    na.strings = na.strings
  )
  if (use.phe) {
    res <- read.table(file.path(base_dir, sprintf("%s.phe", name_pattern)), sep = '\t')
    rownames(res) <- rownames(sample.data$fam)
    cnres <- c(head(colnames(sample.data$fam), 2), rep('NA', ncol(res) - 2))
    colnames(res) <- cnres
    sample.data$fam <- res
  }
  SNPs.full <- sample.data$genotypes
  
  if (!is.na(maf.threshold)) {
    SNPs <- SNPs.full[,col.summary(SNPs.full)$MAF >= maf.threshold]
    if (verbose) {
      cat(sprintf("%d SNPs were eliminated because of MAF < 0.01.\n", 
                  ncol(SNPs.full) - ncol(SNPs)))
    }
  } else {
    SNPs <- SNPs.full
  }
  
  SNPs <- filter_columns(SNPs, use.cols, ignore.cols)
  SNPs.num <- as(SNPs, "numeric") / 2
  
  if (remove.dups)
    SNPs.num <- filter_dups(SNPs.num, verbose)
  
  bnd <- ifelse(use.phe, 2, 5)
  fam.names <- c(colnames(sample.data$fam)[1:bnd], phenotypes)
  colnames(sample.data$fam) <- fam.names
  
  data.m.nona <- SNPs.num
  colnames(data.m.nona) <- sub("^.*:(\\d+)$", "\\1", colnames(data.m.nona), perl = TRUE)
  X <- data.matrix(SNPs.num)
  data.p <- sample.data$fam[(bnd + 1):ncol(sample.data$fam)]
  
  list(X = X, d = data.m.nona, p = data.p)
}

read_pheno_ordering <- function(base_dir, file_name = "fam-col-order.txt")
  strsplit(readLines(file.path(base_dir, file_name)), " ", fixed = TRUE)[[1]]

filter_dups <- function(genotype, verbose) {
  SNPs.nodups <- t(genotype)
  SNPs.nodups <- t(SNPs.nodups[!duplicated(SNPs.nodups),])
  if (verbose) {
    cat(sprintf("%d duplicated columns were removed.\n", 
                ncol(genotype) - ncol(SNPs.nodups)))
  }
  SNPs.nodups 
}

filter_columns <- function(genotype, use.cols, ignore.cols) {
  snsp.names <- colnames(genotype)
  cols.include <- if (is.null(use.cols)) snsp.names else use.cols
  cols.filter <- ignore.cols
  cols.include <- intersect(setdiff(cols.include, cols.filter), snsp.names) 
  genotype[,cols.include]
}

#### SUMMARIZERS ####

#' Update classification accuracy matrix with experiment results
#'
#' @param where existing prediction accuracy matrix or NULL (will be created automatically)
#' @param name row name
#' @param predicted vector of predicted values (should contain no more than 2 unique values)
#' @param actual vector of actual values (should contain no more than 2 unique values)
#' @param neg.value value that should be interpreted as negative outcome
#' @param pos.value value that should be interpreted as positive outcome
#' @return where matrix with appended experiment results

update_classification_stats <- function(where, name, predicted, actual, 
                                        neg.value = "1", pos.value = "2") {
  stopifnot(length(predicted) == length(actual))
  stopifnot(length(actual) > 0)
  stopifnot(length(unique(predicted)) < 3)
  stopifnot(length(unique(actual)) < 3)
  
  levels <- c(neg.value, pos.value)
  ac.f <- factor(actual, levels = levels)
  pr.f <- factor(predicted, levels = levels)
  
  stat <- table(ac.f, pr.f) # stat[i, j], i - actual, j - predicted
  update_classification_table(where, name, stat, neg.value, pos.value)
}

#' Update classification accuracy matrix with experiment results (using tabulated frequency table)
#'
#' @param where existing prediction accuracy matrix or NULL (will be created automatically)
#' @param name row name
#' @param rtable 2 x 2 experiment contingency table 
#' @param neg.value value that should be interpreted as negative outcome
#' @param pos.value value that should be interpreted as positive outcome
#' @return where matrix with appended experiment results

update_classification_table <- function(where, name, rtable, 
                                        neg.value = "1", pos.value = "2") {
  stopifnot(all(dim(rtable) == c(2, 2)))
  
  stat <- rtable # stat[i, j], i - actual, j - predicted
  tp <- stat[pos.value, pos.value]
  fp <- stat[neg.value, pos.value]
  fn <- stat[pos.value, neg.value]
  tn <- stat[neg.value, neg.value]
  
  prec <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  total <- sum(rtable)
  acc <- (tp + tn) / total
  
  row <- list()
  row[["drug"]] <- name
  row[["1"]] <- sum(stat[pos.value,])
  row[["0"]] <- sum(stat[neg.value,])
  row[["1/1"]] <- tp / total
  row[["1/0"]] <- fn / total
  row[["0/1"]] <- fp / total
  row[["0/0"]] <- tn / total
  row[["precission"]] <- prec
  row[["recall"]] <- recall
  row[["f1"]] <- 2 * prec * recall / (prec + recall)
  row[["accuracy"]] <- acc
  
  if (is.null(where)) {
    res <- data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    res <- rbind(where, row)
  }
  res
}

roc_create_data <- function(observed, predicted, pos.value = "2", 
                            neg.value = "1", thresholds = sort(predicted)) {
  observed <- factor(observed, levels = c(neg.value, pos.value))
  res <- sapply(thresholds, function(th) {
    p <- ifelse(predicted < th, neg.value, pos.value)
    p <- factor(p, levels = c(neg.value, pos.value))
    stat <- table(observed, p)
    tp <- stat[pos.value, pos.value]
    fp <- stat[neg.value, pos.value]
    fn <- stat[pos.value, neg.value]
    tn <- stat[neg.value, neg.value]
    tpr <- tp / (tp + fn)
    fpr <- fp / (fp + tn)
    
    c(x = fpr, y = tpr)
  })
  d <- data.frame(t(res))
  d <- d[order(d$x, d$y),]
  d
}