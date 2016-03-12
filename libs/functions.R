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

prepare_data <- function(base_dir, phe = FALSE, 
                         no_dups = TRUE, maf_thr = 0.01, ignore_cols = NULL) {
  pheno_desc <- read_pheno_ordering(base_dir)
  read_snps_plink_binary(
    base_dir,
    sub(".bed", "", list.files(base_dir, pattern = "*.bed")[1]), 
    pheno_desc, 
    use.phe = phe, 
    remove.dups = no_dups,
    maf.threshold = maf_thr,
    ignore.cols = ignore_cols
  )
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

#### FRN ####

includer("igraph")

# E1(x) = |x - p|
frn_phi1 <- function(x, p) abs(rep(x, length(p)) - p) # p - probability of relevance

# E2(x, y) = l * |cor(x, Y)| * I(x != y)
frn_phi2 <- function(r, lmbda) 2 * lmbda * abs(r)

frn_find_significant <- function(mutations, corr.matrix, lambda = 0.2, corr.threshold = 0.7) {
  R <- corr.matrix
  R[abs(R) < corr.threshold] <- 0
  
  E1 <- frn_phi1(1, mutations$relevance) # energy for X = 1
  E0 <- frn_phi1(0, mutations$relevance) # energy for X = 0
  E2 <- frn_phi2(R, lambda)
  
  adjM.rows <- c(rownames(R), "s", "t")
  adjM.cols <- c(rownames(R), "s", "t")
  adjM <- matrix(rep(0, length(adjM.rows) * length(adjM.cols)), 
                 nrow = length(adjM.rows), ncol = length(adjM.cols), 
                 dimnames = list(adjM.rows, adjM.cols))
  
  energComparison <- E0 < E1
  energValue <- abs(E1 - E0)
  
  # connect significant mutations with 's' and not significant with 't'
  adjM["s",] <- c(ifelse(energComparison, energValue, 0), 0, 0) # last 2 els are "s" and "t"
  adjM[,"t"] <- c(ifelse(energComparison, 0, energValue), 0, 0) # last 2 els are "s" and "t"
  #E2[lower.tri(E2, diag = TRUE)] <- 0 # use only one edge between any two variables
  diag(E2) <- 0
  adjM[1:ncol(corr.matrix), 1:ncol(corr.matrix)] <- E2
  
  frn <- graph.adjacency(adjM, weighted = TRUE, mode = "directed")
  #cut <- stMincuts(frn, source = "s", target = "t")
  E(frn)$capacity <- E(frn)$weight
  flow <- graph.maxflow(frn, "s", "t")
  setdiff(V(frn)$name[flow$partition2], "t")
}
