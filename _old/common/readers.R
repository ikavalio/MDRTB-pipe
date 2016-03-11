if(! "snpStats" %in% rownames(installed.packages())) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(c("snpStats"))
}

library(snpStats)

#' Read binary plink files (.bed, .bim, .fam) into single object.
#' Multiple phenotype columns in .fam are allowed (columns 6...Inf).
#'
#' @param base_dir root dir where files are stored
#' @param name_pattern common files name (without extension)
#' @param phenotypes extected phenotype names
#' @param use.phe store real phenotype in separate .phe file
#' @param remove.dups remove duplicated columns of SNP matrix
#' @param maf.threshold colums with 1s frequency < maf.threshold should not be considered
#' @param na.strings array of strings considered as NA
#' @param verbose show output
#' @param use.cols vector of columnames to preserve or NULL (all will be used)
#' @param ignore.cols vector of columnames to remove or NULL (none will be removed)
#' @return R list where corresponding matrixes can be accessed as list$genotypes, list$fam, etc.

read.snps.plink.binary <- function(base_dir, 
                                   name_pattern, 
                                   phenotypes,
                                   use.phe = FALSE, 
                                   remove.dups = FALSE, 
                                   maf.threshold = 0.01,
                                   na.strings = c("-9"),
                                   verbose = TRUE,
                                   use.cols = if(exists("use_cols")) get("use_cols") else NULL,
                                   ignore.cols = if(exists("ignore_cols")) get("use_cols") else NULL) {
  sample.data <- read.plink(
    file.path(base_dir, sprintf("%s.bed", name_pattern)),
    file.path(base_dir, sprintf("%s.bim", name_pattern)),
    file.path(base_dir, sprintf("%s.fam", name_pattern)),
    na.strings = na.strings
  )
  if(use.phe) {
    res <- read.table(file.path(base_dir, sprintf("%s.phe", name_pattern)), sep = '\t')
    rownames(res) <- rownames(sample.data$fam)
    cnres<- c(head(colnames(sample.data$fam), 2), rep('NA', ncol(res) - 2))
    colnames(res) <- cnres
    sample.data$fam <- res
  }
  SNPs.full <- sample.data$genotypes
  
  if(!is.na(maf.threshold)) {
    SNPs <- SNPs.full[,col.summary(SNPs.full)$MAF >= maf.threshold]
    if(verbose) {
      cat(sprintf("%d SNPs were eliminated because of MAF < 0.01.\n", 
                  ncol(SNPs.full) - ncol(SNPs)))
    }
  } else {
    SNPs <- SNPs.full
  }
  
  SNPs <- filter.columns(SNPs, use.cols, ignore.cols)
  SNPs.num <- as(SNPs, "numeric") / 2
  
  if(remove.dups) {
    SNPs.num <- filter.dups(SNPs.num, verbose)
  }
  
  bnd <- ifelse(use.phe, 2, 5)
  fam.names <- c(colnames(sample.data$fam)[1:bnd], phenotypes)
  colnames(sample.data$fam) <- fam.names
  
  data.m.nona <- SNPs.num
  colnames(data.m.nona) <- sub("^.*:(\\d+)$", "\\1", colnames(data.m.nona), perl = TRUE)
  X <- data.matrix(SNPs.num)
  data.p <- sample.data$fam[(bnd+1):ncol(sample.data$fam)]
  
  list(X = X, d = data.m.nona, p = data.p)
}

read.phenotype.ordering <- function(base_dir, file_name = "fam-col-order.txt") {
  strsplit(readLines(file.path(base_dir, file_name)), " ", fixed = TRUE)[[1]]
}

filter.dups <- function(genotype, verbose) {
  SNPs.nodups <- t(genotype)
  SNPs.nodups <- t(SNPs.nodups[!duplicated(SNPs.nodups),])
  if(verbose) {
    cat(sprintf("%d duplicated columns were removed.\n", 
                ncol(genotype) - ncol(SNPs.nodups)))
  }
  SNPs.nodups 
}

filter.columns <- function(genotype, use.cols, ignore.cols) {
  snsp.names <- colnames(genotype)
  cols.include <- if(is.null(use.cols)) snsp.names else use.cols
  cols.filter <- ignore.cols
  cols.include <- intersect(setdiff(cols.include, cols.filter), snsp.names) 
  genotype[,cols.include]
}