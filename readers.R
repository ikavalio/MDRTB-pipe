
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
#' @return R list where corresponding matrixes can be accessed as list$genotypes, list$fam, etc.

read.snps.plink.binary <- function(base_dir, 
                                   name_pattern, 
                                   phenotypes,
                                   use.phe = FALSE, 
                                   remove.dups = FALSE, 
                                   maf.threshold = 0.01,
                                   na.strings = c("-9"),
                                   verbose = TRUE) {
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
  
  SNPs.num <- as(SNPs, "numeric") / 2
  
  if(remove.dups) {
    SNPs.nodups <- t(SNPs.num)
    SNPs.nodups <- t(SNPs.nodups[!duplicated(SNPs.nodups),])
    if(verbose) {
      cat(sprintf("%d duplicated columns were removed.\n", 
                  ncol(SNPs.num) - ncol(SNPs.nodups)))
    }
    SNPs.num <- SNPs.nodups 
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