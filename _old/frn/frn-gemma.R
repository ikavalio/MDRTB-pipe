library(knitr)
library(data.table)
library(glmnet)
library(ggplot2)
library(Cairo)

scriptdir <- "D:\\work\\bio\\rlib\\common"
source(file.path(scriptdir, "config.R"))

source(lib.reader)
source(lib.summary)
source(lib.frn)

dataset <- "7"
drug <- "CAPR"
lambda <- 0.7

inp.base <- home.gemma
out.base <- home.frn
sign.bound <- 0.01
cor.threshold <- 0.8
p.threshold <- 0.01
s.pos <- "2"
s.neg <- "1"
corr.matrix <- "C:\\Users\\Ivan_Kavaliou@epam.com\\Application Data\\Downloads\\ld_from_plink.csv"
ignore_ds <- read.table("D:\\work\\bio\\Exclude_mutations.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
ignore_cols <- sapply(
  ignore_ds[(ignore_ds$Dataset == dataset) & (ignore_ds$Drug == drug), "Position"],
  toString
)

dataset.dirs <- func.filter.datadir(inp.base)
dataset.id <- sprintf("dataset_%s", dataset)
dataset.dir <- dataset.dirs[grepl(dataset.id, dataset.dirs)]
dataset.dir.corr <- file.path(dataset.dir, paste0(drug, "-nodup"))

drug.files <- list.files(dataset.dir.corr, full.names = TRUE)
drug.file <- drug.files[grepl("-lmm.c.assoc.txt$", drug.files)]

plink.files <- file.path(base, dataset.id)
file_pattern <- func.plink.filename(plink.files)

pheno.desc <- read.phenotype.ordering(plink.files)

plink.data <- read.snps.plink.binary(plink.files, 
                                     file_pattern, 
                                     pheno.desc, 
                                     use.phe = phe.as.fam, 
                                     remove.dups = remove.dups,
                                     maf.threshold = maf.thresh)

SNPs.num <- data.table(plink.data$d)
if (length(ignore_cols) > 0) {
  SNPs.num <- SNPs.num[,!colnames(SNPs.num) %in% ignore_cols, with = FALSE]
}
data.p <- data.table(plink.data$p)
y <- data.p[, drug, with = FALSE]
rows.valid <- !is.na(y[[1]])
y <- factor(y[rows.valid][[1]], levels = c(s.neg, s.pos))
SNPs.num <- SNPs.num[rows.valid]

test.data <- read.csv(drug.file,
                      header = TRUE,
                      sep = "\t",
                      check.names = FALSE,
                      na.strings = c("NA"),
                      stringsAsFactors = FALSE)

s <- test.data
common.ps <- intersect(colnames(SNPs.num), s$ps)
SNPs.num <- SNPs.num[, as.character(common.ps), with = FALSE]

snps.cor <- if (is.null(corr.matrix)) {
  cor(SNPs.num)
} else {
  d <- read.table(corr.matrix, header = TRUE, sep = ",", 
                  na.strings = c("NA", "nan"))
  d.m <- as.matrix(d)
  ns <- sub("X", "", colnames(d.m))
  rownames(d.m) <- ns
  colnames(d.m) <- ns
  d.m[!is.finite(d.m)] <- 0
  d.m
}

s <- s[s$ps %in% common.ps,]

stat.observed <- s$p_lrt < p.threshold
stats.table <- data.table(mutations = s$ps, relevance = stat.observed)

ms <- c()
arrf1 <- c()
best.set <- NULL
mutations <- intersect(stats.table$mutations, colnames(snps.cor))
snps.cor <- snps.cor[mutations, mutations]

frn.fit <- function(lambda, dump.to = NULL) {
  signif <- frn.find.significant(stats.table, snps.cor, lambda = lambda, corr.threshold = cor.threshold)
  ms <- length(signif)
  
  if (!is.null(dump.to)) {
    write.table(signif, dump.to, quote = FALSE, sep = ",", row.names = FALSE)
  }
  
  if (length(signif) > 0) {
    X <- as.matrix(SNPs.num[, signif, with = FALSE])
    b <- test.data[test.data$ps %in% signif, "beta"]
    
    p.raw <- -1 * X %*% b
    p <- ifelse(p.raw < 0, s.pos, s.neg)[,1]
    p <- factor(p, levels = c(s.neg, s.pos))
    
    result <- table(p, y)
    tp <- result[s.pos, s.pos]
    fp <- result[s.pos, s.neg]
    fn <- result[s.neg, s.pos]
    m1 <- tp/(tp + fp)
    m2 <- tp/(tp + fn)
    f1 <- 2 * m1 * m2 / (m1 + m2)
    f1 <- if (is.finite(f1)) f1 else 0
    c(f1, ms)
  } else {
    c(0, ms)
  }
}

if (is.na(lambda)) {
  lambdas <- seq(0.05, 1, by = 0.05)
  
  for (lambda in lambdas) {
    res <- frn.fit(lambda)
    ms <- c(ms, res[2])
    arrf1 <- c(arrf1, res[1])
  }
  edata <- data.frame(lambda = lambdas, f1 = arrf1, ms = ms)
  
  ggplot(edata, aes(x = lambda, y = ms, colour = "SNPs")) + 
    geom_line() + geom_point() + ylab(label = "number of significant SNPs")
  
  ggplot(edata, aes(x = lambda, y = f1, colour = "f1")) + 
    geom_line() + geom_point() + ylab(label = "F1 score")
} else {
  out.path <- sprintf("%s/%s/%s", out.base, dataset.id, drug)
  if (!file.exists(out.path)) dir.create(out.path, recursive = TRUE)
  res <- frn.fit(lambda, dump.to = sprintf("%s/%s.gemma.signif.txt", out.path, drug))
  
  stats <- sprintf("Lambda, %f\nF1, %.3f\nTotal, %.3f\n", lambda, res[1], res[2])
  writeLines(stats, con = sprintf("%s/%s.gemma.stat.txt", out.path, drug))
}