#source("http://bioconductor.org/biocLite.R")
#biocLite(c("snpStats"))
library(snpStats)
library(genMOSSplus)

dataset_name <- "dataset_1"
base <- "D:\\work\\bio\\roma_13-march-2015"
scriptdir <- "D:\\work\\bio\\rlib"
extra <- "-no-2155175"
ignore_cols <- c()

plink.files <- file.path(base, dataset_name)
file_pattern <- sub(".bed", "", list.files(plink.files, pattern = "*.bed")[1])
out.dir <- file.path(base, "moss_out", sprintf("%s-%s", dataset_name, file_pattern))

if(!file.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

phe.as.fam <- TRUE
maf.thresh <- 0.01
phe.as.fam <- FALSE
remove.dups <- TRUE

moss.alpha <- 1
moss.c <- 0.1
moss.cPrime <- 1e-5 
moss.q <- 0.1
moss.replicas <- 3 # 4
moss.vars <- 3
moss.cv.k <- 5

source(file.path(scriptdir, "readers.R"))

#pheno.desc <- c(
#  "EMB", "INH", "RIF", "RIFP", 
#  "PZA", "STM", "CYCL", "ETH", 
#  "PARA", "AMIK", "CAPR", "KANA", 
#  "OFLO", "R1-T1", "R1-T2", "R1-T3", 
#  "R2-T1", "R2-T2", "R2-T3", "R2-T4",
#  "R2-T5", "TOT-T1", "TOT-T2"
#  )

pheno.desc <- read.phenotype.ordering(plink.files)

plink.data <- read.snps.plink.binary(plink.files, 
                                     file_pattern, 
                                     pheno.desc, 
                                     use.phe = phe.as.fam, 
                                     remove.dups = remove.dups,
                                     maf.threshold = maf.thresh)

SNPs.num <- plink.data$d
data.p <- plink.data$p

if(length(ignore_cols) > 0) {
  SNPs.num <- SNPs.num[,setdiff(colnames(SNPs.num), ignore_cols)]
}

for(i in 1:ncol(data.p)) {
  drug.test <- data.p[,i]
  moss.data <- data.frame(cbind(SNPs.num, drug.test))
  moss.data <- moss.data[!is.na(drug.test),]
  drug <- colnames(data.p)[i]
  
  moss.r <- MOSS.GWAS(alpha = moss.alpha, 
                      c = moss.c, 
                      cPrime = moss.cPrime,
                      q = moss.q, 
                      replicates = moss.replicas, 
                      maxVars = moss.vars, 
                      moss.data,
                      rep(2, ncol(moss.data)), 
                      k = moss.cv.k)
  
  out.dir.ext <- file.path(out.dir, paste0(drug, extra))
  if(!file.exists(out.dir.ext)) dir.create(out.dir.ext, recursive = FALSE)
  nm <- function(t) file.path(out.dir.ext, sprintf("%s.moss.%s.csv", drug, t))
  
  write.csv(moss.r$topRegressions, nm("regr"))
  write.csv(moss.r$postIncProbs, nm("probs"))
  write.csv(moss.r$interactionModels, nm("mod"))
  write.csv(moss.r$crossValidation, nm("cv"))
}
