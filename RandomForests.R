#source("http://bioconductor.org/biocLite.R")
#biocLite(c("snpStats"))
library(snpStats)
library(randomForest)

dataset_name <- "dataset_4.1"
base <- "D:\\work\\bio\\roma_13-march-2015"
scriptdir <- "D:\\work\\bio\\rlib"
extra <- ""
ignore_cols <- c()

plink.files <- file.path(base, dataset_name)
file_pattern <- sub(".bed", "", list.files(plink.files, pattern = "*.bed")[1])
out.dir <- file.path(base, "rf_out", sprintf("%s-%s", dataset_name, file_pattern))

if(!file.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

maf.thresh <- 0.01
phe.as.fam <- FALSE
remove.dups <- TRUE

maf.thresh <- 0.01
rf.trees <- 5000
rf.vars <- 500

source(file.path(scriptdir, "readers.R"))

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

rf.vars <- min(rf.vars, ncol(SNPs.num))

for(i in 1:ncol(data.p)) {
  drug.test <- data.p[,i]
  rf.data <- SNPs.num[!is.na(drug.test),]
  drug.test <- factor(drug.test[!is.na(drug.test)])
  drug <- colnames(data.p)[i]
  
  res.rf <- randomForest(x = rf.data,
                         y = drug.test,
                         ntree = rf.trees,
                         mtry = rf.vars,
                         importance = TRUE)
  
  out.dir.ext <- file.path(out.dir, paste0(drug, extra))
  if(!file.exists(out.dir.ext)) dir.create(out.dir.ext, recursive = FALSE)
  nm <- function(t) file.path(out.dir.ext, sprintf("%s.rf.%s.csv", drug, t))
  imp <- data.frame(res.rf$importance, check.names = FALSE)
  
  write.csv(imp, nm("imp"))
  write.csv(imp[order(-abs(imp$MeanDecreaseGini)),], nm("imp.ordered"))
  write.csv(res.rf$confusion, nm("confusion"))
}


