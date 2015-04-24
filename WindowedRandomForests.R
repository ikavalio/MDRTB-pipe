# !!! 
# Uncomment bioClite call when running for the first time
# !!!
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("snpStats"))
library(snpStats)
library(randomForest)

dataset_name <- "dataset_1"
base <- "/home/anjenson/BioInformatics/"
scriptdir <- "/home/anjenson/BioInformatics/MDRTB-pipe/"
extra <- ""
ignore_cols <- c()

plink.files <- file.path(base, "datasets", dataset_name)
file_pattern <- sub(".bed", "", list.files(plink.files, pattern = "*.bed")[1])
#out.dir <- file.path(base, "rf_out", sprintf("%s-%s", dataset_name, file_pattern))

#if(!file.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

maf.thresh <- 0.01
phe.as.fam <- FALSE
remove.dups <- TRUE

maf.thresh <- 0.01
rf.trees <- 5000
rf.vars <- 500

source(file.path(scriptdir, "readers.R"))
source(file.path(scriptdir, "pairing.R"))

pheno.desc <- read.phenotype.ordering(plink.files)

plink.data <- read.snps.plink.binary(plink.files, 
                                     file_pattern, 
                                     pheno.desc, 
                                     use.phe = phe.as.fam, 
                                     remove.dups = remove.dups,
                                     maf.threshold = maf.thresh)

SNPs.num <- plink.data$d
data.p <- plink.data$p

window_range <- 10
misses_border <- 3

for (window in 1:window_range) {
  out.dir <- file.path(base, "rf_out", sprintf("window-%s", window),
                       sprintf("%s-%s", dataset_name, file_pattern))
  
  if(!file.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
  
  for(i in 1:ncol(data.p)) {
    tryCatch({
      drug.test <- data.p[,i]
      rf.data <- SNPs.num[!is.na(drug.test),]
      drug.test <- factor(drug.test[!is.na(drug.test)])
      drug <- colnames(data.p)[i]
      
      SNPs.num <- plink.data$d
      data.p <- plink.data$p
      
      pairs.df <- find.resistance.pairs(SNPs.num, data.p[,drug], window, misses_border)
      
      if (nrow(pairs.df) == 0) {
        next
      }
      
      ignore_cols <- c(ignore_cols, as.character(pairs.df[,2]))
      
      for (i in 1:nrow(pairs.df)) {
        snp1 <- as.character(pairs.df[i,1])
        snp2 <- as.character(pairs.df[i,2])
        SNPs.num[,snp1] <- pair.vectors(SNPs.num[,snp1],SNPs.num[,snp2],data.p[,drug])
      }
      
      if(length(ignore_cols) > 0) {
        SNPs.num <- SNPs.num[,setdiff(colnames(SNPs.num), ignore_cols)]
      }
      
      rf.vars <- min(rf.vars, ncol(SNPs.num))
      
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
    },
    error = function(error_condition) {
      message(error_condition)
      print(error_condition)
    })
  }
}
