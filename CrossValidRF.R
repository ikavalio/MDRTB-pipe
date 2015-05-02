#source("http://bioconductor.org/biocLite.R")
#biocLite(c("snpStats"))
library(snpStats)
library(randomForest)

dataset_name <- "datasets/dataset_1"
base <- "/home/anjenson/BioInformatics/"
scriptdir <- "/home/anjenson/BioInformatics/MDRTB-pipe/"
extra <- ""
ignore_cols <- c()

plink.files <- file.path(base, dataset_name)
file_pattern <- sub(".bed", "", list.files(plink.files, pattern = "*.bed")[1])
out.dir <- file.path(base, "rf_out", sprintf("%s-%s", dataset_name, file_pattern))

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

# "EMB"  "INH"  "RIF"  "PZA"  "STM"  "OFLO"


window_range <- 10
misses_border <- 3

out.dir <- file.path(base, "rfcv_out", sprintf("%s-%s", dataset_name, file_pattern))
  
if(!file.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

for(i in 5:ncol(data.p)) {  
  drug.test <- data.p[,i]
  rf.data <- SNPs.num[!is.na(drug.test),]
  drug.test <- factor(drug.test[!is.na(drug.test)])
  drug <- colnames(data.p)[i]
  
  cv.errors.df <- data.frame()
  
  cv.results <- rfcv(trainx = rf.data, trainy = drug.test)
  
  cv.errors.df <- rbind(cv.errors.df, unlist(cv.results$error.cv))
  rownames(cv.errors.df) <- c("original")
  colnames(cv.errors.df) <- cv.results$n.var
  
  for (window in 1:window_range) {
    tryCatch({  
      
      pairs.df <- find.resistance.pairs(rf.data, data.p[,drug], window, misses_border)
      
      if (nrow(pairs.df) == 0) {
        next
      }
      
      ignore_cols <- c(ignore_cols, as.character(pairs.df[,2]))
      
      for (j in 1:nrow(pairs.df)) {
        snp1 <- as.character(pairs.df[j,1])
        snp2 <- as.character(pairs.df[j,2])
        rf.data[,snp1] <- pair.vectors(rf.data[,snp1], rf.data[,snp2], data.p[,drug])
      }
      
      if(length(ignore_cols) > 0) {
        rf.data <- rf.data[,setdiff(colnames(rf.data), ignore_cols)]
      }    
      
      cv.results <- rfcv(trainx = rf.data, trainy = drug.test)
      
      cv.errors.df <- rbind(cv.errors.df, unlist(cv.results$error.cv))
      rownames(cv.errors.df)[nrow(cv.errors.df)] <- window
      
    },
    error = function(error_condition) {
      message(error_condition)
      print(error_condition)
    })
  }
  out.dir.ext <- file.path(out.dir, paste0(drug, extra))
  if(!file.exists(out.dir.ext)) dir.create(out.dir.ext, recursive = FALSE)
  
  out_file <- file.path(out.dir.ext, sprintf("%s.csv", drug))
  write.csv(cv.errors.df, out_file, row.names=TRUE)
}
