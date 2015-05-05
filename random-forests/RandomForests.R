library(randomForest)

scriptdir <- "D:\\work\\bio\\rlib\\common"
source(file.path(scriptdir, "config.R"))
source(lib.reader)

rf.trees <- 5000
rf.vars <- 500

out.base <- home.rf

for(dataset.dir in raw.files) {
  plink.files <- dataset.dir
  dataset_name <- basename(dataset.dir)
  file_pattern <- func.plink.filename(plink.files)
  out.dir <- file.path(out.base, sprintf("%s-%s", dataset_name, file_pattern))
  
  if(!file.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
  
  pheno.desc <- read.phenotype.ordering(plink.files)
  
  plink.data <- read.snps.plink.binary(plink.files,
                                       file_pattern,
                                       pheno.desc,
                                       use.phe = phe.as.fam,
                                       remove.dups = remove.dups,
                                       maf.threshold = maf.thresh)
  
  SNPs.num <- plink.data$d
  data.p <- plink.data$p
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
}

