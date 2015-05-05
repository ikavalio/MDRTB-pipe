library(deepnet)

rbm.hidden <- 300
rbm.chainlen <- 3
rbm.batchsz <- 20
rbm.learnrate <- 0.2
rbm.gibbs <- 5

out.base <- home.rbm

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
  
  for(i in 1:ncol(data.p)) {
    drug.test <- data.p[,i]
    rbm.data <- data.frame(cbind(SNPs.num, drug.test))
    drug <- colnames(data.p)[i]
    
    res <- rbm.train(rbm.data, 
                     hidden = rbm.hidden, 
                     numepochs = rbm.chainlen, 
                     batchsize = rbm.batchsz,
                     learningrate = rbm.learnrate,
                     cd = rbm.gibbs)
  }
}
