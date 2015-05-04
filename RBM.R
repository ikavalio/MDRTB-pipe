#source("http://bioconductor.org/biocLite.R")
#biocLite(c("snpStats"))
library(snpStats)
library(deepnet)

plink.files <- "D:\\work\\bio"
maf.thresh <- 0.01
drug <- "EMB"
rbm.hidden <- 300
rbm.chainlen <- 3
rbm.batchsz <- 20
rbm.learnrate <- 0.2
rbm.gibbs <- 5

pheno.desc <- c(
  "EMB", "INH", "RIF", "RIFP", 
  "PZA", "STM", "CYCL", "ETH", 
  "PARA", "AMIK", "CAPR", "KANA", 
  "OFLO", "R1-T1", "R1-T2", "R1-T3", 
  "R2-T1", "R2-T2", "R2-T3", "R2-T4",
  "R2-T5", "TOT-T1", "TOT-T2"
)

sample.data <- read.plink(
  file.path(plink.files, "merge_v2.bed"),
  file.path(plink.files, "merge_v2.bim"),
  file.path(plink.files, "merge_v2.fam"),
  na.strings = "-9"
)

# Remove SNPs with MAF < maf.thresh

SNPs.full <- sample.data$genotypes
SNPs <- SNPs.full[,col.summary(SNPs.full)$MAF >= maf.thresh]

cat(sprintf("%d SNPs were eliminated because of MAF < 0.01", 
            ncol(SNPs.full) - ncol(SNPs)))

fam.names <- c(colnames(sample.data$fam)[1:5], pheno.desc)
colnames(sample.data$fam) <- fam.names

SNPs.num <- as(SNPs, "numeric") / 2

drug.test <- sample.data$fam[,which(colnames(sample.data$fam) == drug)]
rbm.data <- data.frame(cbind(SNPs.num, drug.test))
rbm.data <- rbm.data[!is.na(drug.test),]

rbm.data.t <- t(rbm.data)
rbm.data <- t(rbm.data.t[!duplicated(rbm.data.t),])

res <- rbm.train(rbm.data, 
                 hidden = rbm.hidden, 
                 numepochs = rbm.chainlen, 
                 batchsize = rbm.batchsz,
                 learningrate = rbm.learnrate,
                 cd = rbm.gibbs)


