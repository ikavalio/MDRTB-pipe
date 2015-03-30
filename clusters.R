library(snpStats)
library(randomForest)

dataset_name <- "dataset_1"
base <- "D:\\work\\bio\\roma_13-march-2015"
scriptdir <- "D:\\work\\bio\\rlib"
extra <- ""
ignore_cols <- c()

plink.files <- file.path(base, dataset_name)
file_pattern <- sub(".bed", "", list.files(plink.files, pattern = "*.bed")[1])

out.file <- file.path(base, "clusters.csv")

maf.thresh <- NA
phe.as.fam <- FALSE
remove.dups <- FALSE

source(file.path(scriptdir, "readers.R"))

pheno.desc <- read.phenotype.ordering(plink.files)

plink.data <- read.snps.plink.binary(plink.files, 
                                     file_pattern, 
                                     pheno.desc, 
                                     use.phe = phe.as.fam, 
                                     remove.dups = remove.dups,
                                     maf.threshold = maf.thresh)

SNPs.num <- plink.data$d
SNPs.num.t <- t(SNPs.num)
SNPs.num.t <- data.table(SNPs.num.t, keep.rownames = TRUE)
seqcols <- tail(colnames(SNPs.num.t), -1)
snps.unique <- SNPs.num.t[!duplicated(SNPs.num.t[,seqcols, with=FALSE])]
setnames(snps.unique, "rn", "cluster")
m <- merge(SNPs.num.t, snps.unique, by = seqcols)
clusters <- m[order(as.integer(rn))][,list("SNP.pos" = rn, "Cluster.Name" = cluster)]

write.csv(clusters, out.file, row.names = FALSE, quote = FALSE)
