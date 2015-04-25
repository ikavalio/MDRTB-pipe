library(data.table)

dataset_name <- "dataset_1"
base <- "D:\\work\\bio\\roma_13-march-2015"
scriptdir <- "D:\\work\\bio\\rlib"
ignore_cols <- c()

maf.thresh <- 0.01
phe.as.fam <- FALSE
remove.dups <- TRUE

plink.files <- file.path(base, dataset_name)
file_pattern <- sub(".bed", "", list.files(plink.files, pattern = "*.bed")[1])

source(file.path(scriptdir, "readers.R"))

pheno.desc <- read.phenotype.ordering(plink.files)

plink.data <- read.snps.plink.binary(plink.files, 
                                     file_pattern, 
                                     pheno.desc, 
                                     use.phe = phe.as.fam, 
                                     remove.dups = remove.dups,
                                     maf.threshold = maf.thresh)

SNPs.num <- data.table(plink.data$d)
data.p <- data.table(plink.data$p)

#### under development

drug <- "EMB"
stat_file <- file.path(base, "gemma_out\\Dataset_1-1stLine_bin", drug, sprintf("%s-lmm.c.assoc.txt", drug))

snps.cor <- cor(SNPs.num)
p.curr <- data.p[,drug, with = FALSE]

s <- data.table(read.table(stat_file, header = TRUE))
s <- s[ps %in% colnames(SNPs.num)]

p.threshold <- 0.01
hard.verdict <- s$p_lrt < p.threshold
