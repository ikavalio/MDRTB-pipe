library(data.table)

scriptdir <- "D:\\work\\bio\\rlib\\common"
source(file.path(scriptdir, "config.R"))
source(lib.reader)

extra <- ""
ignore_cols <- c()

# override defaults
maf.thresh <- NA
remove.dups <- FALSE

for(dataset.dir in raw.files) {
  plink.files <- dataset.dir
  file_pattern <- func.plink.filename(plink.files)
  out.file <- file.path(plink.files, "clusters.csv")
  
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
}


