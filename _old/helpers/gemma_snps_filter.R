
scriptdir <- "D:\\work\\bio\\rlib\\common"
source(file.path(scriptdir, "config.R"))
source(lib.reader)

snp_prefix <- "MT_H37RV_BRD_V5:"

for(dataset.dir in raw.files) {
  plink.files <- dataset.dir
  file_pattern <- func.plink.filename(plink.files)
  out.file <- file.path(plink.files, "nodup_snps.txt")
  
  pheno.desc <- read.phenotype.ordering(plink.files)

  plink.data <- read.snps.plink.binary(plink.files, 
                                     file_pattern, 
                                     pheno.desc, 
                                     use.phe = phe.as.fam, 
                                     remove.dups = remove.dups,
                                     maf.threshold = maf.thresh)

  SNPs.num <- plink.data$d

  ndup.snps <- paste0(snp_prefix, colnames(SNPs.num))
  write(ndup.snps, file = out.file, sep = "\n")
}