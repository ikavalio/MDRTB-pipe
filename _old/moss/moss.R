library(genMOSS)

# Script execution can take a lot of time, so loop over all datasets is not used

scriptdir <- "D:\\work\\bio\\rlib\\common"
source(file.path(scriptdir, "config.R"))
source(lib.reader)

dataset_name <- "dataset_7"
dataset <- sub("dataset_(.+)", "\\1", dataset_name)
out.base <- home.moss

plink.files <- file.path(base, dataset_name)
file_pattern <- func.plink.filename(plink.files)
out.dir <- file.path(out.base, sprintf("%s-%s", dataset_name, file_pattern))

ignore_ds <- read.table("D:\\work\\bio\\Exclude_mutations.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

moss.alpha <- 1
moss.c <- 0.1
moss.cPrime <- 1e-5 
moss.q <- 0.1
moss.replicas <- 3
moss.vars <- 3
moss.cv.k <- 5

if (!file.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

pheno.desc <- read.phenotype.ordering(plink.files)

plink.data <- read.snps.plink.binary(plink.files, 
                                     file_pattern, 
                                     pheno.desc, 
                                     use.phe = phe.as.fam, 
                                     remove.dups = remove.dups,
                                     maf.threshold = maf.thresh)

SNPs.num <- plink.data$d
data.p <- plink.data$p

for (i in 1:ncol(data.p)) {
  drug <- colnames(data.p)[i]
  ignore_cols <- sapply(
    ignore_ds[(ignore_ds$Dataset == dataset) & (ignore_ds$Drug == drug), "Position"],
    toString
  )
  
  d <- if (length(ignore_cols) > 0) {
    SNPs.num[,setdiff(colnames(SNPs.num), ignore_cols)]
  } else {
    SNPs.num
  }
  
  
  drug.test <- data.p[,i]
  moss.data <- data.frame(cbind(d, drug.test))
  moss.data <- moss.data[!is.na(drug.test),]
  
  moss.r <- MOSS_GWAS(alpha = moss.alpha, 
                      c = moss.c, 
                      cPrime = moss.cPrime,
                      q = moss.q, 
                      replicates = moss.replicas, 
                      maxVars = moss.vars, 
                      moss.data,
                      rep(2, ncol(moss.data)), 
                      k = NULL)
  
  out.dir.ext <- file.path(out.dir, paste0(drug, extra))
  if (!file.exists(out.dir.ext)) dir.create(out.dir.ext, recursive = FALSE)
  nm <- function(t) file.path(out.dir.ext, sprintf("%s.moss.%s.csv", drug, t))
  
  write.csv(moss.r$topRegressions, nm("regr"))
  write.csv(moss.r$postIncProbs, nm("probs"))
  write.csv(moss.r$interactionModels, nm("mod"))
  write.csv(moss.r$cvMatrix, nm("cv"))
}
