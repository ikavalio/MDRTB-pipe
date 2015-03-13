#source("http://bioconductor.org/biocLite.R")
#biocLite(c("snpStats"))
library(snpStats)
library(randomForest)

plink.files <- "D:\\work\\bio"
maf.thresh <- 0.01
drug <- "INH"
rf.trees <- 5000
rf.vars <- 500
basedir <- plink.files

pheno.desc <- c(
  "EMB", "INH", "RIF", "RIFP", 
  "PZA", "STM", "CYCL", "ETH", 
  "PARA", "AMIK", "CAPR", "KANA", 
  "OFLO", "R1-T1", "R1-T2", "R1-T3", 
  "R2-T1", "R2-T2", "R2-T3", "R2-T4",
  "R2-T5", "TOT-T1", "TOT-T2"
)

file <- file.path(basedir, "snps.matrix")
pheno <- file.path(basedir, "phenotype.v2.csv")

data.m <- read.csv(file, header=TRUE, sep="\t", check.names = FALSE, na.strings=c("-"))
data.p <- read.csv(pheno, header=TRUE, sep=",", check.names = FALSE, na.strings=c("-"))
rownames(data.p) <- data.p$"snp_pos"

rownames(data.m) <- data.m$snp_pos
data.k <- data.m[,!(names(data.m) %in% "snp_pos")]
data.m.nona <- data.frame(
  apply(data.k, 2, 
        function(col) { 
          ifelse(is.na(col), round(median(col, na.rm = TRUE)), col) 
        }),
  check.names = FALSE)
data.m.t <- t(data.m.nona)
data.m.t <- data.m.t[!duplicated(data.m.t),]
data.m.t <- cbind(data.m.t, rowSums(data.m.t) / ncol(data.m.t) )
data.m.t <- data.m.t[data.m.t[,ncol(data.m.t)] >= maf.thresh,] # MAF threshold
data.m.t <- data.m.t[,1:ncol(data.m.t) - 1]
data.m <- data.frame(t(data.m.t), check.names = FALSE)

drugid <- which(colnames(data.p) == drug)
tM <- data.p[,c(1, drugid)]
tM <- tM[!is.na(tM[,2]),]

# remove strange classes "resistant/susceptible"
tM <- tM[tM[,2] == "resistant" | tM[,2] == "susceptible",]
# correct available factor values
tM[,2] <- as.factor(as.character(tM[,2]))

clPheno <- merge(tM, data.m, all = FALSE, by = "row.names")
X <- clPheno[,4:ncol(clPheno)]
y <- clPheno[,3]

res.rf <- randomForest(x = X,
                       y = y,
                       ntree = rf.trees,
                       mtry = rf.vars,
                       importance = TRUE
                       )

