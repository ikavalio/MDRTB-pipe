library(data.table)

scriptdir <- "D:\\work\\bio\\rlib\\common"
source(file.path(scriptdir, "config.R"))
source(lib.reader)

dataset_name <- "dataset_1"

plink.files <- file.path(base, dataset_name)
file_pattern <- func.plink.filename(plink.files)

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
lmbda <- 0.5
cor.threshold <- 0.6
p.threshold <- 0.01

# E1(x) = |x - p|
# E2(x, y) = l * |cor(x, Y)| * I(x != y)

lphi1 <- function(x, p) abs(rep(x, length(q)) - p) # p - probability of relevance
lphi2 <- function(r, lmbda) 2 * lmbda * abs(r) 

stat_file <- file.path(base, "gemma_out\\Dataset_1-1stLine_bin", drug, sprintf("%s-lmm.c.assoc.txt", drug))

snps.cor <- cor(SNPs.num)
snps.cor[snps.cor < cor.threshold] <- 0
p.curr <- data.p[,drug, with = FALSE]

s <- data.table(read.table(stat_file, header = TRUE))
s <- s[ps %in% colnames(SNPs.num)]

pHard <- s$p_lrt < p.threshold # not significant snps
pExp <- 1 - s$p_lrt # experimental
s.obs <- pHard

E1 <- lphi1(1, s.obs) # log of energy for X = 1
E0 <- lphi1(0, s.obs) # log of energy for X = 0
E2 <- lphi2(snps.cor, lmbda)

adjM.rows <- c(rownames(snps.cor), "s", "t")
adjM.cols <- c(rownames(snps.cor), "s", "t")
adjM <- matrix(rep(0, length(adjM.rows) * length(adjM.cols)), 
               nrow = length(adjM.rows), ncol = length(adjM.cols), 
               dimnames = list(adjM.rows, adjM.cols))

energComparison <- E0 < E1
energValue <- abs(E1 - E0)
adjM["s",] <- c(ifelse(energComparison, energValue, 0), 0, 0) # last 2 els are "s" and "t"
adjM[,"t"] <- c(ifelse(energComparison, 0, energValue), 0, 0) # last 2 els are "s" and "t"
#E2[lower.tri(E2, diag = TRUE)] <- 0
diag(E2) <- 0
adjM[1:ncol(snps.cor), 1:ncol(snps.cor)] <- E2

library(igraph)
frn <- graph.adjacency(adjM, weighted = TRUE, mode = "directed")
#cut <- stMincuts(frn, source = "s", target = "t")
E(frn)$capacity <- E(frn)$weight
flow <- graph.maxflow(frn, "s", "t")
E(frn)[flow$cut]
V(frn)[flow$partition2]
