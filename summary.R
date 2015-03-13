
library(data.table)

drug <- "TOT-T2"

c.p <- "Assoc."
c.m <- "Weak Assoc."
c.n <- "No Assoc."

tests.output <- file.path("D:\\work", "bio", "rawdata", fsep = "\\")
snps.source <- file.path("D:\\work", "bio", "snps.matrix", fsep = "\\")

snps.m <- read.csv(snps.source, 
                   header = TRUE, 
                   sep = "\t", 
                   check.names = FALSE,
                   na.strings = c("-"))
snps.l <- ncol(snps.m) - 1

drug.data <- file.path(tests.output, drug)

res.lasso <- read.csv(file.path(drug.data, sprintf("%s.lm.csv", drug)), header = TRUE,
                      sep = ",", check.names = FALSE)
res.gemma <- read.csv(file.path(drug.data, sprintf("%s.gemma.csv", drug)), header = TRUE,
                      sep = "\t", check.names = FALSE)
res.bf <- read.csv(file.path(drug.data, sprintf("%s.bf.csv", drug)), header = TRUE,
                   sep = "\t", check.names = FALSE, skip = 1)
res.moss <- read.csv(file.path(drug.data, sprintf("%s.moss.csv", drug)), header = TRUE,
                     sep = ",", check.names = FALSE)

gemma.fdr <- p.adjust(res.gemma$p_lrt, c("fdr"))

r <- data.table(pos = tail(colnames(snps.m), snps.l))

r.init <- rep(c.n, snps.l)
r[, lasso := r.init ]
r[, gemma := r.init]
r[, bf := r.init]
r[, moss := r.init]

r[pos %in% res.lasso$bp, "lasso"] <- c.p
r[pos %in% res.gemma[res.gemma$p_lrt < 0.01,]$ps, "gemma"] <- c.m
r[pos %in% res.gemma[gemma.fdr < 0.01,]$ps, "gemma"] <- c.p
r[pos %in% res.bf[res.bf$bf > 0.5,]$pos, "bf"] <- c.m
r[pos %in% res.bf[res.bf$bf > 1,]$pos, "bf"] <- c.p
r[pos %in% res.moss$variable, "moss"] <- c.m
r[pos %in% res.moss[res.moss$postIncProb > 0.2,]$variable, "moss"] <- c.p

r.v <- r[lasso != c.n | gemma != c.n | bf != c.n | moss != c.n]

write.csv(r.v, file.path(drug.data, sprintf("%s.summary.csv", drug)))
