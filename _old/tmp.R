corr.m <- cor(X)

snps <- colnames(X)
res <- data.table(snp1 = character(0), snp2 = character(0), 
                  corr = numeric(0), chisq = numeric(0))

for(i in seq(1, length(snps) - 1)) {
  for(j in seq(i + 1, length(snps))) {
    col1 <- snps[i]
    col2 <- snps[j]
    f1 <- factor(X[, col1], levels = c(0, 1))
    f2 <- factor(X[, col2], levels = c(0, 1))
    tbl <- table(f1, f2)
    test.data <- c(tbl["0", "0"], tbl["0", "1"], tbl["1", "0"], tbl["1", "1"])
    test.res <- chisq.test(test.data)
    res <- rbind(res, list(col1, col2, corr.m[col1, col2], test.res$p.value))
  }
}
