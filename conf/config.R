config <- list(
  common = list(
    use_cols = NULL,
    ignore_cols = NULL,
    phe_as_fam = FALSE,
    remove_dups = TRUE,
    maf_thresh = 0.01
  ),
  lasso = list(
    signif_level = 0.001,
    alpha = 0.7
  ),
  moss = list(
    alpha = 1,
    c = 0.1,
    cPrime = 1e-5,
    q = 0.1,
    replicas = 3,
    vars = 3,
    cv_k = 10
  ),
  gemma = list(
    pheno_pos = "2",
    pheno_neg = "1",
    p_threshold = 0.01,
    p_maxsel = 100,
    p_adj_m = c("bonferroni", "holm", "hochberg", "fdr", "BY")
  ),
  randforests = list(
    trees = 5000,
    vars = 200
  ),
  frn = list(
    corr_matrix = NULL
  )
)
