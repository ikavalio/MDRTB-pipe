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
  )
)
