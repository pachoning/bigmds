test_that("LMDS returns a valid MDS configuration", {
  n <- 1000
  n_cols <- 10
  l <- 100
  r <- 4
  c_points <- 2*r
  diag_mat <- sqrt(diag(c(rep(15, r), rep(1, n_cols - r))))
  x <- matrix(data = rnorm(n_cols*n), nrow = n) %*% diag_mat
  cmds <- pivot_mds(x = x, num_pivots = l, r = r)
  cmds_proc <- perform_procrustes(x = cmds$points,
                                  target = x[, 1:r],
                                  matrix_to_transform = cmds$points, 
                                  translation = FALSE)
  corr_vector <- sapply(1:r, function(i, x, y) cor(x[, i], y[, i]), x = x, y = cmds_proc)
  min_corr <- min(abs(corr_vector))
  expect_gt(min_corr, 0.9)
})