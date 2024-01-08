test_that("Partitions for fast MDS returns a valid partition dataset", {
  partition <- get_partitions_for_fast(n = 1000, l = 100, s_points = 5, r = 10)
  p <- length(partition)
  expect_equal(p, 20)
})

test_that("fast MDS returns a valid MDS configuration when n > l", {
  n <- 1000
  n_cols <- 10
  l <- 100
  r <- 4
  s_points <- 2*r
  diag_mat <- sqrt(diag(c(rep(15, r), rep(1, n_cols - r))))
  x <- matrix(data = rnorm(n_cols*n), nrow = n) %*% diag_mat
  cmds <- fast_mds(x = x, l = l, s_points = s_points, r = r, n_cores = 1)
  cmds_proc <- perform_procrustes(x = cmds$points,
                                  target = x[, 1:r],
                                  matrix_to_transform = cmds$points, 
                                  translation = FALSE)
  corr_vector <- sapply(1:r, function(i, x, y) cor(x[, i], y[, i]), x = x, y = cmds_proc)
  min_corr <- min(abs(corr_vector))
  expect_gt(min_corr, 0.9)
})

test_that("fast MDS returns a valid MDS configuration when n = l", {
  n <- 100
  n_cols <- 10
  l <- 100
  r <- 4
  s_points <- 2*r
  diag_mat <- sqrt(diag(c(rep(15, r), rep(1, n_cols - r))))
  x <- matrix(data = rnorm(n_cols*n), nrow = n) %*% diag_mat
  cmds <- fast_mds(x = x, l = l, s_points = s_points, r = r, n_cores = 1)
  cmds_proc <- perform_procrustes(x = cmds$points,
                                  target = x[, 1:r],
                                  matrix_to_transform = cmds$points, 
                                  translation = FALSE)
  corr_vector <- sapply(1:r, function(i, x, y) cor(x[, i], y[, i]), x = x, y = cmds_proc)
  min_corr <- min(abs(corr_vector))
  expect_gt(min_corr, 0.9)
})

test_that("fast MDS returns a valid MDS configuration when n < l", {
  n <- 90
  n_cols <- 10
  l <- 100
  r <- 4
  s_points <- 2*r
  diag_mat <- sqrt(diag(c(rep(15, r), rep(1, n_cols - r))))
  x <- matrix(data = rnorm(n_cols*n), nrow = n) %*% diag_mat
  cmds <- fast_mds(x = x, l = l, s_points = s_points, r = r)
  cmds_proc <- perform_procrustes(x = cmds$points,
                                  target = x[, 1:r],
                                  matrix_to_transform = cmds$points, 
                                  translation = FALSE)
  corr_vector <- sapply(1:r, function(i, x, y) cor(x[, i], y[, i]), x = x, y = cmds_proc)
  min_corr <- min(abs(corr_vector))
  expect_gt(min_corr, 0.9)
})
