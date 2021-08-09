test_that("Partitions for divide-and-conquer MDS returns a valid partition dataset", {
  partition <- get_partitions_for_interpolation(n = 1000, n_obs = 100, l = 100, r = 10)
  p <- length(partition)
  expect_equal(p, 10)
})

test_that("Number of partition for interpolation MDS is 1 when n = l", {
  partition <- get_partitions_for_interpolation(n = 1000, n_obs = 1000, l = 1000, r = 10)
  p <- length(partition)
  expect_equal(p, 1)
})

test_that("Number of partition for interpolation MDS is 1 when n < l", {
  partition <- get_partitions_for_interpolation(n = 999, n_obs = 1000, l = 1000, r = 10)
  p <- length(partition)
  expect_equal(p, 1)
})

test_that("Partition for interpolation MDS fails when l < r", {
  expect_error(
    get_partitions_for_interpolation(n = 100, n_obs = 10, l = 10, r = 11),
    "\"l\" must be greater than \"r\""
  )
})

test_that("interpolation MDS returns a valid MDS configuration when n > l", {
  n <- 1000
  n_cols <- 10
  l <- 100
  r <- 4
  s_points <- 2*r
  diag_mat <- sqrt(diag(c(rep(15, r), rep(1, n_cols - r))))
  x <- matrix(data = rnorm(n_cols*n), nrow = n) %*% diag_mat
  cmds <- interpolation_mds(x = x, l = l, r = r)
  cmds_proc <- perform_procrustes(x = cmds$points,
                                  target = x[, 1:r],
                                  matrix_to_transform = cmds$points, 
                                  translation = FALSE)
  corr_vector <- sapply(1:r, function(i, x, y) cor(x[, i], y[, i]), x = x, y = cmds_proc)
  min_corr <- min(abs(corr_vector))
  expect_gt(min_corr, 0.9)
})

test_that("interpolation MDS returns a valid MDS configuration when n = l", {
  n <- 100
  n_cols <- 10
  l <- 100
  r <- 4
  s_points <- 2*r
  diag_mat <- sqrt(diag(c(rep(15, r), rep(1, n_cols - r))))
  x <- matrix(data = rnorm(n_cols*n), nrow = n) %*% diag_mat
  cmds <- interpolation_mds(x = x, l = l, r = r)
  cmds_proc <- perform_procrustes(x = cmds$points,
                                  target = x[, 1:r],
                                  matrix_to_transform = cmds$points, 
                                  translation = FALSE)
  corr_vector <- sapply(1:r, function(i, x, y) cor(x[, i], y[, i]), x = x, y = cmds_proc)
  min_corr <- min(abs(corr_vector))
  expect_gt(min_corr, 0.9)
})

test_that("interpolation MDS returns a valid MDS configuration when n < l", {
  n <- 90
  n_cols <- 10
  l <- 100
  r <- 4
  s_points <- 2*r
  diag_mat <- sqrt(diag(c(rep(15, r), rep(1, n_cols - r))))
  x <- matrix(data = rnorm(n_cols*n), nrow = n) %*% diag_mat
  cmds <- interpolation_mds(x = x, l = l, r = r)
  cmds_proc <- perform_procrustes(x = cmds$points,
                                  target = x[, 1:r],
                                  matrix_to_transform = cmds$points, 
                                  translation = FALSE)
  corr_vector <- sapply(1:r, function(i, x, y) cor(x[, i], y[, i]), x = x, y = cmds_proc)
  min_corr <- min(abs(corr_vector))
  expect_gt(min_corr, 0.9)
})
