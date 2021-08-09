test_that("Partitions for divide-and-conquer MDS returns a valid partition dataset", {
  partition <- get_partitions_for_divide_conquer(n = 1000, l = 100, c_points = 5, r = 10)
  p <- length(partition)
  expect_equal(p, 11)
})

test_that("Number of partition for divide-and-conquer MDS is 1 when n = l", {
  partition <- get_partitions_for_divide_conquer(n = 100, l = 100, c_points = 5, r = 10)
  p <- length(partition)
  expect_equal(p, 1)
})

test_that("Number of partition for divide-and-conquer MDS is 1 when n < l", {
  partition <- get_partitions_for_divide_conquer(n = 90, l = 100, c_points = 5, r = 10)
  p <- length(partition)
  expect_equal(p, 1)
})

test_that("Partition for divide-and-conquer MDS fails when l-c_points = 0", {
  expect_error(
    get_partitions_for_divide_conquer(n = 1000, l = 100, c_points = 100, r = 10),
    "\"l\" must be greater than \"c_points\""
    )
})

test_that("Partition for divide-and-conquer MDS fails when l-c_points < 0", {
  expect_error(
    get_partitions_for_divide_conquer(n = 1000, l = 90, c_points = 100, r = 10),
    "\"l\" must be greater than \"c_points\""
  )
})

test_that("Partition for divide-and-conquer MDS fails when l-c_points = c_points", {
  expect_error(
    get_partitions_for_divide_conquer(n = 1000, l = 90, c_points = 45, r = 10),
    "\"l-c_points\" must be greater than \"c_points\""
  )
})

test_that("Partition for divide-and-conquer MDS fails when l-c_points = r", {
  expect_error(
    get_partitions_for_divide_conquer(n = 1000, l = 100, c_points = 10, r = 90),
    "\"l-c_points\" must be greater than \"r\""
  )
})

test_that("Partition for divide-and-conquer MDS fails when l-c_points < r", {
  expect_error(
    get_partitions_for_divide_conquer(n = 1000, l = 100, c_points = 11, r = 90),
    "\"l-c_points\" must be greater than \"r\""
  )
})

test_that("divide-and-conquer MDS returns a valid MDS configuration when n > l", {
  n <- 1000
  n_cols <- 10
  l <- 100
  r <- 4
  c_points <- 2*r
  diag_mat <- sqrt(diag(c(rep(15, r), rep(1, n_cols - r))))
  x <- matrix(data = rnorm(n_cols*n), nrow = n) %*% diag_mat
  cmds <- divide_conquer_mds(x = x, l = l, c_points = c_points, r = r)
  cmds_proc <- perform_procrustes(x = cmds$points,
                                  target = x[, 1:r],
                                  matrix_to_transform = cmds$points, 
                                  translation = FALSE)
  corr_vector <- sapply(1:r, function(i, x, y) cor(x[, i], y[, i]), x = x, y = cmds_proc)
  min_corr <- min(abs(corr_vector))
  expect_gt(min_corr, 0.9)
})

test_that("divide-and-conquer MDS returns a valid MDS configuration when n = l", {
  n <- 100
  n_cols <- 10
  l <- 100
  r <- 4
  c_points <- 2*r
  diag_mat <- sqrt(diag(c(rep(15, r), rep(1, n_cols - r))))
  x <- matrix(data = rnorm(n_cols*n), nrow = n) %*% diag_mat
  cmds <- divide_conquer_mds(x = x, l = l, c_points = c_points, r = r)
  cmds_proc <- perform_procrustes(x = cmds$points,
                                  target = x[, 1:r],
                                  matrix_to_transform = cmds$points, 
                                  translation = FALSE)
  corr_vector <- sapply(1:r, function(i, x, y) cor(x[, i], y[, i]), x = x, y = cmds_proc)
  min_corr <- min(abs(corr_vector))
  expect_gt(min_corr, 0.9)
})

test_that("divide-and-conquer MDS returns a valid MDS configuration when n < l", {
  n <- 90
  n_cols <- 10
  l <- 100
  r <- 4
  c_points <- 2*r
  diag_mat <- sqrt(diag(c(rep(15, r), rep(1, n_cols - r))))
  x <- matrix(data = rnorm(n_cols*n), nrow = n) %*% diag_mat
  cmds <- divide_conquer_mds(x = x, l = l, c_points = c_points, r = r)
  cmds_proc <- perform_procrustes(x = cmds$points,
                                  target = x[, 1:r],
                                  matrix_to_transform = cmds$points, 
                                  translation = FALSE)
  corr_vector <- sapply(1:r, function(i, x, y) cor(x[, i], y[, i]), x = x, y = cmds_proc)
  min_corr <- min(abs(corr_vector))
  expect_gt(min_corr, 0.9)
})
