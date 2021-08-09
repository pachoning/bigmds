test_that("Procrustes parameters are well recovered when translation vector is FALSE", {
  n_rows <- 100
  n_cols <- 2
  x <- matrix(data = rnorm(n_rows*n_cols, sd = 10), nrow = n_rows)
  rot_mat <- matrix(data = c(cos(45), -sin(45), sin(45), cos(45)), ncol = n_cols, nrow = n_cols)
  y <- x %*% rot_mat
  proc <- get_procrustes_parameters(x = x, target = y, translation = FALSE)
  dif_rot <- max(abs(proc$rotation_matrix - rot_mat))
  max_trans <- max(abs(proc$translation_vector))
  max_error <- max(c(dif_rot, max_trans))
  expect_lt(max_error, 1e-6)
})

test_that("Procrustes parameters are well recovered when translation vector is TRUE", {
  n_rows <- 100
  n_cols <- 2
  x <- matrix(data = rnorm(n_rows*n_cols, sd = 10), nrow = n_rows)
  rot_mat <- matrix(data = c(cos(45), -sin(45), sin(45), cos(45)), ncol = n_cols, nrow = n_cols)
  trans_vector <- matrix(data = c(1, -1), nrow = n_cols, ncol = 1)
  ones_vector <- matrix(data = 1, nrow = n_rows, ncol = 1)
  trans_matrix <- ones_vector %*% t(trans_vector)
  y <- x %*% rot_mat + trans_matrix
  proc <- get_procrustes_parameters(x = x, target = y, translation = TRUE)
  dif_rot <- max(abs(proc$rotation_matrix - rot_mat))
  max_trans <- max(abs(proc$translation_vector - trans_vector))
  max_error <- max(c(dif_rot, max_trans))
  expect_lt(max_error, 1e-6)
})

test_that("Procrustes transforms the matrix correctly when translation vector is FALSE", {
  n_rows <- 100
  n_cols <- 2
  x <- matrix(data = rnorm(n_rows*n_cols, sd = 10), nrow = n_rows)
  rot_mat <- matrix(data = c(cos(45), -sin(45), sin(45), cos(45)), ncol = n_cols, nrow = n_cols)
  y <- x %*% rot_mat
  proc <- perform_procrustes(x = x, target = y, matrix_to_transform = x, translation = FALSE)
  max_error <- max(abs(proc - y))
  expect_lt(max_error, 1e-6)
})

test_that("Procrustes transforms the matrix correctly when translation vector is TRUE", {
  n_rows <- 100
  n_cols <- 2
  x <- matrix(data = rnorm(n_rows*n_cols, sd = 10), nrow = n_rows)
  rot_mat <- matrix(data = c(cos(45), -sin(45), sin(45), cos(45)), ncol = n_cols, nrow = n_cols)
  trans_vector <- matrix(data = c(1, -1), nrow = n_cols, ncol = 1)
  ones_vector <- matrix(data = 1, nrow = n_rows, ncol = 1)
  trans_matrix <- ones_vector %*% t(trans_vector)
  y <- x %*% rot_mat + trans_matrix
  proc <- perform_procrustes(x = x, target = y, matrix_to_transform = x, translation = TRUE)
  max_error <- max(abs(proc - y))
  expect_lt(max_error, 1e-6)
})

test_that("Procrustes fails when number of rows does not match (x, target)", {
  x <- matrix(data = c(1, 2, 3, 4, 5, 6), nrow = 3)
  y <- matrix(data = c(-1, -2, -3, -4), nrow = 2)
  expect_error(
    perform_procrustes(x = x, target = y, matrix_to_transform = x, translation = FALSE),
    "\"x\" and \"target\" do not have the same number of rows"
  )
})

test_that("Procrustes fails when number of columns does not match (x, target)", {
  x <- matrix(data = c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
  y <- matrix(data = c(-1, -2, -3, -4), nrow = 2)
  expect_error(
    perform_procrustes(x = x, target = y, matrix_to_transform = x, translation = FALSE),
    "\"x\" and \"target\" do not have the same number of columns"
  )
})

test_that("Procrustes fails when number of columns does not match (x, matrix_to_transform)", {
  x <- matrix(data = c(1, 2, 3, 4, 5, 6), nrow = 3)
  x_t <- t(x)
  y <- matrix(data = c(-1, -2, -3, -4, -5, -6), nrow = 3)
  expect_error(
    perform_procrustes(x = x, target = y, matrix_to_transform = x_t, translation = FALSE),
    "\"x\" and \"matrix_to_transform\" do not have the same number of columns"
  )
})

test_that("Classical MDS returns a valid MDS configuration", {
  x <- matrix(data = rnorm(4*100, sd = 10), nrow = 100)
  cmds <- classical_mds(x = x, k = 4, dist_fn = dist)
  cmds_proc <- perform_procrustes(x = cmds$points, target = x, matrix_to_transform = cmds$points, 
                                  translation = FALSE)
  corr_first <- cor(x[, 1], cmds_proc[, 1])
  expect_gt(corr_first, 0.9)
})

test_that("Classical MDS fails when distance is not provided", {
  x <- matrix(data = rnorm(4*100, sd = 10), nrow = 100)
  expect_error(classical_mds(x = x, k = 4), "argument \"dist_fn\" is missing, with no default")
})

test_that("Classical MDS fails when dist_fn is not a function", {
  x <- matrix(data = rnorm(4*100, sd = 10), nrow = 100)
  expect_error(classical_mds(x = x, k = 4, dist_fn = 3), "\"dist_fn\" must be a function")
})

test_that("Classical MDS fails when there are NA values", {
  x <- matrix(data = rnorm(4*100, sd = 10), nrow = 100)
  x[1,1] <- NA
  expect_error(classical_mds(x = x, k = 4, dist_fn = dist), "There are some \"NA\" values in the data. Please remove them")
})
