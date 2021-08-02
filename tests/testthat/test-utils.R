#test_that("Procrustes aligns matrices", {
#  x <- matrix(data = c(1, 2, 3, 4), nrow = 2)
#  y <- matrix(data = c(-1, -2, -3, -4), nrow = 2)
#  x_proc <- perform_procrustes(x=x, target=y, matrix_to_transform=x, translation=FALSE, dilation=FALSE)
#  max_error <- max(y - x_proc)
#  expect_lt(max_error, 1e-6)
#})

test_that("Classical MDS returns a valid MDS configuration", {
  x <- matrix(data = rnorm(4*100, sd = 10), nrow = 100)
  cmds <- classical_mds(x = x, k = 4, dist_fn = dist)
  cmds_proc <- perform_procrustes(x = cmds$points, target = x, matrix_to_transform = cmds$points, 
                                  translation = FALSE, dilation = FALSE)
  corr_first <- cor(x[, 1], cmds_proc[, 1])
  expect_gt(corr_first, 0.9)
})

test_that("Classical MDS fails when distance is not provided", {
  x <- matrix(data = rnorm(4*100, sd = 10), nrow = 100)
  expect_error(classical_mds(x = x, k = 4), "argument \"dist_fn\" is missing, with no default")
})

test_that("Classical MDS fails when dist_fn is not a function", {
  x <- matrix(data = rnorm(4*100, sd = 10), nrow = 100)
  expect_error(classical_mds(x = x, k = 4, dist_fn = 3), "dist_fn must be a function")
})

test_that("Classical MDS fails when there are NA values", {
  x <- matrix(data = rnorm(4*100, sd = 10), nrow = 100)
  x[1,1] <- NA
  expect_error(classical_mds(x = x, k = 4, dist_fn = dist), "There are some NA values in the data. Please remove them")
})
