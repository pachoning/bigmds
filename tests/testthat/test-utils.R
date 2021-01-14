test_that("Procrustes aligns matrices", {
  x <- matrix(data = c(1, 2, 3, 4), nrow = 2)
  y <- matrix(data = c(-1, -2, -3, -4), nrow = 2)
  x_proc <- perform_procrustes(x=x, target=y, matrix_to_transform=x, translation=FALSE, dilation=FALSE)
  max_error <- max(y - x_proc)
  expect_lt(max_error, 1e-6)
})

test_that("Classical MDS returns a valid MDS configuration", {
  x <- matrix(data = rnorm(4*100, sd = 10), nrow = 100)
  cmds <- classical_mds(x = x, k = 4)
  cmds_proc <- perform_procrustes(x = cmds$points, target = x, matrix_to_transform = cmds$points, 
                                  translation = FALSE, dilation = FALSE)
  corr_first <- cor(x[, 1], cmds_proc[, 1])
  expect_gt(corr_first, 0.9)
})
