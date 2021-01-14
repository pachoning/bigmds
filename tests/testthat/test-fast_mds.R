test_that("Fast MDS returns a valid MDS configuration", {
  x <- matrix(data = rnorm(4*1000, sd = 10), nrow = 1000)
  cmds <- fast_mds(x = x, l = 100, s = 8, k = 4)
  cmds_proc <- perform_procrustes(x = cmds$points, target = x, matrix_to_transform = cmds$points, 
                                  translation = FALSE, dilation = FALSE)
  corr_first <- cor(x[, 1], cmds_proc[, 1])
  expect_gt(corr_first, 0.9)
})
