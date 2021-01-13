test_that("Procrustes works", {
  x <- matrix(data = c(1, 2, 3, 4), nrow = 2)
  y <- matrix(data = c(-1, -2, -3, -4), nrow = 2)
  x_proc <- perform_procrustes(x=x, target=y, matrix_to_transform=x, translation=FALSE, dilation=FALSE)
  max_error <- max(y - x_proc)
  expect_lt(max_error, 1e-6)
})
