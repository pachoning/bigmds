perform_procrustes <- function(x, target, matrix_to_transform, translation = FALSE, dilation = FALSE) {

  n_row <- nrow(x)
  n_col <- ncol(x)

  if (n_row != nrow(target)) {
    stop("x and target do not have same number of rows.\n")
  }

  if (n_col != ncol(target)) {
    stop("x and target do not have same number of columns.\n")
  }

  if (n_col != ncol(matrix_to_transform)) {
    stop("x and matrix_to_transform do not have same number of columns.\n")
  }

  diag_matrix <- diag(n_row)
  if (translation) {
    diag_matrix <- diag(n_row) - 1/n_row * matrix(1, n_row, n_row)
  }

  matrix_prod <- t(target) %*% diag_matrix %*% x
  svd_results <- svd(matrix_prod)
  rotation_matrix <- svd_results$v %*% t(svd_results$u)

  dilation_factor <- 1
  if (dilation) {
    mat1 <- t(target) %*% diag_matrix %*% x %*% rotation_matrix
    mat2 <- t(x) %*% diag_matrix %*% x
    num <- 0
    denom <- 0

    for (i in 1:n_col) {
      num <- num + mat1[i, i]
      denom <- denom + mat2[i, i]
    }

    dilation_factor <- num/denom
  }

  translation_vector <- matrix(0, n_col, 1)
  if (translation) {
    translation_vector <- 1/n_row * t(target - dilation_factor * x %*% rotation_matrix) %*% matrix(1, n_row, 1)
  }

  ones_vector <- matrix(data = 1, nrow = nrow(matrix_to_transform), ncol = 1)
  translation_matrix <- ones_vector %*% t(translation_vector)

  return(dilation_factor * matrix_to_transform %*% rotation_matrix + translation_matrix)
}


classical_mds <- function(x, k, dist_fn = stats::dist, return_distance_matrix = FALSE, ...) {

  if (!is.function(dist_fn)) {
    stop("dist_fn must be a function")
  }

  mds <- list()
  dist_matrix <- dist_fn(x, ...)
  mds_result <- stats::cmdscale(d = dist_matrix, k = k, eig = TRUE)

  mds$points <- mds_result$points
  mds$eigen <- mds_result$eig[1:k]

  if (return_distance_matrix) {
    mds$distance <- as.matrix(dist_matrix)
  }

  return(mds)
}
