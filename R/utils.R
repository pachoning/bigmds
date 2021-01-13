perform_procrustes <- function(x, target, matrix_to_transform, translation, dilation) {

  procrustes_result <- MCMCpack::procrustes(X = x, Xstar = target, translation = translation, dilation = dilation)

  if (translation) {
    trans <- procrustes_result$tt
  } else {
    trans <- matrix(data = 0, nrow = ncol(x), ncol = 1)
  }

  ones_vector <- matrix(data = 1, nrow = nrow(matrix_to_transform), ncol = 1)
  translation_matrix <- ones_vector %*% t(trans)

  if (dilation) {
    dilation_factor <- procrustes_result$s
  } else {
    dilation_factor <- 1
  }

  return(dilation_factor * matrix_to_transform %*% procrustes_result$R + translation_matrix)
}

classical_mds <- function(x, k, return_distance_matrix = FALSE) {

  mds <- list()
  dist_matrix <- stats::dist(x = x)
  mds_result <- stats::cmdscale(d = dist_matrix, k = k, eig = TRUE)

  mds$points <- mds_result$points
  mds$eigen <- mds_result$eig[1:k]

  if (return_distance_matrix) {
    mds$distance <- as.matrix(dist_matrix)
  }

  return(mds)
}
