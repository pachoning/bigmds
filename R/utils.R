get_procrustes_parameters <- function(x, target, translation = FALSE) {

  n_row <- nrow(x)
  n_col <- ncol(x)

  if (translation) {
    c.target <- scale(target, center = TRUE, scale = FALSE)
    m.target <- attr(c.target, "scaled:center")

    c.x <- scale(x, center = TRUE, scale = FALSE)
    m.x <- attr(c.x, "scaled:center")

    matrix_prod <- t(c.target) %*% x
    svd_results <- svd(matrix_prod)
    rotation_matrix <- svd_results$v %*% t(svd_results$u)

    translation_vector <- m.target - t(rotation_matrix) %*% m.x

  } else {
    matrix_prod <- t(target) %*% x
    svd_results <- svd(matrix_prod)
    rotation_matrix <- svd_results$v %*% t(svd_results$u)
    translation_vector <- matrix(data = 0, nrow = n_col, ncol = 1)
  }

  return(list(rotation_matrix = rotation_matrix, translation_vector = translation_vector))
}

perform_procrustes <- function(x, target, matrix_to_transform, translation = FALSE) {

  n_row <- nrow(x)
  n_col <- ncol(x)

  if (n_row != nrow(target)) {
    stop("\"x\" and \"target\" do not have the same number of rows")
  }

  if (n_col != ncol(target)) {
    stop("\"x\" and \"target\" do not have the same number of columns")
  }

  if (n_col != ncol(matrix_to_transform)) {
    stop("\"x\" and \"matrix_to_transform\" do not have the same number of columns")
  }

  procrustes_parameters <- get_procrustes_parameters(x = x, target = target, translation = translation)

  ones_vector <- matrix(data = 1, nrow = nrow(matrix_to_transform), ncol = 1)
  translation_matrix <- matrix(data = 0, nrow = nrow(matrix_to_transform), ncol = ncol(matrix_to_transform))
  translation_matrix <- translation_matrix + ones_vector %*% t(procrustes_parameters$translation_vector)

  return(matrix_to_transform %*% procrustes_parameters$rotation_matrix + translation_matrix)
}

classical_mds <- function(x, k, dist_fn, return_distance_matrix = FALSE, ...) {

  if (!is.function(dist_fn)) {
    stop("\"dist_fn\" must be a function")
  } else if (any(is.na(x))) {
    stop("There are some \"NA\" values in the data. Please remove them")
  }

  mds <- list()
  dist_matrix <- dist_fn(x, ...)
  mds_result <- stats::cmdscale(d = dist_matrix, k = k, eig = TRUE)

  mds$points <- mds_result$points
  mds$eigen <- mds_result$eig[1:k]
  mds$GOF <- mds_result$GOF

  if (return_distance_matrix) {
    mds$distance <- as.matrix(dist_matrix)
  }

  return(mds)
}
