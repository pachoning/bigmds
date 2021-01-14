#'@title MDS based on Gower interpolation formula
#'@description Performs Multidimensional Scaling based on Gower interpolation formula.
#'@param x Data matrix.
#'@param l The highest value where classical MDS can be computed efficiently.
#'@param k Number of principal coordinates.
#'@return Returns MDS based on Fast MDS algorithm as well as the first k eigenvalues.
#' \describe{
#'   \item{points}{MDS}
#'   \item{eigen}{eigenvalues}
#' }
#' @export
#' @examples
#' x <- matrix(data = rnorm(4*10000, sd = 10), nrow = 10000)
#' cmds <- gower_interpolation_mds(x = x, l = 100, k = 2)
#' head(cmds$points)
#' cmds$eigen
#' @seealso
#' \url{https://arxiv.org/abs/2007.11919}
gower_interpolation_mds <- function(x, l, k) {

  nrow_x <- nrow(x)
  p <- ceiling(nrow_x / l)

  if (p < 1) {
    p <- 1
  } 

  if (p > 1) {
    # Do MDS with the first group and then use the Gower interpolation formula
    sample_distribution <- sort(sample(x = p, size = nrow_x, replace = TRUE))
    sample_distribution <- sort(sample_distribution)

    # Get the first group 
    ind_1 <- which(sample_distribution == 1)
    n_1 <- length(ind_1)

    # Do MDS with the first group
    submatrix_data <- x[ind_1, ,drop = FALSE]
    mds_eig <- classical_mds(x = submatrix_data, k = k, return_distance_matrix = TRUE)
    distance_matrix <- mds_eig$distance

    M <- mds_eig$points
    eigen <- mds_eig$eig / nrow(M)
    cum_mds <- M

    # Calculations needed to do Gower interpolation
    delta_matrix <- distance_matrix^2
    In <- diag(n_1)
    ones_vector <- rep(1, n_1)
    J <- In - 1 / n_1 * ones_vector %*% t(ones_vector)
    G <- -1 / 2 * J %*% delta_matrix %*% t(J) 
    g_vector <- diag(G)
    S <- 1 / (nrow(M)-1) * t(M) %*% M
    S_inv <- solve(S)

    # For the rest of the groups, do the interpolation
    for (i_group in 2:p) {
      # Filtering the data
      ind_i_group <- which(sample_distribution == i_group)
      submatrix_data <- x[ind_i_group, ,drop = FALSE]

      # A matrix
      distance_matrix_filter <- pdist::pdist(
        X = submatrix_data,
        Y = x[ind_1, ,drop = FALSE]
      )

      distance_matrix_filter <- as.matrix(distance_matrix_filter)
      A <- distance_matrix_filter^2
      ones_vector <- rep(1, length(ind_i_group))
      MDS_i_group <- 1 / (2 * n_1) * (ones_vector %*% t(g_vector) - A) %*% M %*% S_inv
      cum_mds <- rbind(cum_mds, MDS_i_group)
    }
  } else {
    # It is possible to run MDS directly
    mds_eig <- classical_mds(x = x, k = k, return_distance_matrix = TRUE)
    distance_matrix <- mds_eig$distance
    cum_mds <- mds_eig$points
    eigen <- mds_eig$eig / nrow_x
  }

  return(list(points = cum_mds, eigen = eigen))
}
