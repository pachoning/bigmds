#'@title MDS based on Gower interpolation formula
#'@description Performs *Multidimensional Scaling* for big datasets using Gower interpolation formula. This method can 
#'compute a MDS configuration even when the dataset is so large that classical MDS methods (`cmdscale`) can not be run 
#'due to computational problems.
#'@details *Gower interpolation formula* is the central piece of this algorithm since it allows to add a new set of 
#'points to an existing MDS configuration. 
#'
#'Given the matrix \code{x} with n individuals (rows) and q variables (columns), a random sample of \code{l} individuals
#'is taken and used to compute a MDS configuration. This configuration will be used in order to obtain a MDS
#'configuraton for the entire matrix \code{x}. 
#'
#'In order to obtain a MDS configuration for the remaining part of \code{x}, it is divided into p=\code{(n-l)/l} 
#'submatrices. For every partition, Gower interpolation formula is used to compute and append  a MDS configuration to the
#'exisitng ones.
#'@param x A matrix with n individuals (rows) and q variables (columns).
#'@param l The largest value which allows classical MDS to be computed efficiently, i.e, the largest value which makes 
#'`cmdscale()` be run without any computational issues.
#'@param k Number of principal coordinates to be extracted.
#'@param dist_fn Distance function to be used for obtaining a MDS configuration.
#'@return Returns a list containing the following elements:
#' \describe{
#'   \item{points}{A matrix that consists of n individuals (rows) and \code{k} variables (columns) corresponding to the 
#'   MDS coordinates.}
#'   \item{eigen}{The first \code{k} eigenvalues.}
#' }
#'@examples
#'set.seed(42)
#'x <- matrix(data = rnorm(4*10000), nrow = 10000) %*% diag(c(15, 10, 1, 1))
#'mds <- gower_interpolation_mds(x = x, l = 200, k = 2, dist_fn = stats::dist)
#'head(cbind(mds$points, x[, 1:2]))
#'var(x)
#'var(mds$points)
#'@references
#'Gower, J.C. and D.J, Hand (1995). *Biplots*. Volume 54. CRC Press.
#'
#'Borg I and P. Groenen (1997). *Modern Multidimensional Scaling: Theory and Applications*. New York: Springer. pp. 340-342.
#' @export
gower_interpolation_mds <- function(x, l, k, dist_fn = stats::dist) {

  nrow_x <- nrow(x)
  p <- ceiling(nrow_x / l)

  if (p < 1) {
    p <- 1
  } 

  if (p > 1) {

    initial_row_names <- row.names(x)
    row.names(x) <- 1:nrow(x)
    
    # Do MDS with the first group and then use the Gower interpolation formula
    sample_distribution <- sample(x = p, size = nrow_x, replace = TRUE)

    # Get the first group 
    ind_1 <- which(sample_distribution == 1)
    n_1 <- length(ind_1)

    # Do MDS with the first group
    submatrix_data <- x[ind_1, ,drop = FALSE]
    mds_eig <- classical_mds(x = submatrix_data, k = k, dist_fn = dist_fn, return_distance_matrix = TRUE)
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
      row.names(MDS_i_group) <- row.names(submatrix_data)
      cum_mds <- rbind(cum_mds, MDS_i_group)
    }

    cum_mds <- cum_mds[order(as.numeric(row.names(cum_mds))), , drop = FALSE]
    row.names(cum_mds) <- initial_row_names
    row.names(x) <- initial_row_names

  } else {
    # It is possible to run MDS directly
    mds_eig <- classical_mds(x = x, k = k, dist_fn = dist_fn, return_distance_matrix = TRUE)
    distance_matrix <- mds_eig$distance
    cum_mds <- mds_eig$points
    eigen <- mds_eig$eig / nrow_x
  }

  return(list(points = cum_mds, eigen = eigen))
}
