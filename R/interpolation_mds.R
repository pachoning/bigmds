#'@title Interpolation MDS
#'
#'@description Given that the size of the data set is too large, this algorithm 
#'consists of taking a random sample from it of size 
#'\code{l} \eqn{\leq \bar{l}}, being \eqn{\bar{l}} the limit size for which 
#'classical MDS is applicable, to perform classical MDS to it, and to extend the 
#'obtained results to the rest of the data set by using Gower's 
#'interpolation formula, which allows to add a new set of points 
#'to an existing MDS configuration. 
#'
#'@details *Gower's interpolation formula* is the central piece of this algorithm 
#'since it allows to add a new set of points to an existing MDS configuration 
#'so that the new one has the same coordinate system. 
#'
#'Given the matrix \code{x} with \eqn{n} points (rows) and 
#'and \eqn{k} variables (columns), a first data subsets (based on a random sample) 
#'of size \code{l} is taken and it is used to compute a MDS configuration. 
#'
#'The remaining part of \code{x} is divided into \eqn{p=({n}-}\code{l})\code{/l} 
#'data subsets (randomly). For every data subset, it is obtained a MDS 
#'configuration by means of *Gower's interpolation formula* and the first 
#'MDS configuration obtained previously. Every MDS configuration is appended 
#'to the existing one so that, at the end of the process, a global MDS 
#'configuration for \code{x} is obtained.
#'
#'This method is similar to [landmark_mds()] and [reduced_mds()].
#'
#'@param x A matrix with \eqn{n} individuals (rows) and \eqn{k} variables (columns).
#'@param l The size for which classical MDS can be computed efficiently 
#'(using `cmdscale` function). It means that if \eqn{\bar{l}} is the limit 
#'size for which classical MDS is applicable, then \code{l}\eqn{\leq \bar{l}}.
#'@param r Number of principal coordinates to be extracted.
#'@param n_cores Number of cores wanted to use to run the algorithm.
#'
#'@return Returns a list containing the following elements:
#' \describe{
#'   \item{points}{A matrix that consists of \eqn{n} individuals (rows) 
#'   and \code{r} variables (columns) corresponding to the principal coordinates. Since 
#'   we are performing a dimensionality reduction, \code{r}\eqn{<<k}}
#'   \item{eigen}{The first \code{r} largest eigenvalues:
#'   \eqn{\lambda_i, i \in  \{1, \dots, r\} }, where each \eqn{\lambda_i} is obtained
#'   from applying classical MDS to the first data subset.}
#'}
#' 
#'@references
#'Delicado P. and C. Pachón-García (2021). *Multidimensional Scaling for Big Data*.
#'\url{https://arxiv.org/abs/2007.11919}.
#'
#'Borg, I. and P. Groenen (2005). *Modern Multidimensional Scaling: Theory and Applications*. Springer.
#'
#'Gower JC. (1968). *Adding a point to vector diagrams in multivariate analysis*. Biometrika.
#'
#'@examples
#'set.seed(42)
#'x <- matrix(data = rnorm(4 * 10000), nrow = 10000) %*% diag(c(9, 4, 1, 1))
#'mds <- interpolation_mds(x = x, l = 200, r = 2, n_cores = 1)
#'head(mds$points)
#'mds$eigen
#'
#'@importFrom parallel mclapply
#'@importFrom pracma distmat
#'@importFrom stats cov
#'
#'@export
interpolation_mds <- function(x, l, r, n_cores) {
  
  n <- nrow(x)
  n_row_partition <- l
  indexes_partition <- get_partitions_for_interpolation(n = n, n_obs = n_row_partition, l = l, r = r)
  num_partitions <- length(indexes_partition)
  
  if (num_partitions <= 1) {
    # It is possible to run MDS directly
    mds <- classical_mds(x = x, r = r)
    points <- mds$points
    eigen_v <- mds$eigen/n
    list_to_return <- list(points = points, eigen = eigen_v)
  } else {
    # Get the first group 
    n_row_1 <- length(indexes_partition[[1]])
    
    # Obtain MDS for the first group
    data_1 <- x[indexes_partition[[1]], , drop = FALSE]
    mds_eig <- classical_mds(x = data_1, r = r, return_distance_matrix = TRUE)
    distance_matrix <- as.matrix(mds_eig$distance)
    X_1 <- mds_eig$points
    eigen_v <- mds_eig$eigen/nrow(X_1)
    
    # Get P matrix
    P <- get_P_matrix(n_row = n_row_1)
    Q <- -1/2 * P %*% distance_matrix^2 %*% t(P)
    q_vector <- diag(Q)
    S <- 1 / (n_row_1-1) * t(X_1) %*% X_1
    x_1__s_1__inv <- X_1 %*% solve(S)
    
    # Calculations needed to do Gower interpolation
    mds_others <- parallel::mclapply(
      indexes_partition[2:num_partitions],
      interpolation_mds_main,
      x = x,
      data_1 = data_1,
      x_1 = X_1,
      n_row_1 = n_row_1,
      q_vector = q_vector,
      x_1__s_1__inv = x_1__s_1__inv,
      mc.cores = n_cores
    )
    
    mds_points <- matrix(data = NA, nrow = n, ncol = r)
    mds_points[1:n_row_1, ] <- X_1
    mds_points[(n_row_1 + 1):n, ] <- do.call(rbind, mds_others)
    idx_all <- do.call(c, indexes_partition)
    idx_all <- order(idx_all)
    mds_points <- mds_points[idx_all, , drop = FALSE]
    mds_points <- apply(mds_points, MARGIN = 2, FUN = function(y) y - mean(y))
    mds_points <- mds_points %*% base::eigen(stats::cov(mds_points))$vectors
    list_to_return <- list(points = mds_points, eigen = eigen_v)
  }
  
  return(list_to_return)
}

get_partitions_for_interpolation <- function(n, n_obs, l, r) {
  
  if (l <= r) {
    stop("l must be greater than r")
  }
  
  if (n <= l) {
    p <- 1
  } else {
    p <- 1 + ceiling((n - l)/n_obs)
    n_last <- n - (l + (p - 2) * n_obs)
  }
  
  permutation <- sample(x = n, size = n, replace = FALSE)
  
  if (p > 1) {
    first_part <- permutation[1:l]
    middle_part <- permutation[(l + 1):(n - n_last)]
    last_part <- permutation[(n - n_last + 1):n]
    
    list_index <- split(middle_part, 1:(p - 2))
    names(list_index) <- NULL
    list_index[[p - 1]] <- list_index[[1]]
    list_index[[1]] <- first_part
    list_index[[p]] <- last_part
    
  } else {
    list_index <- list(permutation)
  }
  
  return(list_index)
}

get_P_matrix <- function(n_row) {
  
  identity_matrix <- diag(x = 1, nrow = n_row, ncol = n_row)
  one_vector <- matrix(data = 1, nrow = n_row, ncol = 1)
  P <- identity_matrix - 1/n_row * one_vector %*% t(one_vector)
  return(P)
  
}

interpolation_mds_main <- function(idx, x, data_1, x_1, n_row_1, q_vector, x_1__s_1__inv) {
  
  # Filter the matrix
  x_other <- x[idx, , drop = FALSE]
  n_row_other <- nrow(x_other)
  n_row_1 <- nrow(data_1)
  
  # Get A matrix
  A <- pracma::distmat(X = x_other, Y = data_1)
  
  # Get delta matrix
  A_sq <- A^2
  
  # One vecotr
  one_vector_other <- matrix(data = 1, nrow = n_row_other, ncol = 1)
  
  # Get coordinates
  x_2 <- 1/(2 * n_row_1) * (one_vector_other %*% t(q_vector) - A_sq) %*% x_1__s_1__inv
  
  return(x_2)
}