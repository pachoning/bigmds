#'@title Reduced MDS
#'
#'@description A data subset is selected and classical MDS is performed on it to obtain the 
#'corresponding low dimensional configuration.Then the reaming points are projected 
#'onto this initial configuration.
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
#'data subsets (randomly). For every data point, it is obtained a MDS 
#'configuration by means of *Gower's interpolation formula* and the first 
#'MDS configuration obtained previously. Every MDS configuration is appended 
#'to the existing one so that, at the end of the process, a global MDS 
#'configuration for \code{x} is obtained.
#'
#'#'This method is similar to [landmark_mds()] and [interpolation_mds()].
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
#'Paradis E. (2018). *Multidimensional Scaling With Very Large Datasets*. Journal of Computational and Graphical Statistics.
#'
#'Borg, I. and P. Groenen (2005). *Modern Multidimensional Scaling: Theory and Applications*. Springer.
#'
#'Gower JC. (1968). *Adding a point to vector diagrams in multivariate analysis*. Biometrika.
#'
#'@examples
#'set.seed(42)
#'x <- matrix(data = rnorm(4 * 10000), nrow = 10000) %*% diag(c(9, 4, 1, 1))
#'mds <- reduced_mds(x = x, l = 200, r = 2, n_cores = 1)
#'head(mds$points)
#'mds$eigen
#'
#'@importFrom parallel mclapply
#'@importFrom stats cov
#'
#'@export
reduced_mds <- function(x, l, r, n_cores) {
  n <- nrow(x)
  
  # Take a random sample
  idx <- sample.int(n, size = l)
  
  # Get MDS for the random sample
  x_initial <- x[idx, , drop = FALSE]
  x_initial_mds <- classical_mds(x = x_initial, r = r)
  mds_initial <- x_initial_mds$points
  
  # Prepare to use Gower's formula
  t_mds_initial <- t(mds_initial)
  Q <- mds_initial %*% t_mds_initial
  q <- diag(Q)
  A <- 0.5 * solve(t_mds_initial %*% mds_initial) %*% t_mds_initial
  
  mds_config <- numeric(n * r)
  dim(mds_config) <- c(n, r)
  
  # Get config for the remaining points
  idx_remaining <- (1:n)[-idx]
  mds_remaining_list <- parallel::mclapply(
    idx_remaining,
    partial_reduced_mds,
    x = x,
    l = l,
    x_initial = x_initial,
    A = A,
    q = q,
    mc.cores = n_cores
  )
  mds_remaining <- do.call(rbind, mds_remaining_list)
  # Populate final mds with initial mds
  mds_config[idx, ] <- mds_initial
  # Populate final mds with remaining mds
  mds_config[idx_remaining, ] <- mds_remaining
  # Normalise
  mds_config <- apply(mds_config, MARGIN = 2, FUN = function(y) y - mean(y))
  mds_config <- mds_config %*% base::eigen(stats::cov(mds_config))$vectors
  # Return final object
  conf <- list()
  conf$points <- mds_config
  conf$eigen <- x_initial_mds$eigen/l
  return(conf)
}

partial_reduced_mds <- function(i, x, l, x_initial, A, q) {
  current_point <-  matrix(x[i, , drop = FALSE], l, ncol(x), byrow = TRUE)
  dist_points <- rowSums((x_initial - current_point)^2)
  result <- t(A %*% (q - dist_points))
  return(result)
}
