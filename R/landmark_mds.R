#'@title Landmark MDS
#'
#'@description Landmark MDS (LMDS) algorithm applies first classical MDS to a 
#'subset of the data (*landmark points*) and then the remaining individuals are 
#'projected onto the landmark low dimensional configuration using a 
#'distance-based triangulation procedure.
#'
#'@details LMDS applies first classical MDS to a subset of the data (*landmark points*). Then,
#'it uses a distance-based triangulation procedure to project the non-landmark individuals. This 
#'distance-based triangulation procedure coincides with  *Gower's interpolation formula*. 
#'
#'This method is similar to [interpolation_mds()] and [reduced_mds()].
#'
#'@param x A matrix with \eqn{n} points (rows) and \eqn{k} variables (columns).
#'@param num_landmarks Number of landmark points to obtain an initial MDS configuration. It is
#'equivalent to \code{l} parameter used in [interpolation_mds()], [divide_conquer_mds()] and
#'[fast_mds()]. Therefore, it is the size for which classical MDS can be computed efficiently 
#'(using `cmdscale` function). It means that if \eqn{\bar{l}} is the limit 
#'size for which classical MDS is applicable, then \code{l}\eqn{\leq \bar{l}}.
#'@param r Number of principal coordinates to be extracted.
#'
#'@return Returns a list containing the following elements:
#' \describe{
#'   \item{points}{A matrix that consists of \eqn{n} points (rows) 
#'   and \code{r} variables (columns) corresponding to the principal coordinates. Since 
#'   a dimensionality reduction is performed, \code{r}\eqn{<<k}}
#'   \item{eigen}{The first \code{r} largest eigenvalues:
#'   \eqn{\lambda_i, i \in  \{1, \dots, r\} }, where each \eqn{\lambda_i} is obtained
#'   from applying classical MDS to the first data subset.}
#'}
#'
#'@examples
#'set.seed(42)
#'x <- matrix(data = rnorm(4 * 10000), nrow = 10000) %*% diag(c(9, 4, 1, 1))
#'mds <- landmark_mds(x = x, num_landmarks = 200, r = 2)
#'head(mds$points)
#'mds$eigen
#'
#'@references
#'Delicado P. and C. Pachón-García (2021). *Multidimensional Scaling for Big Data*.
#'\url{https://arxiv.org/abs/2007.11919}.
#'
#'Borg, I. and P. Groenen (2005). *Modern Multidimensional Scaling: Theory and Applications*. Springer.
#'
#'De Silva V. and JB. Tenenbaum (2004). *Sparse multidimensional scaling using landmark points*. Technical Report, Stanford University.
#'
#'Gower JC. (1968). *Adding a point to vector diagrams in multivariate analysis*. Biometrika.
#'
#'@importFrom pracma distmat
#'@importFrom svd trlan.eigen
#'@importFrom stats cov
#'
#'@export
landmark_mds <- function (x, num_landmarks, r) {
  # Select landmark indexes and compute the distance
  idx_lm <- select_landmarks(x = x, num_landmarks = num_landmarks)
  
  # Select landmark points
  x_lm <- x[idx_lm, , drop = FALSE]
  
  # Compute distance with respect those landmark points
  dist_to_lm <- pracma::distmat(X = x_lm, Y = x)
  
  # Get final configuration
  mds_config <- lmds_config(
    dist_landmark = dist_to_lm,
    idx_landmark = idx_lm,
    x_landmark = x_lm,
    r = r,
    rescale = FALSE
  )
  rownames(mds_config$points) <- rownames(x)
  return(mds_config)
}

lmds_config <- function (dist_landmark, idx_landmark, x_landmark, r, rescale = TRUE){
  # Distance between landmark points and themselves
  dist_lm_lm <- dist_landmark[, idx_landmark, drop = FALSE]^2
  
  # Obtain number of landmark points
  num_landmarks <- as.integer(nrow(dist_lm_lm))
  
  # Obtain number of total points
  n_points <- as.integer(ncol(dist_landmark))
  
  # Double center the matrix
  mu_row_lm <- rowMeans(dist_lm_lm)
  mu_total_lm <- mean(dist_lm_lm)
  x_c <- sweep(dist_lm_lm, 1, mu_row_lm, "-")
  x_dc <- sweep(x_c, 2, mu_row_lm, "-") + mu_total_lm
  
  # Obtain SVD
  if (nrow(x_dc) > 10) {
    e <- svd::trlan.eigen(-x_dc/2, neig = r)
    ev <- e$d
    evec <- e$u
    if (ncol(evec) != r) {
      e <- eigen(-x_dc/2)
      ev <- e$values[1:r]
      evec <- e$vectors[, 1:r, drop = FALSE]
    }
  } else {
    e <- eigen(-x_dc/2)
    ev <- e$values[1:r]
    evec <- e$vectors[, 1:r, drop = FALSE]
  }
  
  # Get MDS for all the points
  points_inv <- evec/rep(sqrt(ev), each = num_landmarks)
  mds_config <- (-t(dist_landmark^2 - rep(mu_row_lm, n_points))/2) %*% points_inv
  mds_config <- apply(mds_config, MARGIN = 2, FUN = function(y) y - mean(y))
  mds_config <- mds_config %*% base::eigen(stats::cov(mds_config))$vectors
  mds <- list()
  mds$points <- mds_config
  mds$eigen <- ev/num_landmarks
  return(mds)
}

select_landmarks <- function (x, num_landmarks) {
  n_row_x <- nrow(x)
  if (num_landmarks > n_row_x) {
    stop("num_landmarks cannot be greater than nrow(x)")
  }
  
  idx_landmark <- sample.int(n_row_x, num_landmarks)
  return(idx_landmark)
}
