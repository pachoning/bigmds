#'@title Pivot MDS
#'
#'@description Pivot MDS, introduced in the literature of graph layout algorithms, is similar to 
#'Landmark MDS ([landmark_mds()]) but it uses the distance information between landmark and non-landmark 
#'points to improve the initial low dimensional configuration, 
#'as more relations than just those between landmark points are taken into account.
#'
#'@param x A matrix with \eqn{n} individuals (rows) and \eqn{k} variables (columns).
#'@param num_pivots Number of pivot points to obtain an initial MDS configuration. It is
#'equivalent to \code{l} parameter used in [interpolation_mds()], [divide_conquer_mds()] and
#'[fast_mds()]. Therefore, it is the size for which classical MDS can be computed efficiently 
#'(using `cmdscale` function). It means that if \eqn{\bar{l}} is the limit 
#'size for which classical MDS is applicable, then \code{l}\eqn{\leq \bar{l}}.
#'@param r Number of principal coordinates to be extracted.
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
#'Brandes U. and C. Pich (2007). *Eigensolver Methods for Progressive Multidimensional Scaling of Large Data*. Graph Drawing.
#'
#'Borg, I. and P. Groenen (2005). *Modern Multidimensional Scaling: Theory and Applications*. Springer.
#'
#'Gower JC. (1968). *Adding a point to vector diagrams in multivariate analysis*. Biometrika.
#'
#'@examples
#'set.seed(42)
#'x <- matrix(data = rnorm(4 * 10000), nrow = 10000) %*% diag(c(9, 4, 1, 1))
#'mds <- pivot_mds(x = x, num_pivots = 200, r = 2)
#'head(mds$points)
#'mds$eigen
#'
#'@importFrom pracma distmat
#'@importFrom svd propack.svd
#'
#'@export
pivot_mds <- function(x, num_pivots, r) {
  n_row_x <- nrow(x)
  if (num_pivots > n_row_x) {
    stop("num_pivots cannot be greater than nrow(x)")
  }
  
  # Get indexes for pivot points
  idx_pivots <- sample.int(n_row_x, num_pivots)
  x_pivots <- x[idx_pivots, , drop = FALSE]
  
  # Compute distance matrix
  D <- pracma::distmat(X = x_pivots, Y = x)
  D <- t(D)
  D_sq <- D^2
  
  # Apply formula
  col_mean <- colMeans(D_sq)
  row_mean <- rowMeans(D_sq)
  Dmat <- -(1/2) * (D_sq - outer(row_mean, col_mean, function(x, y) x + y) + mean(D_sq))
  
  # Use svd
  if (n_row_x > 10) {
    svd_results <- svd::propack.svd(Dmat, neig = r)
    u_matrix <- svd_results$u
    singular_values <- svd_results$d
    
    if (ncol(u_matrix) != r) {
      svd_results <- svd(Dmat, nu = r, nv = r)
      u_matrix <- svd_results$u
      singular_values <- svd_results$d
    }
  } else {
    svd_results <- svd(Dmat, nu = r, nv = r)
    u_matrix <- svd_results$u
    singular_values <- svd_results$d
  }
  
  # Compute MDS
  eigen <- singular_values/sqrt((n_row_x * num_pivots))
  points <- u_matrix %*% diag(sqrt(n_row_x * eigen))
  mds <- list(
    points = points,
    eigen = eigen
  )
  return(mds)
}
