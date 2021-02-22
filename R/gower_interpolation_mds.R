get_partitions_for_gower_interpolation <- function(n, l, k) {

  if (l<=k) {
    stop("l must be greater than k")
  }

  p <- ceiling(n/l)
  p <- pmax(1, p)

  list_index <- list()

  idexes <- sample(x = n, size = n)

  for (i in 1:p) {

    ini <- (i-1)*l + 1

    if (i == p) {
      end <- n
    } else{
      end <- i*l
    }

    list_index[[i]] <- idexes[ini:end]
  }

  return(list_index)
}


#'@title MDS based on Gower interpolation formula
#'@description Performs *Multidimensional Scaling* for big datasets using Gower interpolation formula. This method can 
#'compute a MDS configuration even when the dataset is so large that classical MDS methods (`cmdscale`) can not be run 
#'due to computational problems.
#'@details *Gower interpolation formula* is the central piece of this algorithm since it allows to add a new set of 
#'points to an existing MDS configuration so that the new one has the same coordinate system. 
#'
#'Given the matrix \code{x} with n individuals (rows) and q variables (columns), a submatrix based on a random sample 
#'of \code{l} individuals is taken and it is used to compute a MDS configuration. 
#'
#'The remaining part of \code{x} is divided into p=(n-\code{l})/\code{l} submatrices. For every submatrix, it is obtained 
#'a MDS configuration by means of *Gower interpolation formula* and the first (random) submatrix. Every MDS 
#'configuration is appended to the existing one so that, at the end of the process, a MDS configuration for \code{x} is 
#'built.
#'@param x A matrix with n individuals (rows) and q variables (columns).
#'@param l The largest value which allows classical MDS to be computed efficiently, i.e, the largest value which makes 
#'`cmdscale()` be run without any computational issues.
#'@param k Number of principal coordinates to be extracted.
#'@param dist_fn Distance function to be used for obtaining a MDS configuration.
#'@param ... Further arguments passed to \code{dist_fn} function.
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
#'Borg, I. and Groenen, P. (2005). *Modern Multidimensional Scaling: Theory and Applications*. Springer.
#' @export
gower_interpolation_mds <- function(x, l, k, dist_fn = stats::dist, ...) {

  n <- nrow(x)
  idexes_partition <- get_partitions_for_gower_interpolation(n = n, l = l, k = k)
  num_partitions <- length(idexes_partition)

  if (num_partitions <= 1) {

    # It is possible to run MDS directly
    mds <- classical_mds(x = x, k = k, dist_fn = dist_fn, ...)
    points <- mds$points
    eigen <- mds$eigen/n
    list_to_return <- list(points = points, eigen = eigen)

  } else {

    # Get the first group 
    ind_1 <-idexes_partition[[1]]
    n_1 <- length(ind_1)
    idexes_partition[[1]] <- NULL

    # Obtain MDS for the first group
    x_1 <- x[ind_1, , drop = FALSE]
    mds_eig <- classical_mds(x = x_1, k = k, dist_fn = dist_fn, return_distance_matrix = TRUE, ...)
    distance_matrix <- mds_eig$distance

    M <- mds_eig$points
    eigen <- mds_eig$eigen/nrow(M)

    # Calculations needed to do Gower interpolation
    delta_matrix <- distance_matrix^2
    In <- diag(n_1)
    ones_vector <- rep(1, n_1)
    J <- In - 1 / n_1 * ones_vector %*% t(ones_vector)
    G <- -1 / 2 * J %*% delta_matrix %*% t(J) 
    g_vector <- diag(G)
    S <- 1 / (nrow(M)-1) * t(M) %*% M
    S_inv <- solve(S)

    # Get x for each partition
    x_other <- lapply(idexes_partition, function(matrix, idx) { matrix[idx, , drop = FALSE] }, matrix = x)

    # Obtain the distance matrix with respect the first partition
    distance_matrix_filter <- lapply(x_other, function(X, Y){ pdist::pdist(X, Y) }, Y = x_1)
    distance_matrix_filter <- lapply(distance_matrix_filter, as.matrix)

    # A matrix
    A <- lapply(distance_matrix_filter, function(x){ x^2 })
    ones_vector <- lapply(idexes_partition, function(times, x){ rep(x, length(times)) }, x = 1)

    # Get MDS for all the partitions
    MDS <- mapply(function(A, ones_vector) { 1 / (2 * n_1) * (ones_vector %*% t(g_vector) - A) %*% M %*% S_inv  }, 
                  A = A, ones_vector = ones_vector, SIMPLIFY = FALSE)

    # Get cummulative MDS
    cum_mds <- Reduce(rbind, MDS)
    cum_mds <- Reduce(rbind, list(M, cum_mds))

    # Reorder the rows
    idexes_order <- Reduce(c, idexes_partition)
    idexes_order <- Reduce(c, list(ind_1, idexes_order))
    cum_mds <- cum_mds[order(idexes_order), , drop = FALSE]

    # List to return
    list_to_return <- list(points = cum_mds, eigen = eigen)
  }

  return(list_to_return)

}
