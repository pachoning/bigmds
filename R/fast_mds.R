get_partitions_for_fast <- function(n, l, s, k) {

  if (n/l<1.1) {
    p <- 2
  } else {
    p <- ceiling(l/s)
  }

  min_sample_size <- max(k + 2, s)
  size_partition <- floor(n/p)
  last_sample_size <- n - (p-1)*size_partition
  list_indexes <- list()

  if (size_partition < s) {
    stop("nrow(x) must be greater than s")
  } else if (size_partition < k) {
    stop("nrow(x)*s/l must be greater than k")
  }

  if (last_sample_size < min_sample_size & last_sample_size > 0) {
    p <- p - 1
  }

  for (i in 1:p) {
    if (i == 1) {
      ini <- 1
      end <- size_partition
    } else if (i < p) {
      ini <- end + 1
      end <- (ini - 1) + size_partition
    } else {
      ini <- end + 1
      end <- n
    }

    list_indexes[[i]] <- ini:end
  }

  return(list_indexes)
}


#'@title Fast MDS
#'@description Performs *Multidimensional Scaling* for big datasets using a recursive algorithm. This method can compute 
#'a MDS configuration even when the dataset is so large that classical MDS methods (`cmdscale`) can not be run 
#'due to computational problems.
#'@details In order to obtain a MDS configuration for the entire matrix \code{x}, it is partitioned into p submatrices,
#'where p=\code{l/s}.
#'
#'For every partition `cmdscale` is applied if the number of observations is less than \code{l}. Otherwise, `fast_mds`
#'is called. Notice that in this part is where this algorithm becomes a recursive one.
#' 
#'Once every submatrix has its own MDS configuration, \code{s} (random) points are taken from every partition of 
#'\code{x}. These points are put into a matrix *M*. Notice that *M* has \code{s}·p rows and q columns.
#'
#'After that, a MDS configuration for *M* is obtained. So, there are 2 configurations for the \code{s} points: one from 
#'performing MDS over every partition and another one from *M*. This allows to compute Procrustes (alignment method) so 
#'that all the MDS solutions share the same coordinate system.
#'
#'@param x A matrix with n individuals (rows) and q variables (columns).
#'@param l The largest value which allows classical MDS to be computed efficiently, i.e, the largest value which makes 
#'`cmdscale()` be run without any computational issues.
#'@param s Number of points used to align the MDS solutions obtained by the division of \code{x} into p submatrices.
#'Recommended value: \code{2·k}.
#'@param k Number of principal coordinates to be extracted.
#'@param dist_fn Distance function to be used for obtaining a MDS configuration.
#'@param ... Further arguments passed to \code{dist_fn} function.
#'@return Returns a list containing the following elements:
#'\describe{
#'   \item{points}{A matrix that consists of n individuals (rows) and \code{k} variables (columns) corresponding to the 
#'   MDS coordinates.}
#'   \item{eigen}{The first \code{k} eigenvalues.}
#'}
#'@examples
#'set.seed(42)
#'x <- matrix(data = rnorm(4*10000), nrow = 10000) %*% diag(c(15, 10, 1, 1))
#'mds <- fast_mds(x = x, l = 200, s = 2*2, k = 2, dist_fn = stats::dist)
#'head(cbind(mds$points, x[, 1:2]))
#'var(x)
#'var(mds$points)
#'@references
#'Tynia, Y., L. Jinze, M. Leonard, and W. Wei (2006). *A fast approximation to multidimensional scaling*. Proceedings of 
#'the ECCV Workshop on Computation Intensive Methods for Computer Vision (CIMCV).
#'\url{http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.79.2445}
#' 
#'Borg, I. and Groenen, P. (2005). *Modern Multidimensional Scaling: Theory and Applications*. Springer.
#'@export
fast_mds <- function(x, l, s, k, dist_fn = stats::dist, ...) {

  n <- nrow(x)

  if (n <= l) {

    mds <- classical_mds(x = x, k = k, dist_fn = dist_fn, ...)
    mds$eigen <- mds$eigen/nrow(x)

    return(mds)

  } else {

    # Split x
    index_partition <- get_partitions_for_fast(n = n, l = l, s = s, k = k)
    x_partition <- lapply(index_partition, function(idx, matrix) { matrix[idx, , drop = FALSE] }, matrix = x)
    num_partition <- length(index_partition)

    # Apply MDS to all the partitions
    mds_partition <- lapply(x_partition, fast_mds, l = l, s = s, k = k, dist_fn = dist_fn, ...)
    mds_partition_points <- lapply(mds_partition, function(x) x$points)
    mds_partition_eigen <- lapply(mds_partition, function(x) x$eigen)

    # take a subsample for each partition
    length_partition <- lapply(index_partition, length)
    sample_partition <- lapply(length_partition, sample, size = s, replace = FALSE)
    x_partition_sample <- mapply(function(matrix, idx) { matrix[idx, , drop = FALSE] }, 
                                 matrix = x_partition, idx = sample_partition, SIMPLIFY = FALSE)

    # Create two lists: one with initial position inside the matrix and another with end position
    ini_index <- list()
    end_index <- list()

    for (i in 1:num_partition) {
      length_sample_partition_i <- length(sample_partition[[i]])
      ini_index[[i]] <- (i-1)*length_sample_partition_i + 1
      end_index[[i]] <- i*length_sample_partition_i
    }

    # Join each sampled data
    x_M <- Reduce(rbind, x_partition_sample)

    # Apply MDS to the subsampling points
    mds_M <- classical_mds(x = x_M, k = k, dist_fn = stats::dist, ...)
    mds_M_points <- mds_M$points

    # Extract the MDS configuration for the sampling points from mds_M_points 
    mds_M_sampling_points <- mapply(function(matrix, ini, end) {  matrix[ini:end, , drop = FALSE] }, 
                                    ini = ini_index, end = end_index, 
                                    MoreArgs = list(matrix = mds_M_points), SIMPLIFY = FALSE)

    # Extract the MDS configuration for the sampling points from mds_partition_points
    mds_partition_sampling_points <- mapply(function(matrix, idx) { matrix[idx, , drop = FALSE] }, 
                                            matrix = mds_partition_points, idx = sample_partition, SIMPLIFY = FALSE)

    # Apply Procrustes
    procrustes <- mapply(perform_procrustes, x = mds_partition_sampling_points, target = mds_M_sampling_points, 
                         matrix_to_transform = mds_partition_points, translation = FALSE, dilation = FALSE, SIMPLIFY = FALSE)

    # Build the list to be returned
    mds <- Reduce(rbind, procrustes)
    eigen <- Reduce(`+`, mds_partition_eigen)/num_partition

    return(list(points = mds, eigen = eigen))
  }
}
