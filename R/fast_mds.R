get_partitions_for_fast <- function(n, l, s, k) {

  p <- ceiling(l / s)
  min_sample_size <- max(k + 2, s)

  if (ceiling(n / p) < min_sample_size) {
    stop("Too many columns and few observations to perform Fast MDS")
  }

  partition <- sort(rep(x = 1:p, length.out = n, each = ceiling(n / p)))
  p <- max(partition)

  while (p <= n & min(table(partition)) < min_sample_size) {
    p <- p + 1
    partition <- sort(rep(x = 1:p, length.out = n, each = ceiling(n / p)))
  }

  if (min(table(partition)) < min_sample_size) {
    stop("Partitions for Fast MDS suffer from lacking of data")
  }

  return(partition)
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
#'Borg I and P. Groenen (1997). *Modern Multidimensional Scaling: Theory and Applications*. New York: Springer. pp. 340-342.
#'@export
fast_mds <- function(x, l, s, k, dist_fn = stats::dist, ...) {

  has_row_names <- !is.null(row.names(x))
  if (!has_row_names) {
    row.names(x) <- 1:nrow(x)
  }

  #If possible to run classical MDS on the whole matrix, run it
  if (nrow(x) <= l) {
    mds <- classical_mds(x = x, k = k, dist_fn = dist_fn, ...)
    mds$eigen <- mds$eigen / length(mds$eigen)

    if (!has_row_names) {
      row.names(x) <- NULL
      row.names(mds) <- NULL
    }

    return(mds)

    # Otherwise, call it recursively
  } else {
    points <- list()
    eigen <- c()
    min_len <- NA
    sampling_points <- list()

    index_partition <- get_partitions_for_fast(n = nrow(x), l = l, s = s, k = k)
    p <- length(unique(index_partition))

    # For each partition, compute fast MDS
    for (i in 1:p) {
      indexes_partition <- which(index_partition == i)
      x_partition <- x[indexes_partition, ,drop = FALSE]
      mds_partition <- fast_mds(x = x_partition, l = l, s = s, k = k, dist_fn = dist_fn, ...)
      points[[i]] <- mds_partition$points
      row.names(points[[i]]) <- row.names(x_partition)
      sampling_points[[i]] <- sample(x = row.names(x_partition), size = s, replace = FALSE)

      if (i == 1) {
        min_len <- length(mds_partition$eigen)
        eigen <- mds_partition$eigen
      } else {
        min_len <- pmin(min_len, length(mds_partition$eigen))
        eigen <- eigen[1:min_len] + mds_partition$eigen[1:min_len]
      }
    }

    # Perform the mean for the eigenvalues
    eigen <- eigen / p

    # Get M_align by getting the sampled points
    ind <- unlist(sampling_points)
    x_M <- x[ind, ,drop = FALSE]
    row.names(x_M) <- row.names(x[ind, ,drop = FALSE])
    mds_M <- classical_mds(x = x_M, k = k, dist_fn = dist_fn, ...)
    mds_M <- mds_M$points
    row.names(mds_M) <- row.names(x_M)

    #Align the solutions
    for (i in 1:p) {
      mds_i <- points[[i]]
      sampling_points_i <- sampling_points[[i]] 
      mds_M_sampling <- mds_M[sampling_points_i, ,drop = FALSE]
      mds_i_sampling <- mds_i[sampling_points_i, ,drop = FALSE]

      mds_aligned_i <- perform_procrustes(x = mds_i_sampling, target = mds_M_sampling, matrix_to_transform = mds_i, 
                                         translation = FALSE, dilation = FALSE)
      if (i == 1) {
        mds_stitched <- mds_aligned_i
      } else {
        mds_stitched <- rbind(mds_stitched, mds_aligned_i)
      }
    }

    if (!has_row_names) {
      row.names(x) <- NULL
      row.names(mds_stitched) <- NULL
    }

    return(list(points = mds_stitched, eigen = eigen))
  }
}
