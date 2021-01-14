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
#'@description Performs Multidimensional Scaling based on Yang, Tynia & Liu, Jinze & Mcmillan, Leonard & Wang, Wei. (2006).
#'@param x Data matrix.
#'@param l The highest value where classical MDS can be computed efficiently.
#'@param s Number of sampling points. It should be 2 x estimated data dimension.
#'@param k Number of principal coordinates.
#'@return Returns MDS based on Fast MDS algorithm as well as the first k eigenvalues.
#' \describe{
#'   \item{points}{MDS}
#'   \item{eigen}{eigenvalues}
#' }
#' @export
#' @examples
#' x <- matrix(data = rnorm(4*10000, sd = 10), nrow = 10000)
#' cmds <- fast_mds(x = x, l = 100, s = 8, k = 2)
#' head(cmds$points)
#' cmds$eigen
#' @seealso
#' \url{https://arxiv.org/abs/2007.11919}
fast_mds <- function(x, l, s, k) {

  has_row_names <- !is.null(row.names(x))
  if (!has_row_names) {
    row.names(x) <- 1:nrow(x)
  }

  #If possible to run classical MDS on the whole matrix, run it
  if (nrow(x) <= l) {
    mds <- classical_mds(x = x, k = k)
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
      mds_partition <- fast_mds(x = x_partition, l = l, s = s, k = k)
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
    mds_M <- classical_mds(x = x_M, k = k)
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
