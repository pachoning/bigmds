get_partitions_for_divide_conquer <- function(n, l, tie, k) {

  size_partition <- l-tie
  p <- ceiling(n/size_partition)
  min_sample_size <- max(k, tie)
  partition_sample_size <- l-tie
  last_sample_size <- n-(p-1)*partition_sample_size
  list_indexes <- list()

  if (l-tie <= 0) {
    stop("l must be greater than tie")
  } else if(l-tie <= tie) {
    stop("l-tie must be greater than tie")
  } else if(l-tie <= k) {
    stop("l-tie must be greater than k")
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
      end <- (ini-1) + size_partition
    } else {
      ini <- end + 1
      end <- n
    }

    list_indexes[[i]] <- ini:end
  }

  return(list_indexes)

}


divide_matrix <- function(x, long) {

  n_row <- nrow(x)
  x_first <- x[1:long, , drop = FALSE]
  x_rest <- x[(long+1):n_row, , drop = FALSE]
  return(list(first = x_first, rest = x_rest))

}


#'@title Divide and Conquer MDS
#'@description Performs *Multidimensional Scaling* for big datasets using a Divide and Conquer strategy. This method can 
#'compute a MDS configuration even when the dataset is so large that classical MDS methods (`cmdscale`) can not be run 
#'due to computational problems.
#'@details In order to obtain a MDS configuration for the entire matrix \code{x}, it is needed to break the dataset into 
#'p submatrices (*Divide and Conquer strategy*).
#' 
#'In order to obtain p, \code{tie} and \code{l} are taken into account: p=n/\code{(l-tie)}. This allows to use 
#'`cmdscale` function in every submatrix.
#'
#'Taking into account that given a MDS solution, any rotation is another (valid) MDS solution, it is needed a way to 
#'obtain the same coordinate system for all the partitions. 
#'
#'To achieve such a common coordinate system, the algorithm starts by taking the first partition and calculating a MDS 
#'configuration as well as a subsample of size \code{tie} (from the partition, not from its MDS configuration). 
#'These \code{tie} points will be used in order to force the other partitions to have the same coordinate system as the 
#'first one.
#'
#'Given a partition, the \code{tie} points are appended to it. After that, a MDS configuration is obtained. Therefore, 
#'for these \code{tie} points there are two MDS solutions. In order to aligned them, Procrustes parameters are 
#'obtained. These parameters are applied to the MDS configuration of the partition.
#'@param x A matrix with n individuals (rows) and q variables (columns).
#'@param l The largest value which allows classical MDS to be computed efficiently, i.e, the largest value which makes 
#'`cmdscale()` be run without any computational issues.
#'@param tie Number of points used to align the MDS solutions obtained by the division of \code{x} into p submatrices.
#'Recommended value: \code{2·k}.
#'@param k Number of principal coordinates to be extracted.
#'@param dist_fn Distance function to be used for obtaining a MDS configuration.
#'@param ... Further arguments passed to \code{dist_fn} function.
#'@return Returns a list containing the following elements:
#' \describe{
#'   \item{points}{A matrix that consists of n individuals (rows) and \code{k} variables (columns) corresponding to the 
#'   MDS coordinates.}
#'   \item{eigen}{The first \code{k} eigenvalues.}
#'}
#'@examples
#'set.seed(42)
#'x <- matrix(data = rnorm(4*10000), nrow = 10000) %*% diag(c(15, 10, 1, 1))
#'mds <- divide_conquer_mds(x = x, l = 200, tie = 2*2, k = 2, dist_fn = stats::dist)
#'head(cbind(mds$points, x[, 1:2]))
#'var(x)
#'var(mds$points)
#'@references
#'Delicado P. and C. Pachon-Garcia (2020). *Multidimensional Scaling for Big Data*.
#'\url{https://arxiv.org/abs/2007.11919}
#' 
#'Borg, I. and Groenen, P. (2005). *Modern Multidimensional Scaling: Theory and Applications*. Springer.
#'@export
divide_conquer_mds <- function(x, l, tie, k, dist_fn = stats::dist, ...) {

  n_row_x <- nrow(x)

  if (n_row_x <= l) {
    mds_to_return <- classical_mds(x = x, k = k, dist_fn = dist_fn, ...)
    mds_to_return$eigen <- mds_to_return$eigen/nrow(x)

  } else {
    # Generate indexes list. Each element correspond to the index of the partition
    idx <- get_partitions_for_divide_conquer(n = n_row_x, l = l, tie = tie, k = k)
    num_partitions <- length(idx)

    # Get the elements of the first partition and drop from the list
    idx_1 <- idx[[1]]
    idx[[1]] <- NULL 

    # Partition x
    x_partition <- lapply(idx, function(rows, matrix) matrix[rows, , drop = FALSE], matrix = x)

    # Take a sample from the first partition
    sample_first_partition <- sample(x = idx_1, size = tie, replace = FALSE)
    x_1 <- x[idx_1, , drop = FALSE]
    x_sample_1 <- x_1[sample_first_partition, , drop = FALSE]

    # Join each partition with the sample from the first partition
    x_join_1 <- lapply(x_partition, function(m_big, m_small) rbind(m_small, m_big), m_small = x_sample_1)

    # Perform MDS for each partition as well as for the first partition
    mds_1 <- classical_mds(x = x_1, k = k, dist_fn = dist_fn, return_distance_matrix = FALSE, ...)
    mds_1_points <- mds_1$points
    mds_1_eigen <- mds_1$eigen/nrow(x_1)
    mds_1_sample <- mds_1_points[sample_first_partition, ,drop = FALSE]

    mds_join_1 <- lapply(x_join_1, classical_mds, k = k, dist_fn = dist_fn, return_distance_matrix = FALSE)
    mds_join_1_points <- lapply(mds_join_1, function(x) x$points)
    mds_join_1_eigen <- lapply(mds_join_1, function(x) x$eigen)

    # For each partition, divide the matrix into two part: 
    # first corresponding to the sample of the first partition
    # the rest of the matrix corresponding to the observations of each matrix
    mds_division <- lapply(mds_join_1_points, divide_matrix, long = tie)
    mds_division_first <- lapply(mds_division, function(x) x$first)
    mds_division_rest <- lapply(mds_division, function(x) x$rest)

    # Apply Procrustes for each partition
    mds_procrustes <- mapply(FUN = perform_procrustes, 
                             x = mds_division_first, 
                             matrix_to_transform = mds_division_rest,
                             MoreArgs = list(target = mds_1_sample, translation = FALSE, dilation = FALSE))

    # Join all the solutions
    mds_solution <- Reduce(rbind, mds_procrustes)
    mds_solution <- Reduce(rbind, list(mds_1_points, mds_solution))

    # Get eigenvalues
    eigen <- mapply(function(x, y) x/length(y), x = mds_join_1_eigen, y = idx, SIMPLIFY = FALSE)
    eigen <- Reduce(`+`, eigen)
    eigen <- Reduce(`+`, list(mds_1_eigen, eigen))
    eigen <- eigen/num_partitions

    mds_to_return <- list(points = mds_solution, eigen = eigen)
  }

  return(mds_to_return)
}
