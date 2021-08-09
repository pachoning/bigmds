get_partitions_for_divide_conquer <- function(n, l, c_points, r) {

  if (l-c_points <= 0) {
    stop("\"l\" must be greater than \"c_points\"")
  } else if (l-c_points <= c_points) {
    stop("\"l-c_points\" must be greater than \"c_points\"")
  } else if (l-c_points <= r) {
    stop ("\"l-c_points\" must be greater than \"r\"")
  }

  permutation <- sample(x = n, size = n, replace = FALSE)

  if (n<=l) {
    list_indexes <- list(permutation)
  } else {
    min_sample_size <- max(r+2, c_points)
    p <- 1 + ceiling((n-l)/(l-c_points))
    last_partition_sample_size <- n - (l + (p-2) * (l-c_points))

    if (last_partition_sample_size < min_sample_size & last_partition_sample_size > 0) {
      p <- p - 1
      last_partition_sample_size <- n - (l + (p-2) * (l-c_points))
    }

    first_parition <- permutation[1:l]
    last_partition <- permutation[(n-last_partition_sample_size+1):n]
    list_indexes <- split(x = permutation[(l+1):(n-last_partition_sample_size)], f = 1:(p-2))
    names(list_indexes) <- NULL
    list_indexes[[p-1]] <- list_indexes[[1]]
    list_indexes[[p]] <- last_partition
    list_indexes[[1]] <- first_parition
  }

  return(list_indexes)
}

split_matrix <- function(matrix, num_points) {
  x_1 <- matrix[1:num_points, , drop = FALSE]
  x_2 <- matrix[(num_points + 1):nrow(matrix), , drop = FALSE]
  return(list(x_1, x_2))
}

main_divide_conquer_mds <- function(idx, x, x_sample_1, r, original_mds_sample_1, dist_fn, ...) {
  # Filter the matrix
  x_filtered <- x[idx, , drop = FALSE]

  # Join with the sample from the first partition
  x_join_sample_1 <- rbind(x_sample_1, x_filtered)

  # Perform MDS
  mds_all <- classical_mds(x = x_join_sample_1, k = r, dist_fn = dist_fn, ...)
  mds_points <- mds_all$points
  mds_eigen <- mds_all$eigen
  mds_GOF <- mds_all$GOF

  # Split MDS into two parts: the first part is the MDS for the sampling points of the
  # first partition and the second part is the MDS for the points of the partition
  mds_split <- split_matrix(matrix = mds_points, num_points = nrow(x_sample_1))
  mds_sample_1 <- mds_split[[1]]
  mds_partition <- mds_split[[2]]
  mds_procrustes <- perform_procrustes(x = mds_sample_1, 
                                       target = original_mds_sample_1, 
                                       matrix_to_transform = mds_partition, 
                                       translation = FALSE)

  mds_eigen <- mds_eigen/length(idx)
  mds_GOF <- mds_GOF*length(idx)

  return(list(points = mds_procrustes, eigen = mds_eigen, GOF = mds_GOF))
}

#'@title Divide-and-conquer MDS
#'
#'@description Roughly speaking, a large data set, \code{x}, of size \eqn{n}  
#'is divided into parts, then classical MDS is performed over every part and, 
#'finally, the partial configurations are combined so that all the points lie 
#'on the same coordinate system with the aim to obtain a global MDS configuration.
#'
#'@details The divide-and-conquer MDS starts dividing the \eqn{n} points into 
#'\eqn{p} partitions: the first partition contains \code{l} points and the others
#'contain \code{l-c_points} points. Therefore, \eqn{p = 1 + (n-}\code{l)/(l-c_points)}.
#'The partitions are created at random. 
#'
#'Once the partitions are created, \code{c_points} different random 
#'points are taken from the first partition and concatenated to the other 
#'partitions After that, classical MDS is applied to each partition, 
#'with target low dimensional configuration \code{r}.
#'
#'Since all the partitions share \code{c_points}
#'points with the first one, Procrustes can be applied in order to align all
#'the configurations. Finally, all the configurations are
#'concatenated in order to obtain a global MDS configuration.
#'
#'@param x A matrix with \eqn{n} points (rows) and \eqn{k} variables (columns).
#'@param l The size for which classical MDS can be computed efficiently 
#'(using `cmdscale` function). It means that if \eqn{\bar{l}} is the limit 
#'size for which classical MDS is applicable, then \code{l}\eqn{\leq \bar{l}}.
#'@param c_points Number of points used to align the MDS solutions obtained by the 
#'division of \code{x} into \eqn{p} data subsets. Recommended value: \code{2·r}.
#'@param r Number of principal coordinates to be extracted.
#'@param n_cores Number of cores wanted to use to run the algorithm.
#'@param dist_fn Distance function used to compute the distance between the rows.
#'@param ... Further arguments passed to \code{dist_fn} function.
#'
#'@return Returns a list containing the following elements:
#' \describe{
#'   \item{points}{A matrix that consists of \eqn{n} points (rows) 
#'   and \code{r} variables (columns) corresponding to the MDS coordinates. Since 
#'   a dimensionality reduction is performed, \code{r}\eqn{<<k}}
#'   \item{eigen}{The first \code{r} largest eigenvalues: 
#'   \eqn{\bar{\lambda}_i, i \in  \{1, \dots, r\} }, where
#'   \eqn{\bar{\lambda}_i = 1/p \sum_{j=1}^{p}\lambda_j/n_j}, 
#'   being \eqn{n_j} the size of the partition \eqn{j}.}
#'   \item{GOF}{A numeric vector of length 2. 
#'   
#'   The first element corresponds to
#'   \eqn{1/n \sum_{j=1}^{p}n_jG_1^j}, where 
#'   \eqn{G_1^j = \sum_{i = 1}^{r} \lambda_{i}^{j}/ \sum_{i = 1}^{n-1} |\lambda_{i}^{j}|}. 
#'   
#'   The second element corresponds to 
#'    \eqn{1/n \sum_{j=1}^{p}n_jG_2^j} where 
#'    \eqn{G_2^j = \sum_{i = 1}^{r} \lambda_{i}^{j}/ \sum_{i = 1}^{n-1} max(\lambda_{i}^{j}, 0).}}
#'}
#'
#'@examples
#'set.seed(42)
#'x <- matrix(data = rnorm(4*10000), nrow = 10000) %*% diag(c(9, 4, 1, 1))
#'mds <- divide_conquer_mds(x = x, l = 200, c_points = 2*2, r = 2, n_cores = 1, dist_fn = stats::dist)
#'head(mds$points)
#'mds$eigen
#'mds$GOF
#'points <- mds$points
#'plot(x[1:10, 1],
#'      x[1:10, 2],
#'      xlim = range(c(x[1:10,1],points[1:10,1])),
#'      ylim = range(c(x[1:10,2], points[1:10,2])),
#'      pch = 19,
#'      col = "green")
#'text(x[1:10, 1], x[1:10, 2], labels=1:10)
#'points(points[1:10, 1], points[1:10, 2], pch = 19, col = "orange")
#'text(points[1:10, 1], points[1:10, 2], labels=1:10)
#'abline(v = 0, lwd=3, lty=2)
#'abline(h = 0, lwd=3, lty=2)
#'
#'@references
#'Delicado P. and C. Pachon-Garcia (2021). *Multidimensional Scaling for Big Data*.
#'\url{https://arxiv.org/abs/2007.11919}.
#'
#'Borg, I. and P. Groenen (2005). *Modern Multidimensional Scaling: Theory and Applications*. Springer.
#'
#'@export
divide_conquer_mds <- function(x, l, c_points, r, n_cores = 1, dist_fn = stats::dist, ...) {

  n_row_x <- nrow(x)

  if (n_row_x <= l) {
    mds_to_return <- classical_mds(x = x, k = r, dist_fn = dist_fn, ...)
    mds_to_return$eigen <- mds_to_return$eigen/n_row_x
    mds_to_return$GOF <- mds_to_return$GOF

  } else {

    mds_matrix <- matrix(data = NA, nrow = n_row_x, ncol = r)

    # Generate indexes list. Each element corresponds to the index of the partition
    idx <- get_partitions_for_divide_conquer(n = n_row_x, l = l, c_points = c_points, r = r)
    num_partitions <- length(idx)
    length_1 <- length(idx[[1]])

    # Perform MDS for the first partition
    x_1 <- x[idx[[1]], , drop = FALSE]
    mds_1 <- classical_mds(x = x_1, k = r, dist_fn = dist_fn, ...)
    mds_1_points <- mds_1$points
    mds_1_eigen <- mds_1$eigen/length_1
    mds_1_GOF <- mds_1$GOF*length_1

    # Take a sample from the first partition
    sample_1 <- sample(x = length_1, size = c_points, replace = FALSE)
    x_sample_1 <- x_1[sample_1, ,drop = FALSE]
    mds_sample_1 <- mds_1_points[sample_1, , drop = FALSE]

    mds_others_results <- parallel::mclapply(idx[2:num_partitions],
                                             main_divide_conquer_mds,
                                             x = x,
                                             x_sample_1 = x_sample_1,
                                             r = r,
                                             original_mds_sample_1 = mds_sample_1,
                                             dist_fn = dist_fn,
                                             mc.cores = n_cores,
                                             ...)

    # Obtain points
    mds_others_points <- do.call(rbind, parallel::mclapply(mds_others_results, function(x) x$points, mc.cores = n_cores))
    mds_matrix[1:length_1, ] <- mds_1_points
    mds_matrix[(length_1 + 1):n_row_x, ] <- mds_others_points
    order_idx <- do.call(c, idx)
    order_idx <- order(order_idx)
    mds_matrix <- mds_matrix[order_idx, , drop = FALSE]
    mds_matrix <- apply(mds_matrix, MARGIN = 2, FUN = function(x) x - mean(x))
    mds_matrix <- mds_matrix %*% eigen(stats::cov(mds_matrix))$vectors

    # Obtain eigenvalues
    eigen <- parallel::mclapply(mds_others_results, function(x) x$eigen, mc.cores = n_cores)
    eigen[[num_partitions]] <- mds_1_eigen
    eigen <- Reduce(`+`, eigen)
    eigen <- eigen/num_partitions

    # Obtain GOF
    GOF <- parallel::mclapply(mds_others_results, function(x) x$GOF, mc.cores = n_cores)
    GOF[[num_partitions]] <- mds_1_GOF
    GOF <- Reduce(`+`, GOF)
    GOF <- GOF/n_row_x
    mds_to_return <- list(points = mds_matrix, eigen = eigen, GOF = GOF)
  }

  return(mds_to_return)
}
