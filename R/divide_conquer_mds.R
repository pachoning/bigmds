get_partitions_for_divide_conquer <- function(n, l, tie, k) {

  p <- ceiling(n / (l-tie))
  min_sample_size <- max(k + 2, tie)

  index_partition <- sort(rep(x = 1:p, length.out = n, each = ceiling(n / p)))
  p <- max(index_partition)

  if (mean(table(index_partition)) < min_sample_size) {
    stop("Too many columns and few observations to perform Divide and Conquer MDS")
  }

  while ((min(table(index_partition)) < min_sample_size) & (p <= n)) {
    p <- p + 1
    index_partition <- sort(rep(x = 1:p, length.out = n, each = ceiling(n / p)))
  }

  if (min(table(index_partition)) < min_sample_size) {
    stop("Partitions for divide and conquer suffer from lacking of data")
  }

  return(index_partition)
}

#'@title Divide and Conquer MDS
#'@description Performs *Multidimensional Scaling* for big datasets using a Divide and Conquer strategy. This method can 
#'compute a MDS configuration even when the dataset is so large that classical MDS methods (`cmdscale`) can not be run 
#'due to computational problems.
#'@details In order to obtain a MDS configuration for the entire matrix \code{k}, it is needed to break the dataset into 
#'p submatrices (*Divide andCconquer strategy*).
#' 
#'In order to obtain p, \code{tie} as well as \code{l} parameters are taken into account: p=n/\code{(l-tie)}. This
#'allows to use `cmdscale` function in every submatrix.
#'
#'Given a MDS solution, any rotation is another (valid) MDS solution. It means that every partition has its 
#'own coordinate system. 
#'
#'In order to keep the same coordinate system, every two consecutive submatrices are forced to share \code{tie} 
#'points. These points are used to align the corresponding MDS configurations. Such an alignment is performed by means 
#'of Procrustes.
#'@param x A matrix with n individuals (rows) and q variables (columns).
#'@param l The largest value which allows classical MDS to be computed efficiently, i.e, the largest value which makes 
#'`cmdscale()` be run without any computational issues.
#'@param tie Number of points used to align the MDS solutions obtained by the division of \code{x} into p submatrices.
#'Recommended value: \code{2·k}.
#'@param k Number of principal coordinates to be extracted.
#'@param dist_fn Distance function to be used for obtaining a MDS configuration.
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
#'Borg I and P. Groenen (1997). *Modern Multidimensional Scaling: Theory and Applications*. New York: Springer. pp. 340-342.
#'@export
divide_conquer_mds <- function(x, l, tie, k, dist_fn = stats::dist) {

  initial_row_names <- row.names(x)
  row.names(x) <- 1:nrow(x)

  if (nrow(x) <= l) {
    mds_to_return <- classical_mds(x = x, k = k, dist_fn = dist_fn)
    mds_to_return$eigen <- mds_to_return$eigen / length(mds_to_return$eigen)
  } else {
    index_partition <- get_partitions_for_divide_conquer(n = nrow(x), l = l, tie = tie, k = k)
    p <- max(index_partition)

    min_len = NA
    eigen <- c()

    # Calculate mds for each partition and take tie poits from each subsample
    for (i in 1:p) {

      indexes_current <- which(index_partition == i)
      x_current <- x[indexes_current, ,drop = FALSE]
      row_names_current <- row.names(x_current)

      if (i == 1) {

        list_classical_mds <- classical_mds(x = x_current, k = k, dist_fn = dist_fn) 
        cum_mds <- list_classical_mds$points
        eigen <- list_classical_mds$eigen / length(list_classical_mds$eigen)
        min_len <- length(eigen)
      } else {

        list_mds_both <- classical_mds(x = x[c(rn_subsample_previous, row_names_current), ,drop = FALSE], k = k)
        mds_both <- list_mds_both$points
        mds_both_previous <- mds_both[rn_subsample_previous, ,drop = FALSE]
        mds_both_current <- mds_both[row_names_current, ,drop = FALSE]
        cum_mds_previous <- cum_mds[rn_subsample_previous, ,drop = FALSE]
        mds_current_aligned <- perform_procrustes(x = mds_both_previous, target = cum_mds_previous,
                                                  matrix_to_transform = mds_both_current, 
                                                 translation = FALSE, dilation = FALSE)
        cum_mds <- rbind(cum_mds, mds_current_aligned)
        min_len <- pmin(min_len, length(list_mds_both$eigen))
        eigen <- eigen[1:min_len] + (list_mds_both$eigen[1:min_len] / length(list_mds_both$eigen))
      }

      rn_subsample_previous <- sample(x = row_names_current, size = tie, replace = FALSE)

    }

    # Perform the mean for the eigenvalues
    eigen <- eigen / p

    # Divide by the number of observations
    mds_to_return <- list(points = cum_mds, eigen = eigen)
  }

  row.names(x) <- initial_row_names
  row.names(mds_to_return$points) <- initial_row_names

  return(mds_to_return)
}
