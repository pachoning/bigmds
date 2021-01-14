get_partitions_for_divide_conquer <- function(n, l, s, k) {

  p <- ceiling(n / (l-s))
  min_sample_size <- max(k + 2, s)

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
#'@description Performs Multidimensional Scaling based on Delicado and Pachon-Garcia, 2020. 
#'@param x Data matrix.
#'@param l The highest value where classical MDS can be computed efficiently.
#'@param s Number of sampling points. It should be 2 x estimated data dimension.
#'@param k Number of principal coordinates.
#'@return Returns MDS based on Divide and Conquer MDS as well as the first k eigenvalues.
#' \describe{
#'   \item{points}{MDS}
#'   \item{eigen}{eigenvalues}
#' }
#' @export
#' @examples
#' x <- matrix(data = rnorm(4*10000, sd = 10), nrow = 10000)
#' cmds <- divide_conquer_mds(x = x, l = 100, s = 8, k = 2)
#' head(cmds$points)
#' cmds$eigen
#' @seealso
#' \url{https://arxiv.org/abs/2007.11919}
divide_conquer_mds <- function(x, l, s, k) {
  initial_row_names <- row.names(x)
  row.names(x) <- 1:nrow(x)

  if (nrow(x) <= l) {
    mds_to_return <- classical_mds(x = x, k = k)
    mds_to_return$eigen <- mds_to_return$eigen / length(mds_to_return$eigen)
  } else {
    index_partition <- get_partitions_for_divide_conquer(n = nrow(x), l = l, s = s, k = k)
    p <- max(index_partition)

    min_len = NA
    eigen <- c()

    # Calculate mds for each partition and take s poits from each subsample
    for (i in 1:p) {

      indexes_current <- which(index_partition == i)
      x_current <- x[indexes_current, ,drop = FALSE]
      row_names_current <- row.names(x_current)

      if (i == 1) {
        
        list_classical_mds <- classical_mds(x = x_current, k = k) 
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

      rn_subsample_previous <- sample(x = row_names_current, size = s, replace = FALSE)

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
