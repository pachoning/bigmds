get_partitions_for_divide_conquer <- function(n, l, num_stitching_points, k) {

  p <- ceiling(n / (l-num_stitching_points))
  min_sample_size <- max(k + 2, num_stitching_points)

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
#'@description Performs *Multidimensional Scaling* for big datasets. This method can compute a MDS configuration 
#'even when the dataset is so large that classical MDS methods (`cmdscale`) can not be run due to computational 
#'problems.
#'@details In order to obtain a MDS configuration for the entire matrix *x*, it is needed to break the dataset into *p* 
#'submatrices (*divide and conquer strategy*).
#' 
#'In order to obtain *p*, *num_stitching_points* as well as *l* parameters are taken into account. *p* is calculated in
#'such a way that it is possible to use `cmdscale` function in every submatrix.
#'
#'Given a MDS solution, any rotation is another (valid) MDS solution. It means that every partition *1<=k<=p* has its 
#'own coordinate system. 
#'
#'So, in order to keep the same coordinate system, a subsamplig of *num_stitching_points* points are taken from partition 
#'*k-1* and used to align the MDS configuration of partition *k*. Such an alignment is performed by means of Procrustes 
#'transformations. 
#'@param x Dataset.
#'@param l The largest value which allows classical MDS to be computed efficiently, i.e, the larges value which makes 
#'`cmdscale()` be run without any computational issues.
#'@param num_stitching_points Number of stitching points used to align the MDS solutions obtained by the division of *x* 
#'into small groups of matrices. Recommended value: *2·k*.
#'@param k Number of principal coordinates.
#'@param dist_fn Distance function to be used for obtaining a MDS configuration.
#'@return Returns a list containing the following elements:
#' \describe{
#'   \item{points}{A matrix that consists of *k* columns corresponding to the MDS coordinates.}
#'   \item{eigen}{The first *k* eigenvalues.}
#'}
#'@examples
#'x <- matrix(data = rnorm(4*10000, sd = 10), nrow = 10000)
#'cmds <- divide_conquer_mds(x = x, l = 100, num_stitching_points = 8, k = 2, dist_fn = stats::dist)
#'head(cmds$points)
#'cmds$eigen
#'@references
#'Delicado and Pachon-Garcia (2020).
#'\url{https://arxiv.org/abs/2007.11919}
#' 
#'Borg and Groenen (1997). *Modern Multidimensional Scaling*. New York: Springer. pp. 340-342.
#'@export
divide_conquer_mds <- function(x, l, num_stitching_points, k, dist_fn = stats::dist) {

  initial_row_names <- row.names(x)
  row.names(x) <- 1:nrow(x)

  if (nrow(x) <= l) {
    mds_to_return <- classical_mds(x = x, k = k, dist_fn = dist_fn)
    mds_to_return$eigen <- mds_to_return$eigen / length(mds_to_return$eigen)
  } else {
    index_partition <- get_partitions_for_divide_conquer(n = nrow(x), l = l, num_stitching_points = num_stitching_points, k = k)
    p <- max(index_partition)

    min_len = NA
    eigen <- c()

    # Calculate mds for each partition and take num_stitching_points poits from each subsample
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

      rn_subsample_previous <- sample(x = row_names_current, size = num_stitching_points, replace = FALSE)

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
