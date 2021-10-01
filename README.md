
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bigmds

<!-- badges: start -->

<!-- badges: end -->

*MDS* is a statistic tool for reduction of dimensionality, using as
input a distance matrix of dimensions *n × n*. When *n* is large,
classical algorithms suffer from computational problems and *MDS*
configuration can not be obtained.

With this package, we address these problems by means of three
algorithms:

  - Divide-and-conquer MDS.
  - Fast MDS.
  - Interpolation MDS.

The main idea of these methods is based on partitioning the dataset into
small pieces, where classical methods can work. *Fast MDS* was developed
by *Yang, T., J. Liu, L. McMillan, and W. Wang (2006)*, whereas
*divide-and-conquer MDS* and *interpolation MDS* were developed by
*Delicado P. and C. Pachon-Garcia (2021).*

To obtain more information, please read this
[paper](https://arxiv.org/abs/2007.11919).

## Installation

You can install directly from CRAN with:

``` r
# install.packages("bigmds")
```

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("pachoning/bigmds")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
set.seed(42)
library(bigmds)
x <- matrix(data = rnorm(4*10000), nrow = 10000) %*% diag(c(9, 4, 1, 1))

divide_mds_conf <- divide_conquer_mds(x = x, l = 200, c_points = 2*2, r = 2, n_cores = 1, dist_fn = stats::dist)
head(divide_mds_conf$points)
#>             [,1]       [,2]
#> [1,] -12.0029447  4.5482795
#> [2,]   5.3135571 -0.6207096
#> [3,]  -3.0272576 -1.0857873
#> [4,]  -6.5402649 -1.9113426
#> [5,]  -3.3311073  2.8156667
#> [6,]   0.9705889 -6.5670390
divide_mds_conf$eigen
#> [1] 83.26941 16.27533
divide_mds_conf$GOF
#> [1] 0.9795777 0.9795777

fast_mds_conf <- fast_mds(x = x, l = 200, s_points = 2*2, r = 2, n_cores = 1, dist_fn = stats::dist)
head(fast_mds_conf$points)
#>           [,1]       [,2]
#> [1,] 13.439660  5.0882344
#> [2,] -5.180648 -0.6150152
#> [3,]  3.894922 -0.8759228
#> [4,]  5.248688 -1.6144764
#> [5,]  3.520470  3.1887151
#> [6,] -1.329876 -6.7787889
fast_mds_conf$eigen
#> [1] 81.72223 16.07915
fast_mds_conf$GOF
#> [1] 0.9796994 0.9796994

interpolation_mds_conf <- interpolation_mds(x = x, l = 200, r = 2, n_cores = 1, dist_fn = stats::dist)
head(interpolation_mds_conf$points)
#>             [,1]       [,2]
#> [1,] -12.3616929 -4.6878946
#> [2,]   4.9424093  0.7621167
#> [3,]  -3.3580614  1.1415676
#> [4,]  -5.7834592  1.5567990
#> [5,]  -3.6974408 -2.8075217
#> [6,]   0.8118489  6.4465272
interpolation_mds_conf$eigen
#> [1] 80.06032 16.53691
interpolation_mds_conf$GOF
#> [1] 0.9785652 0.9785652
```

With the implementation of *classical MDS*, it takes much more time to
obtain a MDS configuration due to computational problems. Try it
yourself\!

``` r
x <- matrix(data = rnorm(4*10000, sd = 10), nrow = 10000)
dist_matrix <- stats::dist(x = x)
mds_result <- stats::cmdscale(d = dist_matrix, k = 2, eig = TRUE)
```
