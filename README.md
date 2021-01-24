
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bigmds

<!-- badges: start -->

<!-- badges: end -->

We present a set of algorithms for *Multidimensional Scaling (MDS)* to
be used with large datasets. *MDS* is a statistic tool for reduction of
dimensionality, using as input a distance matrix of dimensions *n × n*.
When *n* is large, classical algorithms suffer from computational
problems and *MDS* configuration can not be obtained.

With this package, we address these problems by means of three
algorithms:

  - Divide and Conquer MDS
  - Fast MDS
  - MDS based on Gower interpolation

The main idea of these methods is based on partitioning the dataset into
small pieces, where classical methods can work. *Fast MDS* was developed
by *Yang, Tynia & Liu, Jinze & Mcmillan, Leonard & Wang, Wei. (2006)*,
whereas *Divide and Conquer MDS* and *MDS based on Gower interpolation*
were developed by *Delicado and Pachon-Garcia, 2020.*

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
library(bigmds)
x <- matrix(data = rnorm(4*10000, sd = 10), nrow = 10000)

dc_mds <- divide_conquer_mds(x = x, l = 100, num_stitching_points = 8, k = 2)
head(dc_mds$points)
#>            [,1]       [,2]
#> [1,]  -1.997346   9.963247
#> [2,]   2.660519 -14.227134
#> [3,]  11.192389   5.725120
#> [4,] -18.955549  16.739749
#> [5,]  -9.563736  -1.709172
#> [6,]  -7.986035   1.102856
dc_mds$eigen
#> [1] 6338.125 5330.294

f_mds <- fast_mds(x = x, l = 100, s = 8, k = 2)
head(f_mds$points)
#>             [,1]       [,2]
#> [1,]  -0.3099007  3.9377998
#> [2,]   7.0727473 -9.6146201
#> [3,] -14.7220136  0.4099552
#> [4,]   3.8629467 26.3722092
#> [5,]   6.4291939  5.3061731
#> [6,]   5.6960955  4.1458281
f_mds$eigen
#> [1] 4046.674 3144.964

g_mds <- gower_interpolation_mds(x = x, l = 100, k = 2)
head(g_mds$points)
#>            [,1]         [,2]
#> [1,]   1.659173  10.27853361
#> [2,]  -1.741892 -14.10067271
#> [3,]   9.990853   3.81021207
#> [4,] -18.188629  19.90662840
#> [5,]  -7.151691  -0.02155361
#> [6,]  -4.634496   2.49491540
g_mds$eigen
#> [1] 134.4753 110.6267
```

With the implementation of *classical MDS*, it is not possible to obtain
a MDS configuration due to computational problems. Try it yourself\!

``` r
x <- matrix(data = rnorm(4*10000, sd = 10), nrow = 10000)
dist_matrix <- stats::dist(x = x)
mds_result <- stats::cmdscale(d = dist_matrix, k = 2, eig = TRUE)
```
