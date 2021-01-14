
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

dc_mds <- divide_conquer_mds(x = x, l = 100, s = 8, k = 2)
head(dc_mds$points)
#>            [,1]       [,2]
#> [1,] -11.551444  19.198759
#> [2,]  -3.966821  22.289626
#> [3,]   2.226310  22.271993
#> [4,]  -5.043911  -7.775029
#> [5,]  -2.675354 -18.646642
#> [6,]  -3.138873 -15.864262
dc_mds$eigen
#> [1] 6348.122 5339.048

f_mds <- fast_mds(x = x, l = 100, s = 8, k = 2)
head(f_mds$points)
#>             [,1]       [,2]
#> [1,]   5.8830391  22.912733
#> [2,]  14.8631208  17.788138
#> [3,]  20.6713138  13.208512
#> [4,]   0.2197616   4.891219
#> [5,]  -9.3527003  -4.733618
#> [6,] -16.7064448 -11.082045
f_mds$eigen
#> [1] 4011.772 3203.463

g_mds <- gower_interpolation_mds(x = x, l = 100, k = 2)
head(g_mds$points)
#>             [,1]       [,2]
#> [1,] -12.8725810  18.006569
#> [2,]  -6.6191531  23.155624
#> [3,]  -1.0917355  24.675625
#> [4,]  -6.0497757  -4.674271
#> [5,]  -1.4437130 -17.599594
#> [6,]  -0.5971475 -18.598584
g_mds$eigen
#> [1] 118.1387 105.1294
```

With the implementation of *classical MDS*, it is not possible to obtain
a MDS configuration due to computational problems. Try it yourself\!

``` r
x <- matrix(data = rnorm(4*10000, sd = 10), nrow = 10000)
dist_matrix <- stats::dist(x = x)
mds_result <- stats::cmdscale(d = dist_matrix, k = 2, eig = TRUE)
```
