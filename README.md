
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
#> [1,]   4.880834  -6.517863
#> [2,]  -6.884757 -13.301085
#> [3,]  -2.329795   8.285799
#> [4,] -12.511260   1.877379
#> [5,] -13.757347 -24.448850
#> [6,]   3.904790   3.752449
dc_mds$eigen
#> [1] 6332.318 5301.475

f_mds <- fast_mds(x = x, l = 100, s = 8, k = 2)
head(f_mds$points)
#>             [,1]        [,2]
#> [1,]  -8.4268457  -0.1388874
#> [2,]  -6.2877264 -12.6013831
#> [3,]   6.2530668   1.9676543
#> [4,]  11.1738334  -9.7664017
#> [5,] -13.0800450 -22.8738623
#> [6,]  -0.9024505   2.2110943
f_mds$eigen
#> [1] 4000.572 3179.285

g_mds <- gower_interpolation_mds(x = x, l = 100, k = 2)
head(g_mds$points)
#>            [,1]       [,2]
#> [1,] -4.9745759  -7.049725
#> [2,]  5.7384657 -15.775527
#> [3,]  0.8572678  13.091812
#> [4,] 13.1097118  -1.342500
#> [5,]  9.5546334 -24.905714
#> [6,] -1.6554319   1.444606
g_mds$eigen
#> [1] 122.2643 101.3566
```

With the implementation of *classical MDS*, it is not possible to obtain
a MDS configuration due to computational problems. Try it yourself\!

``` r
library(bigmds)
x <- matrix(data = rnorm(4*10000, sd = 10), nrow = 10000)
dist_matrix <- stats::dist(x = x)
mds_result <- stats::cmdscale(d = dist_matrix, k = 2, eig = TRUE)
```
