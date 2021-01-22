
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

dc_mds <- divide_conquer_mds(x = x, l = 100, s = 8, k = 2)
head(dc_mds$points)
#>             [,1]       [,2]
#> [1,] -12.4731741   9.326620
#> [2,]  20.4632999  -4.628658
#> [3,] -13.1688431  19.977839
#> [4,]  15.0430206 -20.701338
#> [5,]  -1.9341617 -17.129690
#> [6,]   0.7434526   3.108340
dc_mds$eigen
#> [1] 6370.328 5289.808

f_mds <- fast_mds(x = x, l = 100, s = 8, k = 2)
head(f_mds$points)
#>            [,1]       [,2]
#> [1,]   9.869926  13.318321
#> [2,] -18.670312  -9.068459
#> [3,]  13.753474  16.834234
#> [4,]  -4.860855 -25.716088
#> [5,]   6.821159 -15.108753
#> [6,]   6.845598  -2.794526
f_mds$eigen
#> [1] 4022.707 3185.568

g_mds <- gower_interpolation_mds(x = x, l = 100, k = 2)
head(g_mds$points)
#>             [,1]       [,2]
#> [1,]  11.3808773  10.250360
#> [2,] -18.8838869  -6.512792
#> [3,]  11.8943066  21.672535
#> [4,] -11.0437163 -22.356115
#> [5,]   4.9306077 -16.938222
#> [6,]   0.1313901   2.712069
g_mds$eigen
#> [1] 112.4474 106.2335
```

With the implementation of *classical MDS*, it is not possible to obtain
a MDS configuration due to computational problems. Try it yourself\!

``` r
x <- matrix(data = rnorm(4*10000, sd = 10), nrow = 10000)
dist_matrix <- stats::dist(x = x)
mds_result <- stats::cmdscale(d = dist_matrix, k = 2, eig = TRUE)
```
