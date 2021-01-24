
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
#>            [,1]       [,2]
#> [1,] -20.067096 15.0594659
#> [2,]  -2.466720 -0.3565643
#> [3,]   3.015789 -7.9333540
#> [4,] -13.484393  4.1503936
#> [5,]  -4.579958  0.1910763
#> [6,]   8.640089 17.0066808
dc_mds$eigen
#> [1] 6404.142 5308.511

f_mds <- fast_mds(x = x, l = 100, s = 8, k = 2)
head(f_mds$points)
#>           [,1]       [,2]
#> [1,] 17.654401  4.9166479
#> [2,]  3.478005  0.5466267
#> [3,] -2.976424 -5.9361843
#> [4,]  9.144669 -4.8461553
#> [5,]  6.040185  2.1664650
#> [6,] -5.028190 19.4034958
f_mds$eigen
#> [1] 4012.743 3186.626

g_mds <- gower_interpolation_mds(x = x, l = 100, k = 2)
head(g_mds$points)
#>            [,1]        [,2]
#> [1,] -22.233968  -3.5649938
#> [2,]  -2.229975  -0.9943472
#> [3,]   5.156879   5.1874497
#> [4,] -13.465149   6.0804409
#> [5,]  -5.665802  -2.3634714
#> [6,]   3.667514 -19.0890417
g_mds$eigen
#> [1] 125.2177 109.8476
```

With the implementation of *classical MDS*, it is not possible to obtain
a MDS configuration due to computational problems. Try it yourself\!

``` r
x <- matrix(data = rnorm(4*10000, sd = 10), nrow = 10000)
dist_matrix <- stats::dist(x = x)
mds_result <- stats::cmdscale(d = dist_matrix, k = 2, eig = TRUE)
```
