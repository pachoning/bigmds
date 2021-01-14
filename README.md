
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
#>            [,1]        [,2]
#> [1,]  -3.452944  -2.1672275
#> [2,]  -3.071801  -0.8931906
#> [3,] -10.445148  -5.3594996
#> [4,]   8.399469   0.1979912
#> [5,]  24.613314 -11.4814874
#> [6,] -10.013429   4.1016702
dc_mds$eigen
#> [1] 6367.130 5287.151

f_mds <- fast_mds(x = x, l = 100, s = 8, k = 2)
head(f_mds$points)
#>            [,1]       [,2]
#> [1,]   4.771132  2.7682618
#> [2,]   3.262440 -0.8731155
#> [3,]   9.649698 -5.1563154
#> [4,]  -7.678745  4.3893413
#> [5,] -11.095901 22.7803413
#> [6,]   3.209353 -7.1795478
f_mds$eigen
#> [1] 3992.765 3204.905

g_mds <- gower_interpolation_mds(x = x, l = 100, k = 2)
head(g_mds$points)
#>            [,1]        [,2]
#> [1,]   1.179854  -1.3243016
#> [2,]   2.945878  -0.8688787
#> [3,]  12.567168  -6.2436956
#> [4,]  -7.887559  -0.3552226
#> [5,] -23.318461 -12.2386295
#> [6,]   8.234146   4.4834985
g_mds$eigen
#> [1] 109.7817 104.8187
```

With the implementation of *classical MDS*, it is not possible to abtain
a MDS configuration due to computation problems. Try it yourself\!

``` r
library(bigmds)
x <- matrix(data = rnorm(4*10000, sd = 10), nrow = 10000)
dist_matrix <- stats::dist(x = x)
mds_result <- stats::cmdscale(d = dist_matrix, k = 2, eig = TRUE)
```
