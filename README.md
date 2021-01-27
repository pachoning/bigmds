
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

  - Divide and Conquer MDS
  - Fast MDS
  - MDS based on Gower interpolation

The main idea of these methods is based on partitioning the dataset into
small pieces, where classical methods can work. *Fast MDS* was developed
by *Tynia, Y., L. Jinze, M. Leonard, and W. Wei, (2006)*, whereas
*Divide and Conquer MDS* and *MDS based on Gower interpolation* were
developed by *Delicado and Pachon-Garcia, 2020.*

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

divide_mds <- divide_conquer_mds(x = x, l = 200, tie = 2*2, k = 2, dist_fn = stats::dist)
head(cbind(divide_mds$points, x[, 1:2]))
#>            [,1]      [,2]        [,3]       [,4]
#> [1,]   4.440107 15.820560  12.5152440  0.5256734
#> [2,] -26.693325  5.184695 -22.5675594 -3.9144910
#> [3,]  -2.493868 11.504705   0.4493810  0.6594500
#> [4,]  -5.827052 -7.525502  -7.3791648 -4.2781657
#> [5,]  -7.524394 -7.878567  -9.3121468 16.7501551
#> [6,]  -4.900891  1.647671   0.3161213 -3.4711479
var(x)
#>             [,1]        [,2]        [,3]         [,4]
#> [1,] 100.5818680  -1.7505342  0.88875629  -0.59582966
#> [2,]  -1.7505342 100.2481097 -0.72941647   2.29347065
#> [3,]   0.8887563  -0.7294165 97.27741476   0.06190491
#> [4,]  -0.5958297   2.2934706  0.06190491 101.25899778
var(divide_mds$points)
#>             [,1]        [,2]
#> [1,] 112.4372211   0.2722196
#> [2,]   0.2722196 112.0091896

fast_mds <- fast_mds(x = x, l = 200, s = 2*2, k = 2, dist_fn = stats::dist)
head(cbind(fast_mds$points, x[, 1:2]))
#>            [,1]       [,2]        [,3]       [,4]
#> [1,]  8.4226897  13.717244  12.5152440  0.5256734
#> [2,] 20.7389506 -14.542903 -22.5675594 -3.9144910
#> [3,] 10.4851020   5.552435   0.4493810  0.6594500
#> [4,] -2.3324524  -8.585684  -7.3791648 -4.2781657
#> [5,]  0.8381321 -13.321215  -9.3121468 16.7501551
#> [6,]  3.3011076  -1.505433   0.3161213 -3.4711479
var(x)
#>             [,1]        [,2]        [,3]         [,4]
#> [1,] 100.5818680  -1.7505342  0.88875629  -0.59582966
#> [2,]  -1.7505342 100.2481097 -0.72941647   2.29347065
#> [3,]   0.8887563  -0.7294165 97.27741476   0.06190491
#> [4,]  -0.5958297   2.2934706  0.06190491 101.25899778
var(fast_mds$points)
#>              [,1]         [,2]
#> [1,] 111.60852604  -0.02264179
#> [2,]  -0.02264179 112.73061992

gower_mds <- gower_interpolation_mds(x = x, l = 200, k = 2, dist_fn = stats::dist)
head(cbind(gower_mds$points, x[, 1:2]))
#>            [,1]        [,2]        [,3]       [,4]
#> [1,] -4.2237053  -9.0632917  12.5152440  0.5256734
#> [2,] 11.6017817 -23.4996323 -22.5675594 -3.9144910
#> [3,] -4.4930638  -9.5259722   0.4493810  0.6594500
#> [4,]  3.4803793  -0.1826961  -7.3791648 -4.2781657
#> [5,]  0.1445507  14.3240928  -9.3121468 16.7501551
#> [6,]  5.2694817  -6.7434663   0.3161213 -3.4711479
var(x)
#>             [,1]        [,2]        [,3]         [,4]
#> [1,] 100.5818680  -1.7505342  0.88875629  -0.59582966
#> [2,]  -1.7505342 100.2481097 -0.72941647   2.29347065
#> [3,]   0.8887563  -0.7294165 97.27741476   0.06190491
#> [4,]  -0.5958297   2.2934706  0.06190491 101.25899778
var(gower_mds$points)
#>           [,1]     [,2]
#> [1,] 100.09401  2.54281
#> [2,]   2.54281 99.06275
```

With the implementation of *classical MDS*, it is not possible to obtain
a MDS configuration due to computational problems. Try it yourself\!

``` r
x <- matrix(data = rnorm(4*10000, sd = 10), nrow = 10000)
dist_matrix <- stats::dist(x = x)
mds_result <- stats::cmdscale(d = dist_matrix, k = 2, eig = TRUE)
```
