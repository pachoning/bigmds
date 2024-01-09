
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bigmds

<!-- badges: start -->
<!-- badges: end -->

*MDS* is a statistic tool for reduction of dimensionality, using as
input a distance matrix of dimensions *n × n*. When *n* is large,
classical algorithms suffer from computational problems and MDS
configuration can not be obtained. With this package, we address these
problems by means of six algorithms, being two of them original
proposals:

- Landmark MDS proposed by De Silva V. and JB. Tenenbaum (2004).
- Interpolation MDS proposed by Delicado P. and C. Pachón-García (2021)
  \<arXiv:2007.11919\> (original proposal).
- Reduced MDS proposel by Paradis E (2018).
- Pivot MDS prposed by Brandes U. and C. Pich (2007)
- Divide-and-conquer MDS proposed by Delicado P. and C.
  Pachón-García (2021) \<arXiv:2007.11919\> (original proposal).
- Fast MDS, proposed by Yang, T., J. Liu, L. McMillan and W. Wang
  (2006).

To obtain more information, please read this
[paper](https://arxiv.org/abs/2007.11919).

## Installation

You can install the development version of bigmds from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("pachoning/bigmds")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(bigmds)
x <- matrix(data = rnorm(4 * 10000), nrow = 10000) %*% diag(c(9, 4, 1, 1))

landmark_mds_conf <- landmark_mds(x = x, num_landmarks = 200, r = 2)
head(landmark_mds_conf$points)
#>           [,1]       [,2]
#> [1,] -5.365249  0.7337736
#> [2,]  3.306016 -4.5927014
#> [3,] -7.294365  8.3014463
#> [4,] 15.699140  4.6770914
#> [5,] 14.089091  8.1901709
#> [6,]  9.108796  2.6281555
landmark_mds_conf$eigen
#> [1] 93.04856 16.24862

interpolation_mds_conf <- interpolation_mds(x = x, l = 200, r = 2, n_cores = 1)
head(interpolation_mds_conf$points)
#>            [,1]       [,2]
#> [1,]   5.318333  0.7392544
#> [2,]  -3.298290 -4.5693209
#> [3,]   7.265958  8.2536291
#> [4,] -15.606589  4.6378296
#> [5,] -14.014565  8.1367691
#> [6,]  -9.061479  2.6155918
interpolation_mds_conf$eigen
#> [1] 86.61081 16.94258

reduced_mds_conf <- reduced_mds(x = x, l = 200, r = 2, n_cores = 1)
head(reduced_mds_conf$points)
#>           [,1]       [,2]
#> [1,] -5.351389  0.7611677
#> [2,]  3.304342 -4.5616023
#> [3,] -7.302640  8.3011158
#> [4,] 15.673700  4.7006125
#> [5,] 14.069313  8.2312053
#> [6,]  9.110989  2.6164607
reduced_mds_conf$eigen
#> [1] 81.83749 17.31284

pivot_mds_conf <- pivot_mds(x = x, num_pivots = 200, r = 2)
head(pivot_mds_conf$points)
#>            [,1]      [,2]
#> [1,]   5.394381  0.653598
#> [2,]  -3.538906 -4.580966
#> [3,]   7.730890  8.216092
#> [4,] -15.382590  4.823825
#> [5,] -13.580535  8.303857
#> [6,]  -8.941668  2.738167
pivot_mds_conf$eigen
#> [1] 80.62164 16.28707

divide_mds_conf <- divide_conquer_mds(x = x, l = 200, c_points = 5 * 2, r = 2, n_cores = 1)
head(divide_mds_conf$points)
#>            [,1]      [,2]
#> [1,]   5.600940  1.016067
#> [2,]  -3.697173 -4.910147
#> [3,]   5.899926  8.319883
#> [4,] -15.286603  2.617717
#> [5,] -13.690950  8.371893
#> [6,]  -9.687884  2.328566
divide_mds_conf$eigen
#> [1] 84.41986 16.86691

fast_mds_conf <- fast_mds(x = x, l = 200, s_points = 2*2, r = 2, n_cores = 1)
head(fast_mds_conf$points)
#>            [,1]       [,2]
#> [1,]   4.802321  0.5703175
#> [2,]  -3.555341 -4.0507039
#> [3,]   6.766592  8.0851116
#> [4,] -16.051768  4.3472176
#> [5,] -14.437049  8.4028915
#> [6,]  -9.422427  2.4428193
fast_mds_conf$eigen
#> [1] 80.70308 16.11917
```
