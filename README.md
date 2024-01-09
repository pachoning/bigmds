
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
#>            [,1]       [,2]
#> [1,] -22.652391 -2.9637642
#> [2,]  -4.028894  0.8557315
#> [3,]  -2.053481  5.0792924
#> [4,] -12.298362  0.1943286
#> [5,]  -7.464274  6.6271097
#> [6,]  -7.440109 -2.4224585
landmark_mds_conf$eigen
#> [1] 80.74139 16.03830

interpolation_mds_conf <- interpolation_mds(x = x, l = 200, r = 2, n_cores = 1)
head(interpolation_mds_conf$points)
#>           [,1]       [,2]
#> [1,] 22.549026 -2.9628124
#> [2,]  4.028638  0.8427597
#> [3,]  2.036448  5.0659826
#> [4,] 12.229514  0.2076230
#> [5,]  7.433366  6.5968027
#> [6,]  7.404454 -2.4194837
interpolation_mds_conf$eigen
#> [1] 99.26074 13.40565

reduced_mds_conf <- reduced_mds(x = x, l = 200, r = 2, n_cores = 1)
head(reduced_mds_conf$points)
#>            [,1]       [,2]
#> [1,] -22.643870 -2.9838078
#> [2,]  -4.040102  0.7672776
#> [3,]  -2.063242  5.0807164
#> [4,] -12.311794  0.1926956
#> [5,]  -7.475948  6.5787120
#> [6,]  -7.428051 -2.4057501
reduced_mds_conf$eigen
#> [1] 84.92534 18.47403

pivot_mds_conf <- pivot_mds(x = x, num_pivots = 200, r = 2)
head(pivot_mds_conf$points)
#>           [,1]       [,2]
#> [1,] 23.735663  2.8842654
#> [2,]  4.231575 -0.7958750
#> [3,]  2.191797 -4.9042669
#> [4,] 12.918229 -0.1648954
#> [5,]  7.871017 -6.3850176
#> [6,]  7.781056  2.3363735
pivot_mds_conf$eigen
#> [1] 89.25984 14.78125

divide_mds_conf <- divide_conquer_mds(x = x, l = 200, c_points = 5 * 2, r = 2, n_cores = 1)
head(divide_mds_conf$points)
#>            [,1]       [,2]
#> [1,] -23.297597 -2.9259512
#> [2,]  -2.400426  0.7265173
#> [3,]  -3.290848  4.6427289
#> [4,] -11.887703  0.2840322
#> [5,]  -7.756640  6.8074297
#> [6,]  -7.920026 -2.3979162
divide_mds_conf$eigen
#> [1] 87.52543 16.17188

fast_mds_conf <- fast_mds(x = x, l = 200, s_points = 2*2, r = 2, n_cores = 1)
head(fast_mds_conf$points)
#>           [,1]       [,2]
#> [1,] 21.358785  2.7025729
#> [2,]  2.745327 -0.9114227
#> [3,]  1.094829 -5.2416839
#> [4,] 12.034316 -0.3889722
#> [5,]  7.691530 -6.3699432
#> [6,]  8.328885  2.1651105
fast_mds_conf$eigen
#> [1] 80.82792 15.63501
```
