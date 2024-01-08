
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
#> [1,] -11.434061  4.8926684
#> [2,]  -4.838360  2.4399620
#> [3,]  -7.324822  0.9891919
#> [4,]  -3.620186 -4.0093853
#> [5,]  19.399027  4.6847212
#> [6,]  -1.046009 -8.3715539
landmark_mds_conf$eigen
#> [1] 84.45246 14.90989

interpolation_mds_conf <- interpolation_mds(x = x, l = 200, r = 2, n_cores = 1)
head(interpolation_mds_conf$points)
#>            [,1]       [,2]
#> [1,]  11.394169 -4.8289864
#> [2,]   4.830291 -2.3886329
#> [3,]   7.304446 -0.9396153
#> [4,]   3.606957  3.9901913
#> [5,] -19.316366 -4.7114725
#> [6,]   1.041535  8.3013801
interpolation_mds_conf$eigen
#> [1] 82.87377 19.46736

reduced_mds_conf <- reduced_mds(x = x, l = 200, r = 2, n_cores = 1)
head(reduced_mds_conf$points)
#>             [,1]       [,2]
#> [1,]  11.4971097 -4.5882164
#> [2,]   4.8428897 -2.2911743
#> [3,]   7.2950662 -0.7591714
#> [4,]   3.4539138  4.0999787
#> [5,] -19.3412024 -5.1280972
#> [6,]   0.7574504  8.3214323
reduced_mds_conf$eigen
#> [1] 70.86710 15.85381

pivot_mds_conf <- pivot_mds(x = x, num_pivots = 200, r = 2)
head(pivot_mds_conf$points)
#>            [,1]      [,2]
#> [1,]  10.752574 -5.096018
#> [2,]   4.552134 -2.542650
#> [3,]   6.933063 -1.059457
#> [4,]   3.503738  4.117936
#> [5,] -18.465853 -4.736696
#> [6,]   1.139950  8.631945
pivot_mds_conf$eigen
#> [1] 71.56621 17.08411

divide_mds_conf <- divide_conquer_mds(x = x, l = 200, c_points = 5 * 2, r = 2, n_cores = 1)
head(divide_mds_conf$points)
#>             [,1]       [,2]
#> [1,] -11.7233116  5.0459639
#> [2,]  -6.2774397  2.2368830
#> [3,]  -7.3041934  0.9251586
#> [4,]  -3.0532926 -3.9483349
#> [5,]  20.2630250  5.0752004
#> [6,]  -0.9987763 -7.9887308
divide_mds_conf$eigen
#> [1] 83.94578 16.48651

fast_mds_conf <- fast_mds(x = x, l = 200, s_points = 2*2, r = 2, n_cores = 1)
head(fast_mds_conf$points)
#>             [,1]      [,2]
#> [1,]  11.5501977  4.287541
#> [2,]   5.3589039  2.485099
#> [3,]   5.9829092  1.310702
#> [4,]   3.4068226 -3.990230
#> [5,] -19.3640550  4.825610
#> [6,]   0.6678306 -8.477065
fast_mds_conf$eigen
#> [1] 79.53543 15.93030
```
