
<!-- README.md is generated from README.Rmd. Please edit that file -->

# qwalkr <a href="https://vitormarquesr.github.io/qwalkr/"><img src="man/figures/logo.png" align="right" height="139" alt="qwalkr website" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/vitormarquesr/qwalkr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vitormarquesr/qwalkr/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## Overview

qwalkr is a numerical suite for investigating quantum walks, providing
estimates of matrices of interests that help you obtain insight into the
evolution of such systems:

**Quantum Walks**

- `ctqwalk()` creates a continuous-time quantum walk.

**Investigate the Hamiltonian**

- `get_eigspace()` obtains the eigenvectors associated with an
  eigenspace.
- `get_eigproj()` obtains the orthogonal projector associated with an
  eigenspace.
- `get_eigschur()`obtains the Schur product of orthogonal projectors.
- `act_eigfun()` applies a function to the Hamiltonian.

**Time Evolution**

- `unitary_matrix()` returns the unitary time evolution operator at a
  given time.
- `mixing_matrix()` returns the mixing matrix at a given time.

**Average Evolution**

- `avg_matrix()` returns the average mixing matrix.
- `gavg_matrix()` returns the generalized average mixing matrix under a
  probability distribution.

## Installation

Currently, qwalkr is only available on GitHub. You can install the
development version of it like so:

``` r
# install.packages("devtools")
devtools::install_github("vitormarquesr/qwalkr")
```

## Usage

``` r
library(qwalkr)

K3 <- rbind(c(0, 1, 1),
            c(1, 0, 1),
            c(1, 1, 0))

w <- ctqwalk(hamiltonian = K3)

w
#> Continuous-Time Quantum Walk
#> 
#> [+]Order: 3 
#> 
#> [+]Spectrum of the Hamiltonian:
#>                   
#> Eigenvalue:   2 -1
#> Multiplicity: 1  2


get_eigproj(w, id=2)
#>            [,1]       [,2]       [,3]
#> [1,]  0.6666667 -0.3333333 -0.3333333
#> [2,] -0.3333333  0.6666667 -0.3333333
#> [3,] -0.3333333 -0.3333333  0.6666667

unitary_matrix(w, t=pi/3)
#>                       [,1]                  [,2]                  [,3]
#> [1,]  0.1666667-0.2886751i -0.3333333+0.5773503i -0.3333333+0.5773503i
#> [2,] -0.3333333+0.5773503i  0.1666667-0.2886751i -0.3333333+0.5773503i
#> [3,] -0.3333333+0.5773503i -0.3333333+0.5773503i  0.1666667-0.2886751i

mixing_matrix(w, t=pi/3)
#>           [,1]      [,2]      [,3]
#> [1,] 0.1111111 0.4444444 0.4444444
#> [2,] 0.4444444 0.1111111 0.4444444
#> [3,] 0.4444444 0.4444444 0.1111111

avg_matrix(w)
#>           [,1]      [,2]      [,3]
#> [1,] 0.5555556 0.2222222 0.2222222
#> [2,] 0.2222222 0.5555556 0.2222222
#> [3,] 0.2222222 0.2222222 0.5555556
```

## Getting Help

For further reference on the usability, check the vignette or the
website of the package.

If you happen to encounter a bug, please file an issue on GitHub.
