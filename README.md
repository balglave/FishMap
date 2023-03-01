
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FishMap

<!-- badges: start -->

[![.github/workflows/R-CMD-check](https://github.com/balglave/FishMap/actions/workflows/.github/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/balglave/FishMap/actions/workflows/.github/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/balglave/FishMap/branch/main/graph/badge.svg)](https://app.codecov.io/gh/balglave/FishMap?branch=main)
<!-- badges: end -->

The goal of FishMap is to …

## Installation

### Dependencies

You must [install INLA](https://www.r-inla.org/download-install).
Checkout [INLA doc](https://www.r-inla.org/download-install) to upgrade
INLA.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("graph", "Rgraphviz"), dep=TRUE)

install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
```

### FishMap

You can install the development version of FishMap from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("balglave/FishMap")
```

## Documentation

Full documentation website on: <https://balglave.github.io/FishMap>

## Example

This is a basic example which shows you how to solve a common problem:

``` r
# library(FishMap)
## basic example code
```

## Code of Conduct

Please note that the FishMap project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
