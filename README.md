
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FishMap <img src="man/figures/logo.png" align="right" alt="" width="120" />

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

### Vignettes

You can install the development version of FishMap together with the
vignettes with:

``` r
# install.packages("remotes")
remotes::install_github("balglave/FishMap", build_vignettes = TRUE)
```

Installing with `build_vignettes = TRUE` will allow you to access the
vignettes locally. Once you have installed FishMap with vignettes, you
can access them using the command `vignette(package = "FishMap")`.

## Documentation

Full documentation website on: <https://balglave.github.io/FishMap>  
In particular, you can follow this vignette:

``` r
vignette("user-running-fishmap", package = "FishMap")
```

## Workflow

<img src="man/figures/FishMap_user.png" width="100%" />

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(FishMap)

# Read internal data as Rdata
survey_data_file <- system.file("original_data",
                                "Solea_solea",
                                "survey_data.Rdata",
                                package = "FishMap"
                                )

vmslogbook_data_file <- system.file("original_data",
                                "Solea_solea",
                                "vmslogbook_data.Rdata",
                                package = "FishMap"
                                )

study_domain_file <- system.file("original_data",
                                "Solea_solea",
                                "study_domain.Rdata",
                                package = "FishMap"
                                )


# prepare and load model inputs
fm_data_inputs <- fm_load_data(species = "Solea_solea",
                         fleet = c("OTB_DEF_>=70_0","OTB_CEP_>=70_0","OTT_DEF_>=70_0"),
                         fitted_data = "biomass",
                         survey_data_file = survey_data_file,
                         vmslogbook_data_file = vmslogbook_data_file,
                         study_domain_file = study_domain_file,
                         year_start = 2018,
                         year_end = 2018,
                         month_start = 11,
                         month_end = 11,
                         time_step = "Month",
                         k = 0.25,
                         grid_xmin = -6,
                         grid_xmax = 0,
                         grid_ymin = 42,
                         grid_ymax = 48)

# Fit the model
fm_model_results <- fm_fit_model(fm_data_inputs = fm_data_inputs,
                                 SE = 1,
                                 data_source = 1,
                                 data_obs = 2,
                                 samp_process = 0,
                                 b_constraint = 2,
                                 cov_samp_process = 0,
                                 biomass_temporal = 1,
                                 sampling_temporal = 0,
                                 lf_link = 0,
                                 ref_data = "com",
                                 EM = "est_b",
                                 month_ref = 1)

# Generate figure outputs
fm_generate_graphs(fm_model_results = fm_model_results)
```

## Code of Conduct

Please note that the FishMap project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
