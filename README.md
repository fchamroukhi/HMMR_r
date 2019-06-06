
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Overview

<!-- badges: start -->

<!-- badges: end -->

User-friendly and flexible algorithm for time series **segmentation**
with a Regression model with a Hidden Markov Model Regression (HMMR).

Hidden Markov Model Regression (HMMR) for segmentation of time series
with regime changes. The model assumes that the time series is governed
by a sequence of hidden discrete regimes/states, where each regime/state
has Gaussian regressors as observations. The model parameters are
estimated by MLE via the EM algorithm.

## Installation

You can install the development version of HMMR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fchamroukhi/HMMR_r")
```

To build *vignettes* for examples of usage, type the command below
instead:

``` r
# install.packages("devtools")
devtools::install_github("fchamroukhi/HMMR_r", 
                         build_opts = c("--no-resave-data", "--no-manual"), 
                         build_vignettes = TRUE)
```

Use the following command to display vignettes:

``` r
browseVignettes("HMMR")
```

## Usage

``` r
library(HMMR)

data("simulatedtimeserie")

K <- 5 # number of regimes (states)
p <- 3 # dimension of beta (order of the polynomial regressors)
variance_type <- variance_types$heteroskedastic

n_tries <- 1
max_iter <- 1500
threshold <- 1e-6
verbose <- TRUE

solution <- emHMMR(simulatedtimeserie$X, t(simulatedtimeserie[, 2:ncol(simulatedtimeserie)]), K, p, variance_type, n_tries, max_iter, threshold, verbose)
#> [1] "HMM_regression | EM   : Iteration : 1  Log-likelihood :  -1556.39696825601"
#> [1] "HMM_regression | EM   : Iteration : 2  Log-likelihood :  -1022.47935723687"
#> [1] "HMM_regression | EM   : Iteration : 3  Log-likelihood :  -1019.51830707432"
#> [1] "HMM_regression | EM   : Iteration : 4  Log-likelihood :  -1019.51780361388"

solution$plot()
```

<img src="man/figures/README-unnamed-chunk-5-1.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-2.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-3.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-4.png" style="display: block; margin: auto;" />
