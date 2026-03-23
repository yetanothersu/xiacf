# xiacf: Chatterjee's Rank Correlation for Time Series Analysis

[![CRAN status](https://www.r-pkg.org/badges/version/xiacf)](https://CRAN.R-project.org/package=xiacf)
[![R-CMD-check](https://github.com/[Your-GitHub-Username]/xiacf/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/[Your-GitHub-Username]/xiacf/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
The **xiacf** package provides a robust framework for detecting non-linear, deterministic "If-Then" dependence in time series data. Traditional linear metrics, such as the standard Autocorrelation Function (ACF), often conflate stochastic noise with structural decay. This package overcomes these limitations by integrating **Chatterjee's rank correlation coefficient ($\xi$)** with high-performance C++ backend computations.

## Key Features

* **Fast Computation**: Core engine is written in C++ (`Rcpp` and `RcppArmadillo`) for lightning-fast calculations, even on high-frequency financial data.
* **$\xi$-ACF Test**: Computes Chatterjee's $\xi$ alongside the standard linear ACF for multiple lags (`xi_test()`).
* **Robust Significance Testing**: Automatically generates Iterative Amplitude Adjusted Fourier Transform (IAAFT) surrogate datasets to establish 95% confidence thresholds for non-linear dependence.
* **Rolling Analysis**: Fully parallelized rolling window analysis (`run_rolling_xi_analysis()`) to monitor structural changes over time.
* **Beautiful Visualizations**: Built-in `ggplot2` integration (`autoplot()`) for publication-ready correlograms.

## Installation

You can install the development version of xiacf from [GitHub](https://github.com/) with:

```r
# install.packages("remotes")
remotes::install_github("yetanothersu/xiacf")