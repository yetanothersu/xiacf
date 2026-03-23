
<!-- README.md is generated from README.Rmd. Please edit that file -->

# xiacf: Chatterjee’s Rank Correlation for Time Series Analysis

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/xiacf)](https://CRAN.R-project.org/package=xiacf)
[![R-CMD-check](https://github.com/yetanothersu/xiacf/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yetanothersu/xiacf/actions/workflows/R-CMD-check.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end --> The xiacf package provides a robust framework for
detecting complex non-linear and functional dependence in time series
data. Traditional linear metrics, such as the standard Autocorrelation
Function (ACF), often conflate stochastic noise with structural decay.
This package overcomes these limitations by integrating Chatterjee’s
rank correlation coefficient
(![\\xi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cxi
"\\xi")) with high-performance C++ backend computations. Traditional
linear metrics, such as the standard Autocorrelation Function (ACF),
often conflate stochastic noise with structural decay. This package
overcomes these limitations by integrating **Chatterjee’s rank
correlation coefficient
(![\\xi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cxi
"\\xi"))** with high-performance C++ backend computations.

## Key Features

  - **Fast Computation**: Core engine is written in C++ (`Rcpp` and
    `RcppArmadillo`) for lightning-fast calculations, even on
    high-frequency financial data.
  - **![\\xi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cxi
    "\\xi")-ACF Test**: Computes Chatterjee’s
    ![\\xi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cxi
    "\\xi") alongside the standard linear ACF for multiple lags
    (`xi_test()`).
  - **Robust Significance Testing**: Automatically generates Iterative
    Amplitude Adjusted Fourier Transform (IAAFT) surrogate datasets to
    establish 95% confidence thresholds for non-linear dependence.
  - **Rolling Analysis**: Fully parallelized rolling window analysis
    (`run_rolling_xi_analysis()`) to monitor structural changes over
    time.
  - **Beautiful Visualizations**: Built-in `ggplot2` integration
    (`autoplot()`) for publication-ready correlograms.

## Installation

You can install the development version of xiacf from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("yetanothersu/xiacf")
```

*(Note: CRAN submission is currently pending. Once accepted, you can
install it via `install.packages("xiacf")`)*

## Quick Start

Here is a basic example showing how to compute and visualize the
![\\xi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cxi
"\\xi")-ACF against a standard linear ACF.

``` r
library(xiacf)
library(ggplot2)

# Generate a chaotic Logistic Map: x_{t+1} = r * x_t * (1 - x_t)
set.seed(42)
n <- 500
x <- numeric(n)
x[1] <- 0.1 # Initial condition
r <- 4.0 # Fully chaotic regime

for (t in 1:(n - 1)) {
  x[t + 1] <- r * x[t] * (1 - x[t])
}

# 1. Run the Xi-ACF test
# Computes up to 20 lags with 100 IAAFT surrogates for significance testing
results <- xi_test(x, max_lag = 20, n_surr = 100)

# Print summary
print(results)
#> 
#>  Chatterjee's Xi-ACF Test
#> 
#> Data length:   500 
#> Max lag:       20 
#> Surrogates:    100  (IAAFT)
#> 
#>  Lag          ACF        Xi Xi_Threshold_95
#>    1 -0.094245571 0.9880120      0.04806988
#>    2 -0.002595258 0.9760366      0.04447587
#>    3  0.022361912 0.9523173      0.04603879
#>    4  0.014398212 0.9065301      0.05262789
#>    5 -0.031941140 0.8207033      0.05262117
#> ... with 15 more lags

# 2. Visualize the results
# The autoplot method automatically generates a ggplot2 object
autoplot(results)
```

<div class="figure">

<img src="man/figures/README-xi-acf-test-1.png" alt="A correlogram comparing linear and non-linear dependence." width="100%" />

<p class="caption">

Comparison between standard linear ACF and Chatterjee’s Xi-ACF.

</p>

</div>

## Rolling Window Analysis

For advanced market microstructure or structural break detection, you
can run a rolling
![\\xi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cxi
"\\xi")-ACF analysis. This function supports parallel processing using
the `future` ecosystem.

``` r
# Run rolling analysis with a window size of 100 and step size of 10
rolling_res <- run_rolling_xi_analysis(
  x = x,
  window_size = 100,
  step_size = 10,
  max_lag = 5,
  n_surr = 50,
  n_cores = 2 # Set to NULL for sequential execution
)

head(rolling_res)
#>   Window_ID Window_Start_Idx Window_End_Idx Lag Xi_Original Xi_Threshold_95
#> 1         1                1            100   1   0.9403061      0.07725510
#> 2         1                1            100   2   0.8819119      0.10907529
#> 3         1                1            100   3   0.7720026      0.12075893
#> 4         1                1            100   4   0.5930548      0.11375475
#> 5         1                1            100   5   0.3085106      0.11791888
#> 6         2               11            110   1   0.9403061      0.06684694
#>   Xi_Excess
#> 1 0.8630510
#> 2 0.7728366
#> 3 0.6512436
#> 4 0.4793001
#> 5 0.1905918
#> 6 0.8734592
```

## References

  - Chatterjee, S. (2021). A new coefficient of correlation. *Journal of
    the American Statistical Association*, 116(536), 2009-2022.

## License

This project is licensed under the MIT License - see the LICENSE file
for details.
