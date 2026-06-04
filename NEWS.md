# xiacf 0.6.4
* Fixed a severe memory bottleneck by transforming `xi_coefficient()` arguments into const references, eliminating millions of redundant matrix deep copies during intensive loops.
* Implemented memory-safe, in-place surrogate generation for univariate IAAFT (`generate_single_iaaft_worker`), preventing OOM-killer crashes in Linux parallel worker environments (`future.apply`).
* Enhanced `run_rolling_xi_acf` and `run_rolling_xi_ccf` with strict input validation for `time_index` and explicit `warning()` tracking for failed windows.
* Removed the memory-hazardous `surrogate_miaaft_cpp` C++ backend and replaced it with a safe R-level deprecated wrapper.
 
# xiacf 0.6.3
* Fixed a critical bug where massive Monte Carlo simulations could cause a segmentation fault (exit code -11) due to floating-point rounding errors in C++.
* Ensured MIAAFT surrogate generation strictly respects R's global RNG state (e.g., `set.seed()`), enabling bit-for-bit reproducibility in parallel environments.

# xiacf 0.6.2

## CRAN Policy Compliance
* **Strict CPU Core Limits:** Updated the default core allocation in rolling functions (`run_rolling_xi_acf`, `run_rolling_xi_ccf`) to strictly comply with CRAN server policies (e.g., respecting `MC_CORES` and `OMP_THREAD_LIMIT`). Replaced `parallel::detectCores()` with `max(1L, parallelly::availableCores() - 1L)` to ensure safe execution without exceeding permitted limits or causing 0-core crashes.
* **Sequential Testing:** Enforced strict sequential execution (`future::plan(sequential)`) via `setup.R` across the entire `testthat` suite to adhere to CRAN's 2-core testing limit.

## Bug Fixes & Refactoring
* Fixed a bug in `xi_matrix()` where variable names (`var_names`) were inadvertently omitted from the output object, resulting in empty variable lists during `print()`.
* Unified the surrogate count FWER validation in `xi_ccf()` to utilize the centralized internal helper `check_surrogate_count()`.
* Cleaned up redundant `"Error:"` and `"Warning:"` prefixes in custom message strings across the package to ensure native R console formatting.
* Removed unreachable dead code regarding duplicate `sig_level` validation in `xi_matrix()`.

## Enhancements
* Exposed the `max_iter` argument in `run_rolling_xi_acf()` and `run_rolling_xi_ccf()`, allowing users to explicitly tune the maximum number of iterations for the IAAFT/MIAAFT surrogate convergence.
* Simplified and standardized the linear ACF confidence interval mathematical formula in `xi_acf()`.

# xiacf 0.6.1

## Major Changes and Bug Fixes
* **Core Math Fix:** Corrected the denominator in the internal C++ implementation of Chatterjee's Xi coefficient to exactly $n^2 - 1$. This fixes an issue where the baseline for independent variables was slightly inflated (approx. 0.33) rather than asymptotically converging to 0. (Thanks to user reports and rigorous auditing).
* **Dual-Family FWER Control:** Introduced a mathematically rigorous separation of FWER control families in `xi_ccf` and `xi_matrix`. The Max-Statistic null distribution is now strictly evaluated independently for Contemporaneous (Lag 0) and Temporal (Lag > 0) dependencies. This resolves the "Lag-0 Masking Effect," significantly improving the statistical power to detect delayed causal propagation.
* **Autocorrelation Exclusion:** Excluded completely diagonal elements (autocorrelations) from the FWER search space in `xi_matrix` to prevent arbitrary threshold inflation.

## Minor Improvements
* Improved parallel backend safety using `%dofuture%` with explicit package loading.
* Restored the `Window_ID` tracking column in the outputs of `run_rolling_xi_acf` and `run_rolling_xi_ccf`.
* Upgraded internal surrogate generation checks (`check_surrogate_count`) to rigorously validate user-provided `n_surr` against the dynamic size of FWER test families.

# xiacf 0.6.0

* **Directional Separation:** `xi_ccf()` now explicitly separates causal directions (`direction = "both"`, `"x_leads"`, `"y_leads"`) and returns Tidy data frames for easier downstream EDA.
* **Enhanced FWER Control:** Excluded self-loops (autocorrelations) from the `xi_matrix()` Max-Statistic empirical null distribution, restoring statistical power for detecting true cross-edge pathways.
* **Lag 0 Synchronization:** Included contemporaneous effects (Lag 0) in `xi_matrix()` C++ calculations to synchronize behavior with `xi_ccf()`.
* **On-demand Thresholds:** Updated extraction functions (`extract_xi_acf()`, `extract_xi_ccf()`) to dynamically recompute exact FWER thresholds using preserved raw data (`data_raw`).
* **Unified Plotting Aesthetics:** Redesigned `autoplot()` methods for `xi_ccf` (vertical faceting) and `xi_matrix` (diagonal variable labels) to provide publication-ready, unified Tidyverse compatibility.
* **Bug Fixes:** Fixed a C++ segmentation fault related to Lag 0 memory allocation, and resolved a Non-Standard Evaluation (NSE) bug in extractors that caused messy facet labels.
* **Test Suite Overhaul:** Expanded and strictly validated the new tidy structures and directional logic (58 tests passing).

# xiacf 0.5.0

* **FWER Control:** Implemented strict Family-Wise Error Rate ('FWER') control via Max-statistic across all core functions (`xi_acf`, `xi_ccf`, `xi_matrix`).
* **Performance:** Refactored core engines using C++ (RcppArmadillo) for MIAAFT and IAAFT surrogate generation, significantly improving computational speed.
* **Visualization:** Enhanced `autoplot` methods to produce publication-ready `ggplot2` charts utilizing base R `expression()` for native math rendering.
* **Compatibility:** Improved backward compatibility to smoothly transition legacy `sig_level = 0.95` (confidence level) inputs to significance levels, and systematically deprecated older functions with proper warnings.

# xiacf 0.4.1

* Fixed a potential C++ RNG synchronization issue by delegating RNG state management entirely to Rcpp::RNGScope.