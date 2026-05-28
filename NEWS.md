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