# xiacf 0.5.0

## Major Features & Improvements
* **FWER Control via Max-Statistic:** Implemented rigorous Family-Wise Error Rate (FWER) control for `xi_acf`, `xi_ccf`, and `xi_matrix` using the Max-statistic approach. This prevents data snooping and false positives across multiple lags and variable pairs.
* **C++ Engine Refactoring:** Consolidated and heavily optimized the backend C++ engines to handle simultaneous testing efficiently.
* **Rolling Analysis Modernization:** Updated `run_rolling_xi_acf` and `run_rolling_xi_ccf` to use the modern `%dofuture%` operator and ensured statistically safe parallel RNG with `.options.future = list(seed = TRUE)`.
* **Enhanced Visualizations:** Upgraded all `autoplot` methods to include the new "Global FWER Threshold" (gray ribbon) and to elegantly contrast non-linear dependencies (firebrick triangles) with linear baselines (steelblue points) using the `patchwork` package.
* **Extractor Functions:** Fixed `extract_xi_acf` and `extract_xi_ccf` to correctly populate standard linear confidence intervals (CIs) for seamless integration with plotting methods.

# xiacf 0.4.1

* Fixed a potential C++ RNG synchronization issue by delegating RNG state management entirely to Rcpp::RNGScope.