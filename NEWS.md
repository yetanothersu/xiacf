# xiacf 0.4.0

## New Features
* **Multivariate Network Analysis**: Added `xi_matrix()` to compute directional Chatterjee's Xi and significance thresholds across all pairs in a multivariate dataset simultaneously, powered by a highly optimized n-dimensional MIAAFT C++ engine.
* **Dynamic Significance Levels**: Added the `sig_level` argument to `xi_acf()`, `xi_ccf()`, and `xi_matrix()`, as well as their rolling window counterparts. Users can now easily adjust the hypothesis testing threshold (e.g., 0.95 to 0.99).

## UI/UX Improvements
* **Enhanced Autoplot Aesthetics**: The `autoplot()` methods for all functions have been unified. Statistically significant points (exceeding the dynamic threshold) are now elegantly highlighted with filled red triangles, while non-significant points remain as empty triangles.

## Bug Fixes & Under the Hood
* Resolved NSE (Non-Standard Evaluation) notes in `R CMD check` related to `ggplot2` variables.