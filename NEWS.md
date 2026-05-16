# xiacf 0.5.0

* **FWER Control:** Implemented strict Family-Wise Error Rate ('FWER') control via Max-statistic across all core functions (`xi_acf`, `xi_ccf`, `xi_matrix`).
* **Performance:** Refactored core engines using C++ (RcppArmadillo) for MIAAFT and IAAFT surrogate generation, significantly improving computational speed.
* **Visualization:** Enhanced `autoplot` methods to produce publication-ready `ggplot2` charts utilizing base R `expression()` for native math rendering.
* **Compatibility:** Improved backward compatibility to smoothly transition legacy `sig_level = 0.95` (confidence level) inputs to significance levels, and systematically deprecated older functions with proper warnings.

# xiacf 0.4.1

* Fixed a potential C++ RNG synchronization issue by delegating RNG state management entirely to Rcpp::RNGScope.