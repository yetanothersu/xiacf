# AI Assistant Instructions for `xiacf` Package Development

Welcome, AI Assistant! Before proposing any code modifications, refactoring, or generating new files for the `xiacf` repository, you MUST read and strictly adhere to the guidelines and specifications in this document. 

## 1. Core Philosophy

The `xiacf` package is designed to measure non-linear autocorrelation and directional causal spillover using Chatterjee's Xi coefficient, strictly controlling the Family-Wise Error Rate (FWER).

| Aspect | Principle |
| :--- | :--- |
| **Primary Goal** | Provide robust, non-parametric FWER control for time-series network and autocorrelation analysis. |
| **Scientific Rigor** | **Pursue scientific precision and statistical rigor above all.** Features and algorithms must be mathematically and econometrically valid. Never compromise statistical correctness for the sake of simpler or more elegant code. |
| **Target User** | Researchers and practitioners conducting Exploratory Data Analysis (EDA). |
| **Trade-offs** | Prioritize User Experience (UX) and analytical convenience **only within the strict boundaries of scientific validity.** |

## 2. Rules for AI Behavior

* **No Destructive Refactoring:** Do NOT delete existing functions, arguments, or data frame columns just to make the code "cleaner" or "shorter."
* **Preserve Baselines:** Never remove linear metrics (Pearson's ACF/CCF and their confidence intervals). Users absolutely need these as a baseline to compare against non-linear signals.
* **Propose Delta Changes:** Keep existing logic and English ROxygen comments intact. Propose modifications only for the requested scope.
* **Ask Before Deleting:** If you believe a feature or column is truly obsolete, you must explain your reasoning and get explicit permission from the user before removing it.

## 3. The "Holy" Specifications (Crucial Features to Maintain)

### A. Immutable Data Structures (Tidy Format)
The output of all core functions must always return a tidy `data.frame` (or `tibble`) containing exactly the following columns. Do not rename, retype, or remove them:

* **`xi_acf` Output Format:**
  * `Lag`: The integer lag value (strictly greater than 0).
  * `ACF`: Linear autocorrelation coefficient (Mandatory for baseline comparison).
  * `Xi`: Chatterjee's Xi value for the given lag.
  * `Pointwise_Threshold`: The pointwise significance threshold.
  * `Global_Threshold`: The FWER-controlled max-statistic threshold.
  * `ACF_CI`: Confidence interval for the linear autocorrelation.
  * `Xi_Excess`: The excess dependence (`pmax(0, Xi - Global_Threshold)`).

* **`xi_ccf` and `xi_matrix` Output Format:**
  * `Lead_Var`: The leading variable name (Candidate cause).
  * `Lag_Var`: The lagging variable name (Candidate effect).
  * `Lag`: The integer lag value.
  * `Xi`: Chatterjee's Xi value.
  * `Global_Threshold`: The FWER-controlled max-statistic threshold.
  * `Xi_Excess`: The excess dependence (`pmax(0, Xi - Global_Threshold)`).
  * `CCF`: Linear cross-correlation coefficient (Mandatory for baseline comparison).
  * `CCF_CI`: Confidence interval for the linear correlation.

### B. Core Functions & Boundaries
* **C++ Engine (`src/testing_engine.cpp`)**:
  * `compute_xi_acf_maxstat_cpp`: Must evaluate single time-series non-linear autocorrelation thresholds under the null hypothesis using proper surrogate/shuffling methods.
  * `compute_xi_matrix_maxstat_cpp`: Must SKIP self-loops (`i == j`) in the FWER max-statistic search to isolate pure cross-edges.
  * `compute_xi_ccf_maxstat_cpp`: Must independently evaluate X -> Y and Y -> X based on the `both_directions` flag. Do not revert to the legacy continuous negative-to-positive lag implementation.
  * **Memory Management:** Surrogate generation (`generate_single_miaaft` or equivalent) must use an **in-place update** structure inside the FWER loop. Do NOT rewrite this to generate massive 3D arrays upfront, as it will cause critical memory overflow.
* **R Wrappers (`R/xi_acf.R`, `R/xi_ccf.R`, `R/xi_matrix.R`)**:
  * `xi_acf` must calculate linear ACF using `stats::acf` to populate the baseline `ACF` and `ACF_CI` columns.
  * `xi_ccf` must support the `direction = c("both", "x_leads")` argument and dynamically scale the family size (`num_tests`) for FWER penalty calculations.
  * **Strict NA Handling:** Do NOT use `na.omit()` or implicit omissions as it distorts time-series lag structures. If `NA` values are present anywhere in the input vectors, throw a hard `stop()`.

### C. Domain-Specific Mathematical Rules
* **Lag 0 Exclusion in `xi_acf`:** Lag 0 represents perfect self-identity and carries no statistical meaning for autocorrelation testing. It MUST be excluded from both the calculation and the FWER family size (`num_tests = max_lag`) to avoid inflating the global threshold unnecessarily.
* **Asymmetry in CCF:** Chatterjee's Xi is strictly asymmetric. Even at Lag 0 (contemporaneous), X -> Y and Y -> X must be treated, stored, and plotted as independent values.
* **Anti p-hacking:** Do not alter the family size calculations used for FWER thresholds (e.g., p(p-1) for matrix analysis, or `max_lag` for ACF analysis) to simplify code. These are mathematically required to maintain statistical validity.

### D. UI & Visualization Rules (Print and Autoplot)
* **`autoplot` Stabilization (DO NOT TOUCH VISUALS WITHOUT PERMISSION):**
  * The current chart layouts, custom themes, color palettes, and ribbon layers (`geom_ribbon` for confidence intervals, `geom_hline`/`geom_line` for significance thresholds) are carefully designed for EDA.
  * Do NOT modify the theme settings, font scaling, axis boundaries, or layer ordering in `R/autoplot.R`.
  * For `autoplot.xi_acf` and `autoplot.xi_ccf`, the visual comparison between linear metrics (with their CIs) and non-linear `Xi` (with its `Global_Threshold`) must remain clear and integrated as currently implemented.
  * For `autoplot.xi_matrix`, the multi-panel matrix grid or network plot layout must preserve its custom mapping aesthetics.
* **`print` Methods (S3 Method Consistency):**
  * Console print output (e.g., `print.xi_matrix`, `print.xi_acf`, `print.xi_ccf`) must cleanly summarize key metadata (e.g., number of series `p`, time-series length `n`, `max_lag`, `n_surr`, and `sig_level`).
  * Do NOT rewrite print methods to output raw data frame dumps or truncate important analytical configurations unless requested.

### E. Parallelization & Reproducibility
* **Parallel Framework (`doFuture`):** The package strictly uses the `doFuture` framework (with `foreach` and `%dofuture%`) for parallel processing (e.g., in rolling analyses). Do NOT replace this with `parallel::mclapply`, `pbapply`, or `furrr` unless explicitly requested. Ensure cross-platform compatibility.
* **Strict RNG Control:** Scientific reproducibility is an absolute requirement.
  * **In R:** When executing parallel loops, always ensure proper random seed handling by explicitly passing `.options.future = list(seed = TRUE)` to `foreach`.
  * **In C++:** Strictly rely on R's random number generator (e.g., via `Rcpp::RNGScope` and R's native C API like `R::runif`). Do NOT use standard C++ `<random>` libraries that bypass R's RNG state. Calling `set.seed()` in the R session must guarantee identical, bit-for-bit reproducible results across all surrogate generations and statistical tests.

**Acknowledgment:** Upon reading this file at the start of a session, reply ONLY with: "Guidelines acknowledged. What are the requirements for the current task?"