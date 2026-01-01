#' Rolling Xi-ACF Analysis
#'
#' Performs a rolling window analysis using Chatterjee's Xi coefficient
#' and IAAFT surrogates to detect non-linear dependencies.
#'
#' @useDynLib xiacf, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @param ts_vec A numeric vector of the time series (e.g., log-returns).
#' @param window_size Integer. The size of the rolling window.
#' @param step_size Integer. The step size for moving the window.
#' @param max_lag Integer. Maximum lag to calculate Xi for.
#' @param n_surr Integer. Number of surrogates for the null hypothesis test.
#' @param n_cores Integer (optional). Number of cores for parallel processing.
#'   If NULL, detects automatically.
#'
#' @return A data.frame containing:
#'   \item{Window_Start_Idx}{Start index of the window}
#'   \item{Lag}{Lag (tau)}
#'   \item{Xi_Original}{Observed Xi value}
#'   \item{Xi_Threshold_95}{95th percentile of surrogate Xi values}
#'   \item{Xi_Excess}{Excess amount (Xi_Original - Xi_Threshold_95)}
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom doFuture registerDoFuture
#' @importFrom future plan multisession sequential
#' @importFrom progressr progressor with_progress
#' @importFrom dplyr bind_rows mutate
#' @export
run_rolling_xi_analysis <- function(
    ts_vec,
    window_size,
    step_size,
    max_lag,
    n_surr = 100,
    n_cores = NULL
) {
    # --- 1. Parallel Setup ---
    ov_plan <- future::plan()

    on.exit(
        {
            future::plan(ov_plan) # futureプランを元に戻す
            foreach::registerDoSEQ() # foreachをデフォルト（直列）に戻す
        },
        add = TRUE
    )

    if (is.null(n_cores)) {
        use_cores <- max(1, parallel::detectCores(logical = FALSE) - 2)
    } else {
        use_cores <- n_cores
    }

    if (use_cores > 1) {
        message(sprintf("Running on %d cores...", use_cores))
        future::plan(future::multicore, workers = use_cores)
        doFuture::registerDoFuture()
    } else {
        future::plan(future::sequential)
    }

    # --- 2. Prepare Windows ---
    n_total <- length(ts_vec)
    starts <- seq(1, n_total - window_size + 1, by = step_size)
    n_windows <- length(starts)

    # Progress bar
    p <- progressr::progressor(steps = n_windows)

    # --- 3. Execution ---
    # Note: Requires internal C++ function `run_xi_test_cpp` to be exported
    results_df <- foreach::foreach(
        i = seq_along(starts),
        .combine = dplyr::bind_rows,
        .packages = c("dplyr", "xiacf"), # 自分自身をロード
        .options.future = list(seed = TRUE) # Robust RNG
    ) %dopar%
        {
            p()

            idx_start <- starts[i]
            idx_end <- idx_start + window_size - 1
            y_sub <- ts_vec[idx_start:idx_end]

            # Call C++ Engine
            res <- run_xi_test_cpp(y_sub, max_lag, n_surr)

            # Calculate Threshold
            null_95 <- apply(res$xi_surrogates, 1, function(x) {
                quantile(x, 0.95)
            })

            data.frame(
                Window_Start_Idx = idx_start,
                Window_End_Idx = idx_end,
                Lag = 1:max_lag,
                Xi_Original = res$xi_original,
                Xi_Threshold_95 = null_95
            ) %>%
                dplyr::mutate(
                    Is_Significant = Xi_Original > Xi_Threshold_95,
                    Xi_Excess = pmax(0, Xi_Original - Xi_Threshold_95)
                )
        }

    return(results_df)
}
