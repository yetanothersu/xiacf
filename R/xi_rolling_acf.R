#' Rolling Window Analysis for Xi-ACF
#'
#' @param x A numeric vector representing the time series data.
#' @param time_index An optional vector representing timestamps.
#' @param window_size An integer specifying the size of the rolling window.
#' @param step_size An integer specifying the step size. Default is 1.
#' @param max_lag An integer specifying the maximum lag to compute.
#' @param n_surr An integer specifying the number of IAAFT surrogate datasets. Default is 399.
#' @param sig_level A numeric value specifying the significance level (FWER). Default is 0.05.
#' @param max_iter An integer specifying the maximum iterations for the IAAFT algorithm. Default is 100.
#' @param n_cores An integer specifying the number of cores for parallel execution.
#' @param save_dir A character string specifying the directory path to save intermediate results.
#'
#' @return A \code{data.frame} containing the rolling window results.
#'
#' @importFrom foreach foreach
#' @importFrom doFuture registerDoFuture %dofuture%
#' @importFrom future plan multisession sequential
#' @importFrom progressr progressor with_progress
#' @importFrom dplyr bind_rows
#' @importFrom stats quantile
#' @importFrom parallelly availableCores
#' @export
run_rolling_xi_acf <- function(
    x,
    time_index = NULL,
    window_size,
    step_size = 1,
    max_lag,
    n_surr = 399,
    sig_level = 0.05,
    max_iter = 100,
    n_cores = NULL,
    save_dir = NULL
) {
    # --- Backward Compatibility Check for sig_level ---
    if (sig_level <= 0 || sig_level >= 1) {
        stop("Parameter 'sig_level' must be strictly between 0 and 1.")
    }
    if (sig_level > 0.5) {
        warning(
            sprintf(
                "The interpretation of 'sig_level' has changed from Confidence Level to Significance Level. Your input %g has been automatically converted to %g. Please use sig_level = %g in the future.",
                sig_level,
                1 - sig_level,
                1 - sig_level
            ),
            call. = FALSE
        )
        sig_level <- 1 - sig_level
    }

    n <- length(x)
    if (n < window_size) {
        stop("window_size must be smaller than the length of x.")
    }

    # --- Validate time_index length ---
    if (!is.null(time_index)) {
        if (length(time_index) != length(x)) {
            stop(
                "The length of 'time_index' must exactly match the length of the input time series."
            )
        }
    }

    if (n < window_size) {
        stop("window_size must be smaller than the length of the time series.")
    }

    if (window_size < max_lag + 2) {
        stop("window_size must be strictly greater than max_lag + 1.")
    }

    starts <- seq(1, n - window_size + 1, by = step_size)
    num_windows <- length(starts)

    if (is.null(n_cores)) {
        allowed_cores <- parallelly::availableCores()
        n_cores <- max(1L, allowed_cores - 1L)
    }
    future::plan(future::multisession, workers = n_cores)
    on.exit(future::plan(future::sequential))

    if (!is.null(save_dir) && !dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
    }

    run_with_progress <- function() {
        p <- progressr::progressor(steps = num_windows)

        results_list <- foreach::foreach(
            i = 1:num_windows,
            .options.future = list(
                packages = c("xiacf", "stats"),
                seed = TRUE
            ),
            .errorhandling = "pass"
        ) %dofuture%
            {
                p()

                idx_start <- starts[i]
                idx_end <- idx_start + window_size - 1
                x_win <- x[idx_start:idx_end]

                # Check for zero variance
                if (stats::var(x_win) == 0) {
                    return(NULL)
                }

                # Execute C++ engine for the window
                cpp_res <- compute_xi_acf_maxstat_cpp(
                    x = as.numeric(x_win),
                    max_lag = as.integer(max_lag),
                    n_surr = as.integer(n_surr),
                    max_iter = as.integer(max_iter)
                )

                global_threshold <- stats::quantile(
                    cpp_res$max_statistic_dist,
                    probs = 1 - sig_level,
                    names = FALSE
                )

                df_window <- data.frame(
                    Window_ID = i,
                    Lag = 1:max_lag,
                    Xi = cpp_res$xi_empirical,
                    Global_Threshold = rep(global_threshold, max_lag),
                    stringsAsFactors = FALSE
                )

                df_window$Xi_Excess <- pmax(
                    0,
                    df_window$Xi - df_window$Global_Threshold
                )

                if (!is.null(time_index)) {
                    df_window$Window_Start_Time <- time_index[idx_start]
                    df_window$Window_End_Time <- time_index[idx_end]
                }

                if (!is.null(save_dir)) {
                    saveRDS(
                        df_window,
                        file = file.path(
                            save_dir,
                            sprintf("window_%06d.rds", i)
                        )
                    )
                }

                return(df_window)
            }
        return(results_list)
    }

    new_results_list <- progressr::with_progress(run_with_progress())

    # --- Error Handling & Reporting ---
    error_items <- Filter(
        function(res) inherits(res, "error"),
        new_results_list
    )
    if (length(error_items) > 0) {
        warning(sprintf(
            "[%d rolling windows failed] The computation failed for some windows and they were removed. First error message: %s",
            length(error_items),
            conditionMessage(error_items[[1]])
        ))
    }

    new_results_list <- Filter(is.data.frame, new_results_list)
    final_df <- dplyr::bind_rows(new_results_list)

    return(final_df)
}
