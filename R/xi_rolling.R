#' Rolling Xi-ACF Analysis
#'
#' Performs a rolling window analysis using Chatterjee's Xi coefficient to assess
#' the time-varying non-linear dependence structure of a time series.
#'
#' @param x A numeric vector representing the time series (e.g., log-returns).
#' @param window_size An integer specifying the size of the rolling window.
#' @param step_size An integer specifying the step size by which the window is shifted. Default is 1.
#' @param max_lag An integer specifying the maximum lag to compute Chatterjee's Xi for.
#' @param n_surr An integer specifying the number of surrogate datasets for the null hypothesis test.
#' @param n_cores An integer specifying the number of cores for parallel execution. If \code{NULL}, runs sequentially.
#' @param save_dir A character string specifying the directory path to save intermediate window results as RDS files. If \code{NULL} (default), results are not saved to disk.
#'
#' @return A \code{data.frame} containing the rolling window results, including window indices, lags, computed Xi values, surrogate thresholds, and the excess Xi.
#'
#' @importFrom foreach foreach
#' @importFrom doFuture registerDoFuture %dofuture%
#' @importFrom future plan multisession sequential
#' @importFrom progressr progressor
#' @importFrom dplyr bind_rows
#' @importFrom stats quantile sd
#' @export
run_rolling_xi_analysis <- function(
    x,
    window_size,
    step_size = 1,
    max_lag = 20,
    n_surr = 100,
    n_cores = NULL,
    save_dir = NULL
) {
    # --- 1. Input Validation ---
    n_total <- length(x)
    if (window_size > n_total) {
        stop("window_size cannot be larger than the time series length.")
    }

    # Validate core count specification
    if (!is.null(n_cores)) {
        if (!is.numeric(n_cores) || n_cores <= 0 || n_cores %% 1 != 0) {
            stop("n_cores must be a positive integer or NULL.")
        }
    }

    # --- 2. Safe Parallel Setup (Polite Programming) ---
    # Modify the future plan only if the user explicitly specified n_cores,
    # and ensure it is restored upon exit to comply with CRAN policies.
    old_opts <- options(doFuture.rng.onMisuse = "ignore")
    old_plan <- future::plan()

    on.exit(
        {
            options(old_opts)
            future::plan(old_plan)
        },
        add = TRUE
    )

    doFuture::registerDoFuture()
    if (!is.null(n_cores) && n_cores > 1) {
        future::plan(future::multisession, workers = n_cores)
    } else {
        future::plan(future::sequential)
    }

    # --- 3. Prepare Windows ---
    starts <- seq(1, n_total - window_size + 1, by = step_size)
    n_windows <- length(starts)

    # --- 4. Optional Directory Creation for Saving Results ---
    if (!is.null(save_dir)) {
        if (!dir.exists(save_dir)) {
            dir.create(save_dir, recursive = TRUE)
        }
    }

    p <- progressr::progressor(steps = n_windows)

    # --- 5. Execution ---
    results_df <- foreach::foreach(
        i = seq_along(starts),
        .combine = dplyr::bind_rows,
        .options.future = list(seed = TRUE, packages = c("xiacf", "stats")),
        .errorhandling = "remove"
    ) %dofuture%
        {
            p()

            out_file <- NULL
            if (!is.null(save_dir)) {
                out_file <- file.path(save_dir, sprintf("window_%06d.rds", i))
                if (file.exists(out_file)) {
                    return(readRDS(out_file))
                }
            }
            tryCatch(
                {
                    idx_start <- starts[i]
                    idx_end <- idx_start + window_size - 1

                    # Extract the sub-series for the current window
                    y_sub <- x[idx_start:idx_end]

                    # Check for constant series (zero variance)
                    if (stats::sd(y_sub, na.rm = TRUE) == 0) {
                        return(NULL)
                    }

                    # Call the C++ engine
                    res <- compute_xi_acf_iaaft(y_sub, max_lag, n_surr)

                    # Calculate thresholds (NA removal is required since uncomputable lags return NaN in C++)
                    xi_threshold <- rep(NA, max_lag)
                    if (n_surr > 0) {
                        xi_threshold <- apply(
                            res$xi_surrogates,
                            1,
                            function(r) {
                                stats::quantile(r, 0.95, na.rm = TRUE)
                            }
                        )
                    }

                    # Construct the result data frame for the current window
                    df_window <- data.frame(
                        Window_ID = i,
                        Window_Start_Idx = idx_start,
                        Window_End_Idx = idx_end,
                        Lag = 1:max_lag,
                        Xi_Original = as.numeric(res$xi_original),
                        Xi_Threshold_95 = xi_threshold,
                        # Excess Xi (storing the raw difference without flooring to zero)
                        Xi_Excess = pmax(
                            0,
                            as.numeric(res$xi_original) - xi_threshold
                        )
                    )

                    if (!is.null(out_file)) {
                        saveRDS(df_window, file = out_file)
                    }

                    # return
                    return(df_window)
                },
                error = function(e) {
                    # On error, issue a warning and return NULL (which is safely ignored by bind_rows)
                    warning(paste("Error in window", i, ":", e$message))
                    return(NULL)
                }
            )
        }

    return(results_df)
}
