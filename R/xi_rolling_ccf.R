#' Rolling Multivariate Xi-CCF Analysis
#'
#' Performs a rolling window analysis using Chatterjee's Xi cross-correlation to assess
#' the time-varying non-linear lead-lag relationship between two time series.
#'
#' @param x A numeric vector representing the first time series (predictor/lead candidate).
#' @param y A numeric vector representing the second time series (response/lag candidate).
#' @param time_index Optional vector of timestamps (e.g., Date, POSIXct) corresponding to x and y.
#' @param window_size An integer specifying the size of the rolling window.
#' @param step_size An integer specifying the step size by which the window is shifted. Default is 1.
#' @param max_lag An integer specifying the maximum positive lag to compute.
#' @param n_surr An integer specifying the number of MIAAFT surrogate datasets for the null hypothesis test.
#' @param bidirectional Logical. If TRUE (default), computes both "X leads Y" and "Y leads X".
#' @param sig_level A numeric value specifying the significance level for the confidence intervals. Default is 0.95.
#' @param n_cores An integer specifying the number of cores for parallel execution. If \code{NULL}, runs sequentially.
#' @param save_dir A character string specifying the directory path to save intermediate window results as RDS files. If \code{NULL} (default), results are not saved to disk.
#'
#' @return A \code{data.frame} containing the rolling window results in a tidy long-format.
#'
#' @importFrom foreach foreach
#' @importFrom doFuture registerDoFuture %dofuture%
#' @importFrom future plan multisession sequential
#' @importFrom progressr progressor with_progress
#' @importFrom dplyr bind_rows
#' @importFrom stats quantile sd ccf qnorm na.pass
#' @export
run_rolling_xi_ccf <- function(
    x,
    y,
    time_index = NULL,
    window_size,
    step_size = 1,
    max_lag = 20,
    n_surr = 100,
    bidirectional = TRUE,
    sig_level = 0.95,
    n_cores = NULL,
    save_dir = NULL
) {
    # --- 1. Robust Input Validation ---
    if (length(x) != length(y)) {
        stop("Time series 'x' and 'y' must have the exact same length.")
    }
    n <- length(x)
    if (window_size > n || window_size <= 2 * max_lag + 5) {
        stop(
            "Invalid window_size. Must be <= length(x) and reasonably larger than max_lag."
        )
    }
    if (!is.null(time_index)) {
        if (length(time_index) != n) {
            stop("time_index must have the exact same length as x and y.")
        }
    }

    # Calculate starting indices for each window
    start_indices <- seq(1, n - window_size + 1, by = step_size)
    n_windows <- length(start_indices)

    # --- 2. Setup Checkpointing ---
    completed_windows <- integer(0)
    if (!is.null(save_dir)) {
        if (!dir.exists(save_dir)) {
            dir.create(save_dir, recursive = TRUE)
        }
        existing_files <- list.files(
            save_dir,
            pattern = "^window_.*\\.rds$",
            full.names = TRUE
        )
        if (length(existing_files) > 0) {
            # Extract window IDs from filenames
            completed_windows <- as.integer(gsub(
                ".*window_([0-9]+)\\.rds",
                "\\1",
                existing_files
            ))
            message(sprintf(
                "Found %d completed windows in %s. Resuming...",
                length(completed_windows),
                save_dir
            ))
        }
    }

    # Identify windows that still need to be processed
    windows_to_process <- setdiff(1:n_windows, completed_windows)

    # If all windows are already processed, load and return them
    if (length(windows_to_process) == 0) {
        message("All windows are already processed. Loading from disk...")
        all_files <- list.files(
            save_dir,
            pattern = "^window_.*\\.rds$",
            full.names = TRUE
        )
        res_list <- lapply(all_files, readRDS)
        return(dplyr::bind_rows(res_list))
    }

    # --- 3. Parallel Backend Setup ---
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

    # --- 4. Main Rolling Execution ---
    # Define a wrapper for progress reporting
    run_with_progress <- function() {
        p <- progressr::progressor(steps = length(windows_to_process))

        results_list <- foreach::foreach(
            i = windows_to_process,
            .errorhandling = 'pass'
        ) %dofuture%
            {
                p() # Update progress bar
                idx_start <- start_indices[i]
                idx_end <- idx_start + window_size - 1

                x_window <- x[idx_start:idx_end]
                y_window <- y[idx_start:idx_end]

                # If variance is zero, return NULL to skip
                if (stats::sd(x_window) == 0 || stats::sd(y_window) == 0) {
                    return(NULL)
                }

                # Standard CCF Confidence Interval
                ccf_ci <- stats::qnorm((1 + sig_level) / 2) / sqrt(window_size)

                # Call C++ Engine
                res <- compute_xi_ccf_miaaft(
                    x_window,
                    y_window,
                    max_lag,
                    n_surr
                )

                # Helper to calculate 95% threshold
                calc_threshold <- function(surr_matrix) {
                    if (n_surr == 0) {
                        return(rep(NA, length(res$lags)))
                    }
                    apply(surr_matrix, 1, function(r) {
                        stats::quantile(r, sig_level, na.rm = TRUE)
                    })
                }

                # --- Build Long-Format DataFrame for the Window ---

                # (A) Forward: X leads Y
                ccf_fwd_full <- stats::ccf(
                    y_window,
                    x_window,
                    lag.max = max_lag,
                    plot = FALSE,
                    na.action = stats::na.pass
                )
                ccf_fwd <- as.numeric(ccf_fwd_full$acf[
                    (max_lag + 1):(2 * max_lag + 1)
                ])

                df_fwd <- data.frame(
                    Direction = "X leads Y",
                    Lag = as.numeric(res$lags),
                    CCF = ccf_fwd,
                    Xi = as.numeric(res$xi_original_forward),
                    Xi_Threshold = calc_threshold(res$xi_surrogates_forward),
                    CCF_CI = ccf_ci
                )

                # (B) Backward: Y leads X
                if (bidirectional) {
                    ccf_bwd_full <- stats::ccf(
                        x_window,
                        y_window,
                        lag.max = max_lag,
                        plot = FALSE,
                        na.action = stats::na.pass
                    )
                    ccf_bwd <- as.numeric(ccf_bwd_full$acf[
                        (max_lag + 1):(2 * max_lag + 1)
                    ])

                    df_bwd <- data.frame(
                        Direction = "Y leads X",
                        Lag = as.numeric(res$lags),
                        CCF = ccf_bwd,
                        Xi = as.numeric(res$xi_original_backward),
                        Xi_Threshold = calc_threshold(
                            res$xi_surrogates_backward
                        ),
                        CCF_CI = ccf_ci
                    )
                    df_window <- rbind(df_fwd, df_bwd)
                } else {
                    df_window <- df_fwd
                }

                # Add Excess Xi and Window Metadata
                df_window$Xi_Excess <- pmax(
                    0,
                    df_window$Xi - df_window$Xi_Threshold
                )
                df_window$Window_ID <- i
                df_window$Window_Start_Idx <- idx_start
                df_window$Window_End_Idx <- idx_end

                # Map the actual timestamps if time_index was provided
                if (!is.null(time_index)) {
                    df_window$Window_Start_Time <- time_index[idx_start]
                    df_window$Window_End_Time <- time_index[idx_end]
                }

                # Save checkpoint
                if (!is.null(save_dir)) {
                    out_file <- file.path(
                        save_dir,
                        sprintf("window_%06d.rds", i)
                    )
                    saveRDS(df_window, file = out_file)
                }

                return(df_window)
            }
        return(results_list)
    }

    # Execute with progress bar
    new_results_list <- progressr::with_progress(run_with_progress())

    # --- 5. Final Result Compilation ---
    # Filter out potential errors/NULLs from parallel execution
    new_results_list <- Filter(function(x) is.data.frame(x), new_results_list)
    final_df <- dplyr::bind_rows(new_results_list)

    # If resumed from disk, combine newly calculated results with old ones
    if (length(completed_windows) > 0 && !is.null(save_dir)) {
        old_files <- file.path(
            save_dir,
            sprintf("window_%06d.rds", completed_windows)
        )
        old_results_list <- lapply(old_files, readRDS)
        old_df <- dplyr::bind_rows(old_results_list)
        final_df <- dplyr::bind_rows(old_df, final_df)
    }

    # Sort correctly taking Direction into account
    if (nrow(final_df) > 0) {
        final_df <- final_df[
            order(final_df$Window_ID, final_df$Direction, final_df$Lag),
        ]
    }

    return(final_df)
}
