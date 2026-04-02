#' Rolling Xi-ACF Analysis
#'
#' Performs a rolling window analysis using Chatterjee's Xi coefficient to assess
#' the time-varying non-linear dependence structure of a time series.
#'
#' @param x A numeric vector representing the time series (e.g., log-returns).
#' @param time_index Optional vector of timestamps (e.g., Date, POSIXct) corresponding to x.
#' @param window_size An integer specifying the size of the rolling window.
#' @param step_size An integer specifying the step size by which the window is shifted. Default is 1.
#' @param max_lag An integer specifying the maximum lag to compute Chatterjee's Xi for.
#' @param n_surr An integer specifying the number of surrogate datasets for the null hypothesis test.
#' @param n_cores An integer specifying the number of cores for parallel execution. If \code{NULL}, runs sequentially.
#' @param save_dir A character string specifying the directory path to save intermediate window results as RDS files. If \code{NULL} (default), results are not saved to disk.
#'
#' @return A \code{data.frame} containing the rolling window results, including timestamps if provided.
#'
#' @importFrom foreach foreach
#' @importFrom doFuture registerDoFuture %dofuture%
#' @importFrom future plan multisession sequential
#' @importFrom progressr progressor with_progress
#' @importFrom dplyr bind_rows
#' @importFrom stats quantile sd
#' @export
run_rolling_xi_analysis <- function(
    x,
    time_index = NULL,
    window_size,
    step_size = 1,
    max_lag = 20,
    n_surr = 100,
    n_cores = NULL,
    save_dir = NULL
) {
    # --- 1. Robust Input Validation ---
    n <- length(x)
    if (window_size > n || window_size <= max_lag) {
        stop("Invalid window_size. Must be <= length(x) and > max_lag.")
    }

    # ★ time_index のバリデーション追加
    if (!is.null(time_index)) {
        if (length(time_index) != n) {
            stop("time_index must have the exact same length as x.")
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

    windows_to_process <- setdiff(1:n_windows, completed_windows)

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
    run_with_progress <- function() {
        p <- progressr::progressor(steps = length(windows_to_process))

        results_list <- foreach::foreach(
            i = windows_to_process,
            .errorhandling = 'pass'
        ) %dofuture%
            {
                p()
                idx_start <- start_indices[i]
                idx_end <- idx_start + window_size - 1

                x_window <- x[idx_start:idx_end]

                # If variance is zero, return NULL to skip
                if (stats::sd(x_window) == 0) {
                    return(NULL)
                }

                # ★ Call C++ Engine (最新の関数名に更新)
                res <- compute_xi_acf_iaaft(
                    x_window,
                    max_lag,
                    n_surr
                )

                # 95% surrogate threshold
                xi_threshold <- rep(NA, max_lag)
                if (n_surr > 0 && !is.null(res$xi_surrogates)) {
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
                    # Excess Xi (storing the raw difference)
                    Xi_Excess = pmax(
                        0,
                        as.numeric(res$xi_original) - xi_threshold
                    )
                )

                # ★ Map the actual timestamps if time_index was provided
                if (!is.null(time_index)) {
                    df_window$Window_Start_Time <- time_index[idx_start]
                    df_window$Window_End_Time <- time_index[idx_end]
                }

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

    new_results_list <- progressr::with_progress(run_with_progress())

    # --- 5. Final Result Compilation ---
    new_results_list <- Filter(function(x) is.data.frame(x), new_results_list)
    final_df <- dplyr::bind_rows(new_results_list)

    if (length(completed_windows) > 0 && !is.null(save_dir)) {
        old_files <- file.path(
            save_dir,
            sprintf("window_%06d.rds", completed_windows)
        )
        old_results_list <- lapply(old_files, readRDS)
        old_df <- dplyr::bind_rows(old_results_list)
        final_df <- dplyr::bind_rows(old_df, final_df)
    }

    if (nrow(final_df) > 0) {
        final_df <- final_df[order(final_df$Window_ID, final_df$Lag), ]
    }

    return(final_df)
}
