#' Rolling Directional Xi-CCF Analysis
#'
#' Performs a rolling window analysis using Chatterjee's Xi cross-correlation to assess
#' the time-varying non-linear lead-lag relationship between two time series with FWER control.
#'
#' @param x A numeric vector representing the first time series (predictor/lead candidate).
#' @param y A numeric vector representing the second time series (response/lag candidate).
#' @param time_index Optional vector of timestamps.
#' @param window_size An integer specifying the size of the rolling window.
#' @param step_size An integer specifying the step size. Default is 1.
#' @param max_lag An integer specifying the maximum positive lag to compute.
#' @param n_surr An integer specifying the number of MIAAFT surrogate datasets. Default is 399.
#' @param sig_level A numeric value specifying the significance level (FWER). Default is 0.05.
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
#' @export
run_rolling_xi_ccf <- function(
    x,
    y,
    time_index = NULL,
    window_size,
    step_size = 1,
    max_lag,
    n_surr = 399,
    sig_level = 0.05,
    n_cores = NULL,
    save_dir = NULL
) {
    n <- length(x)
    if (n != length(y)) {
        stop("x and y must have the same length.")
    }
    if (n < window_size) {
        stop("Time series is shorter than the window size.")
    }

    starts <- seq(1, n - window_size + 1, by = step_size)
    n_windows <- length(starts)

    # Create directory for saving intermediate results if specified
    if (!is.null(save_dir)) {
        if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    }

    # Setup parallel backend
    doFuture::registerDoFuture()
    if (!is.null(n_cores)) {
        future::plan(future::multisession, workers = n_cores)
    } else {
        future::plan(future::sequential)
    }

    run_with_progress <- function() {
        p <- progressr::progressor(steps = n_windows)

        results_list <- foreach::foreach(
            i = 1:n_windows,
            .options.future = list(packages = c("xiacf", "stats"), seed = TRUE)
        ) %dofuture%
            {
                p()
                idx_start <- starts[i]
                idx_end <- idx_start + window_size - 1
                x_window <- x[idx_start:idx_end]
                y_window <- y[idx_start:idx_end]

                # Call the latest C++ function (bidirectional FWER control engine)
                cpp_res <- compute_xi_ccf_maxstat_cpp(
                    x = as.numeric(x_window),
                    y = as.numeric(y_window),
                    max_lag = as.integer(max_lag),
                    n_surr = as.integer(n_surr),
                    max_iter = 100L
                )

                # Calculate the global threshold from the Max-statistic distribution
                global_threshold <- stats::quantile(
                    cpp_res$max_statistic_dist,
                    probs = 1 - sig_level,
                    names = FALSE
                )

                num_tests <- 2 * max_lag + 1
                df_window <- data.frame(
                    Window_ID = i,
                    Lag = seq(-max_lag, max_lag, by = 1),
                    Xi = cpp_res$xi_empirical,
                    Global_Threshold = rep(global_threshold, num_tests)
                )

                # Calculate excess Xi above the threshold
                df_window$Xi_Excess <- pmax(
                    0,
                    df_window$Xi - df_window$Global_Threshold
                )

                # Map timestamps if available
                if (!is.null(time_index)) {
                    df_window$Window_Start_Time <- time_index[idx_start]
                    df_window$Window_End_Time <- time_index[idx_end]
                }

                # Save checkpoint
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

    # Execute with progress bar
    new_results_list <- progressr::with_progress(run_with_progress())

    # Filter valid data frames and combine
    new_results_list <- Filter(is.data.frame, new_results_list)
    final_df <- dplyr::bind_rows(new_results_list)

    return(final_df)
}
