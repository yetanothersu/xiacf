#' Directional Xi-CCF Test for Multivariate Time Series
#'
#' Calculates Chatterjee's Xi cross-correlation and the standard Cross-Correlation Function (CCF)
#' across positive lags to evaluate directional lead-lag relationships.
#'
#' @useDynLib xiacf, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats na.omit sd ccf quantile qnorm
#'
#' @param x A numeric vector representing the first time series (potential cause / lead).
#' @param y A numeric vector representing the second time series (potential effect / lag).
#' @param max_lag An integer specifying the maximum positive lag to compute.
#' @param n_surr An integer specifying the number of MIAAFT surrogate datasets.
#' @param bidirectional Logical. If TRUE (default), computes both "X leads Y" and "Y leads X" using the same surrogates for zero extra computational cost.
#' @return An object of class \code{"xi_ccf"} containing the computed statistics in a tidy long-format.
#' @rdname xi_ccf
#' @export
xi_ccf <- function(x, y, max_lag = 20, n_surr = 100, bidirectional = TRUE) {
    # --- 1. Robust Input Validation ---
    if (!is.numeric(x) || !is.numeric(y)) {
        stop("Inputs 'x' and 'y' must be numeric vectors.")
    }

    if (length(x) != length(y)) {
        stop("Time series 'x' and 'y' must have the exact same length.")
    }

    if (any(is.na(x)) || any(is.na(y))) {
        stop(
            "Inputs contain NA values. Please handle missing values (e.g., via interpolation or na.omit) before running xi_ccf() to ensure data alignment integrity."
        )
    }

    if (stats::var(x) == 0 || stats::var(y) == 0) {
        stop(
            "Time series 'x' or 'y' has zero variance. Correlation cannot be computed."
        )
    }

    n <- length(x)

    # --- 2. C++ Engine Call (The magic of simultaneous bidirectional computation) ---
    # NOTE: Function name updated to match the new naming convention!
    xi_res <- compute_xi_ccf_miaaft(x, y, max_lag, n_surr)

    # Helper function to calculate the 95% threshold from surrogate matrix rows
    calc_threshold <- function(surr_matrix) {
        apply(surr_matrix, 1, function(row) {
            stats::quantile(row, probs = 0.95, na.rm = TRUE)
        })
    }

    # Standard CCF Confidence Interval
    ccf_ci <- stats::qnorm((1 + 0.95) / 2) / sqrt(n)

    # --- 3. Build Long-Format DataFrames ---

    # (A) Forward: X leads Y
    # In base R stats::ccf(y, x), positive lags measure correlation between y[t+k] and x[t], meaning x leads y.
    ccf_fwd_full <- stats::ccf(
        y,
        x,
        lag.max = max_lag,
        plot = FALSE,
        na.action = stats::na.pass
    )
    ccf_fwd <- as.numeric(ccf_fwd_full$acf[(max_lag + 1):(2 * max_lag + 1)])

    df_fwd <- data.frame(
        Direction = "X leads Y",
        Lag = as.numeric(xi_res$lags),
        CCF = ccf_fwd,
        Xi = as.numeric(xi_res$xi_original_forward),
        Xi_Threshold_95 = calc_threshold(xi_res$xi_surrogates_forward),
        CCF_CI = ccf_ci
    )

    # (B) Backward: Y leads X
    if (bidirectional) {
        # stats::ccf(x, y) positive lags mean y leads x.
        ccf_bwd_full <- stats::ccf(
            x,
            y,
            lag.max = max_lag,
            plot = FALSE,
            na.action = stats::na.pass
        )
        ccf_bwd <- as.numeric(ccf_bwd_full$acf[(max_lag + 1):(2 * max_lag + 1)])

        df_bwd <- data.frame(
            Direction = "Y leads X",
            Lag = as.numeric(xi_res$lags),
            CCF = ccf_bwd,
            Xi = as.numeric(xi_res$xi_original_backward),
            Xi_Threshold_95 = calc_threshold(xi_res$xi_surrogates_backward),
            CCF_CI = ccf_ci
        )

        # Combine into a tidy long-format dataframe
        df_res <- rbind(df_fwd, df_bwd)
    } else {
        df_res <- df_fwd
    }

    # --- 4. Construct S3 Object ---
    structure(
        list(
            data = df_res,
            n = n,
            max_lag = max_lag,
            n_surr = n_surr,
            bidirectional = bidirectional
        ),
        class = "xi_ccf"
    )
}

#' Print method for xi_ccf objects
#'
#' @param x An object of class \code{"xi_ccf"}.
#' @param ... Additional arguments.
#' @return Invisibly returns the original object.
#' @importFrom utils head
#' @export
print.xi_ccf <- function(x, ...) {
    cat("\n=== Directional Chatterjee's Xi Cross-Correlation ===\n")
    cat("Time series length (n):", x$n, "\n")
    cat("Maximum positive lag:  ", x$max_lag, "\n")
    cat("Number of surrogates:  ", x$n_surr, "\n")
    cat("Bidirectional evaluation:", x$bidirectional, "\n\n")

    print(head(x$data, 10), row.names = FALSE)

    if (nrow(x$data) > 10) {
        cat("... (Showing first 10 rows. Use `object$data` to see all)\n")
    }
    cat("\n")

    invisible(x)
}
