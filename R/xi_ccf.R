#' Xi-CCF Test for Multivariate Time Series
#'
#' Calculates Chatterjee's Xi cross-correlation and the standard Cross-Correlation Function (CCF)
#' across positive and negative lags, along with their respective MIAAFT significance thresholds.
#'
#' @useDynLib xiacf, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats na.omit sd ccf quantile qnorm
#'
#' @param x A numeric vector representing the first time series (predictor/lead candidate).
#' @param y A numeric vector representing the second time series (response/lag candidate).
#' @param max_lag An integer specifying the maximum number of lags to compute (computes from -max_lag to +max_lag).
#' @param n_surr An integer specifying the number of MIAAFT surrogate datasets to generate for the null hypothesis test.
#' @return An object of class \code{"xi_ccf"} containing the computed statistics and metadata.
#' @rdname xi_ccf
#' @export
xi_ccf <- function(x, y, max_lag = 20, n_surr = 100) {
    # --- 1. Robust Input Validation ---
    if (!is.numeric(x) || !is.numeric(y)) {
        stop("Inputs 'x' and 'y' must be numeric vectors.")
    }

    if (length(x) != length(y)) {
        stop("Time series 'x' and 'y' must have the exact same length.")
    }

    if (any(is.na(x)) || any(is.na(y))) {
        stop(
            "Inputs contain NA values. Please handle missing values (e.g., imputation or removal) before running xi_ccf."
        )
    }

    n <- length(x)
    if (n < 5) {
        stop("Time series length is too short (n < 5).")
    }
    if (stats::sd(x) == 0 || stats::sd(y) == 0) {
        stop("One or both time series have zero variance.")
    }

    # --- 2. Call C++ Engine for Multivariate Xi Analysis ---
    xi_res <- compute_xi_ccf_cpp(x, y, max_lag, n_surr)

    # Calculate the 95% significance threshold for each lag from the surrogates
    xi_threshold <- apply(
        xi_res$xi_surrogates,
        1,
        function(r) {
            if (n_surr > 0) stats::quantile(r, 0.95, na.rm = TRUE) else NA
        }
    )

    # --- 3. Compute Standard CCF ---
    # Note: R's stats::ccf(x, y) calculates correlation between x[t+k] and y[t].
    # Our C++ convention is x[t] and y[t+k] (Positive lag = x leads y).
    # To perfectly align the lags, we call stats::ccf(y, x).
    ccf_res <- stats::ccf(
        y,
        x,
        lag.max = max_lag,
        plot = FALSE,
        na.action = stats::na.pass
    )
    ccf_ci <- stats::qnorm((1 + 0.95) / 2) / sqrt(n)

    # --- 4. Construct S3 Object ---
    df_res <- data.frame(
        Lag = as.numeric(xi_res$lags),
        CCF = as.numeric(ccf_res$acf),
        Xi = as.numeric(xi_res$xi_original),
        Xi_Threshold_95 = xi_threshold,
        CCF_CI = ccf_ci
    )

    structure(
        list(data = df_res, n = n, max_lag = max_lag, n_surr = n_surr),
        class = "xi_ccf"
    )
}

#' Print method for xi_ccf objects
#'
#' @param x An object of class \code{"xi_ccf"}.
#' @param ... Additional arguments.
#' @return Invisibly returns the original object.
#' @export
print.xi_ccf <- function(x, ...) {
    cat("\n=== Multivariate Chatterjee's Xi Cross-Correlation (CCF) ===\n")
    cat("Time series length (n):", x$n, "\n")
    cat("Maximum lag:", x$max_lag, "\n")
    cat("Number of MIAAFT surrogates:", x$n_surr, "\n")
    cat("\nLags around 0 (Lead-Lag relationship):\n")

    # Show lags around 0 (e.g., -2, -1, 0, 1, 2)
    zero_idx <- which(x$data$Lag == 0)
    show_idx <- max(1, zero_idx - 2):min(nrow(x$data), zero_idx + 2)

    print(
        x$data[show_idx, c("Lag", "CCF", "Xi", "Xi_Threshold_95")],
        row.names = FALSE
    )
    cat("...\n")
    invisible(x)
}
