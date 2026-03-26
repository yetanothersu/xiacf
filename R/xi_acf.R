#' Xi-ACF Test for Time Series
#'
#' Calculates Chatterjee's Xi and the standard Autocorrelation Function (ACF)
#' along with their respective significance thresholds.
#'
#' @useDynLib xiacf, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats na.omit sd acf quantile qnorm
#'
#' @param x A numeric vector representing the time series data.
#' @param max_lag An integer specifying the maximum number of lags to compute.
#' @param n_surr An integer specifying the number of surrogate datasets to generate for the IAAFT test.
#' @return An object of class \code{"xi_acf"} containing the computed statistics and metadata.
#' @rdname xi_acf
#' @export
xi_acf <- function(x, max_lag = 20, n_surr = 100) {
    # --- 1. Robust Input Validation ---
    if (!is.numeric(x)) {
        stop("Input 'x' must be a numeric vector.")
    }

    if (any(is.na(x))) {
        warning("'x' contains NA values. Removing them before analysis.")
        x <- as.numeric(stats::na.omit(x))
    }

    n <- length(x)
    if (n < 5) {
        stop("Time series length is too short (n < 5).")
    }
    if (stats::sd(x) == 0) {
        stop("Time series has zero variance.")
    }

    # --- 2. Call C++ Engine for Xi Analysis ---
    xi_res <- compute_xi_lags(x, max_lag, n_surr)

    xi_threshold <- rep(NA, max_lag)
    if (n_surr > 0) {
        xi_threshold <- apply(xi_res$xi_surrogates, 1, function(row) {
            stats::quantile(row, 0.95, na.rm = TRUE)
        })
    }

    # --- 3. Compute Standard ACF ---
    acf_res <- stats::acf(
        x,
        lag.max = max_lag,
        plot = FALSE,
        na.action = stats::na.pass
    )
    acf_ci <- stats::qnorm((1 + 0.95) / 2) / sqrt(n)

    # --- 4. Construct S3 Object ---
    df_res <- data.frame(
        Lag = 1:max_lag,
        ACF = as.numeric(acf_res$acf[-1]),
        Xi = as.numeric(xi_res$xi_original),
        Xi_Threshold_95 = xi_threshold,
        ACF_CI = acf_ci
    )

    structure(
        list(data = df_res, n = n, max_lag = max_lag, n_surr = n_surr),
        class = "xi_acf"
    )
}

#' @rdname xi_acf
#' @export
xi_test <- function(x, max_lag = 20, n_surr = 100) {
    # Wrapper function for backward compatibility (alias)
    warning(
        "`xi_test()` is deprecated and will be removed in a future version. Please use `xi_acf()` instead."
    )
    xi_acf(x, max_lag = max_lag, n_surr = n_surr)
}

#' Print method for xi_acf objects
#'
#' @param x An object of class \code{"xi_acf"}.
#' @param ... Additional arguments.
#' @return Invisibly returns the original object.
#' @export
print.xi_acf <- function(x, ...) {
    cat("\n\tChatterjee's Xi-ACF Test\n\n")
    cat("Data length:  ", x$n, "\n")
    cat("Max lag:      ", x$max_lag, "\n")
    cat("Surrogates:   ", x$n_surr, " (IAAFT)\n\n")

    show_cols <- c("Lag", "ACF", "Xi", "Xi_Threshold_95")
    cols_exist <- intersect(show_cols, names(x$data))
    print(x$data[, cols_exist, drop = FALSE], row.names = FALSE)

    invisible(x)
}
