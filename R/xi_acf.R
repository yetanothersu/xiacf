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
#' @param sig_level A numeric value between 0 and 1 specifying the significance level. Default is 0.95.
#' @return An object of class \code{"xi_acf"} containing the computed statistics and metadata.
#' @rdname xi_acf
#' @export
xi_acf <- function(x, max_lag = 20, n_surr = 100, sig_level = 0.95) {
    if (!is.numeric(x)) {
        stop("Input 'x' must be a numeric vector.")
    }
    if (any(is.na(x))) {
        stop(
            "Inputs contain NA values. Please handle missing values before running."
        )
    }
    if (sig_level <= 0 || sig_level >= 1) {
        stop("'sig_level' must be strictly between 0 and 1.")
    }

    n <- length(x)
    if (n < 5) {
        stop("Time series length is too short (n < 5).")
    }
    if (stats::sd(x) == 0) {
        stop("Time series has zero variance.")
    }

    xi_res <- compute_xi_acf_iaaft(x, max_lag, n_surr)

    xi_threshold <- rep(NA, max_lag)
    if (n_surr > 0) {
        surr_mat <- as.matrix(xi_res$xi_surrogates)
        xi_threshold <- apply(surr_mat, 1, function(row) {
            stats::quantile(row, probs = sig_level, na.rm = TRUE)
        })
    }

    acf_res <- stats::acf(
        x,
        lag.max = max_lag,
        plot = FALSE,
        na.action = stats::na.pass
    )
    # Dynamically link ACF confidence interval to sig_level
    acf_ci <- stats::qnorm((1 + sig_level) / 2) / sqrt(n)

    df_res <- data.frame(
        Lag = 1:max_lag,
        ACF = as.numeric(acf_res$acf[-1]),
        Xi = as.numeric(xi_res$xi_original),
        Xi_Threshold = xi_threshold,
        ACF_CI = acf_ci
    )
    # Flag for statistical significance
    df_res$Xi_Excess <- pmax(0, df_res$Xi - df_res$Xi_Threshold)

    structure(
        list(
            data = df_res,
            n = n,
            max_lag = max_lag,
            n_surr = n_surr,
            sig_level = sig_level
        ),
        class = "xi_acf"
    )
}

#' @rdname xi_acf
#' @export
xi_test <- function(x, max_lag = 20, n_surr = 100) {
    warning("`xi_test()` is deprecated. Please use `xi_acf()` instead.")
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
    cat(sprintf(
        "Significance: %g%% (IAAFT, n_surr = %d)\n\n",
        x$sig_level * 100,
        x$n_surr
    ))

    show_cols <- c("Lag", "ACF", "Xi", "Xi_Threshold", "Xi_Excess")
    cols_exist <- intersect(show_cols, names(x$data))
    print(x$data[, cols_exist, drop = FALSE], row.names = FALSE)
    invisible(x)
}
