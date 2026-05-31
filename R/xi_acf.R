#' @importFrom Rcpp evalCpp
#' @useDynLib xiacf, .registration = TRUE
NULL

#' Compute empirical Xi-ACF and its significance via IAAFT surrogates
#'
#' @param x A numeric vector representing the time series data. Must not contain missing values (NA) or be a constant.
#' @param max_lag An integer specifying the maximum lag to compute. Default is 10.
#' @param n_surr An integer specifying the number of surrogate datasets to generate. Default is 399.
#' @param sig_level A numeric value between 0 and 1 specifying the significance level. Default is 0.05.
#' @param max_iter An integer specifying the maximum iterations for the IAAFT algorithm. Default is 100.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An object of class \code{xi_acf} containing the empirical ACF,
#'   pointwise thresholds, global threshold, and metadata.
#' @export
xi_acf <- function(
    x,
    max_lag = 10,
    n_surr = 399,
    sig_level = 0.05,
    max_iter = 100,
    ...
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

    if (!is.numeric(x)) {
        stop("Input 'x' must be a numeric vector.")
    }
    if (any(is.na(x))) {
        stop(
            "Input 'x' contains NA values. Please impute or remove them before running xi_acf()."
        )
    }

    n <- length(x)
    if (n < (max_lag + 2)) {
        stop(sprintf(
            "Time series length (%d) is too short for max_lag = %d. Minimum length required is %d.",
            n,
            max_lag,
            max_lag + 2
        ))
    }
    if (stats::var(x) == 0) {
        stop(
            "Input 'x' has zero variance (constant series). Xi coefficient cannot be computed."
        )
    }

    check_surrogate_count(
        n_surr = n_surr,
        sig_level = sig_level,
        num_tests = max_lag
    )

    acf_vals <- as.numeric(stats::acf(x, lag.max = max_lag, plot = FALSE)$acf[
        -1
    ])

    cpp_res <- compute_xi_acf_maxstat_cpp(
        x = as.numeric(x),
        max_lag = as.integer(max_lag),
        n_surr = as.integer(n_surr),
        max_iter = as.integer(max_iter)
    )

    pointwise_threshold <- stats::quantile(
        cpp_res$pointwise_dist,
        probs = 1 - sig_level,
        names = FALSE,
        na.rm = TRUE
    )
    global_threshold <- stats::quantile(
        cpp_res$max_statistic_dist,
        probs = 1 - sig_level,
        names = FALSE,
        na.rm = TRUE
    )

    acf_ci <- stats::qnorm((1 + (1 - sig_level)) / 2) / sqrt(n)

    df_res <- data.frame(
        Lag = 1:max_lag,
        ACF = acf_vals,
        Xi = cpp_res$xi_empirical,
        Pointwise_Threshold = pointwise_threshold,
        Global_Threshold = rep(global_threshold, max_lag),
        ACF_CI = rep(acf_ci, max_lag)
    )

    df_res$Xi_Excess <- pmax(0, df_res$Xi - df_res$Global_Threshold)

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
print.xi_acf <- function(x, ...) {
    cat("\n=== Univariate Xi-Autocorrelation Function ===\n")
    cat(sprintf("Time series length: %d\n", x$n))
    cat(sprintf("Max Lag: %d\n", x$max_lag))
    cat(sprintf("Surrogates (IAAFT): %d\n", x$n_surr))
    cat(sprintf("Significance Level: %g (FWER controlled)\n", x$sig_level))
    cat("==============================================\n")

    sig_data <- x$data[
        x$data$Xi_Excess > 0,
        c("Lag", "Xi", "Global_Threshold", "Xi_Excess")
    ]

    if (nrow(sig_data) == 0) {
        cat(
            "No significant autocorrelations found above the global threshold.\n"
        )
    } else {
        cat("Significant Lags:\n")
        print(sig_data, row.names = FALSE)
    }
    cat("\n")
    invisible(x)
}
