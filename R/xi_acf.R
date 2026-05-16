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
    # 1. Strict Input Validation
    if (!is.numeric(x)) {
        stop("Input 'x' must be a numeric vector.")
    }
    if (any(is.na(x))) {
        stop(
            "Input 'x' contains NA values. Please impute or remove them before running xi_acf()."
        )
    }

    n <- length(x)
    if (n < max_lag + 2) {
        stop(sprintf(
            "Time series length (%d) is too short for max_lag = %d.",
            n,
            max_lag
        ))
    }

    # Check for zero variance (constant series)
    if (stats::var(x) == 0) {
        stop(
            "Input 'x' has zero variance (it is a constant). Cannot compute rank-based correlation."
        )
    }

    if (sig_level <= 0 || sig_level >= 1) {
        stop("'sig_level' must be strictly between 0 and 1 (e.g., 0.05).")
    }

    # 2. Check surrogate count for stable Max-Statistic
    check_surrogate_count(n_surr, sig_level, max_lag)

    # 3. Call the highly optimized C++ engine
    cpp_res <- compute_xi_acf_maxstat_cpp(
        x = as.numeric(x),
        max_lag = as.integer(max_lag),
        n_surr = as.integer(n_surr),
        max_iter = as.integer(max_iter)
    )

    # 4. Process results and calculate thresholds
    empirical_xi <- cpp_res$xi_empirical
    pointwise_dist <- cpp_res$pointwise_dist
    maxstat_dist <- cpp_res$max_statistic_dist

    # Calculate pointwise threshold (per lag)
    pointwise_threshold <- apply(pointwise_dist, 1, function(row) {
        stats::quantile(row, probs = 1 - sig_level, names = FALSE)
    })

    # Calculate GLOBAL threshold for FWER control
    global_threshold <- stats::quantile(
        maxstat_dist,
        probs = 1 - sig_level,
        names = FALSE
    )

    # Standard linear ACF confidence interval (for comparison)
    linear_acf <- stats::acf(
        x,
        lag.max = max_lag,
        plot = FALSE,
        na.action = stats::na.pass
    )
    acf_vals <- as.numeric(linear_acf$acf)[-1] # Remove lag 0
    acf_ci <- stats::qnorm(1 - sig_level / 2) / sqrt(n)

    # Prepare output data frame
    df_res <- data.frame(
        Lag = 1:max_lag,
        ACF = acf_vals,
        Xi = cpp_res$xi_empirical,
        Pointwise_Threshold = pointwise_threshold,
        Global_Threshold = rep(global_threshold, max_lag),
        ACF_CI = rep(acf_ci, max_lag)
    )

    # Flag for statistical significance (using Global Threshold for FWER control)
    df_res$Xi_Excess <- pmax(0, df_res$Xi - df_res$Global_Threshold)

    # 5. Attach attributes for smart plotting and printing
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

    # Extract only significant lags (Xi_Excess > 0)
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
        # 出力を少し綺麗に整えて表示
        print(sig_data, row.names = FALSE)
    }
    cat("\n")
    invisible(x)
}
