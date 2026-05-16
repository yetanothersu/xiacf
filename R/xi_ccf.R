#' Directional Xi-CCF Test for Bivariate Time Series
#'
#' Computes the empirical Cross-Correlation Function (CCF) based on Chatterjee's Xi,
#' and evaluates its statistical significance using MIAAFT surrogates.
#'
#' @param x A numeric vector. Must not contain missing values (NA) or be a constant.
#' @param y A numeric vector of the same length as x. Must not contain missing values or be a constant.
#' @param max_lag An integer specifying the maximum lag. Default is 10.
#' @param n_surr An integer specifying the number of MIAAFT surrogate datasets to generate. Default is 399.
#' @param sig_level A numeric value between 0 and 1 specifying the significance level for FWER control. Default is 0.05.
#' @param max_iter An integer specifying the maximum iterations for the MIAAFT algorithm. Default is 100.
#' @param ... Additional arguments.
#'
#' @return An object of class \code{"xi_ccf"} containing the empirical correlations and Max-statistic thresholds.
#' @rdname xi_ccf
#' @export
xi_ccf <- function(
    x,
    y,
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

    if (!is.numeric(x) || !is.numeric(y)) {
        stop("Inputs 'x' and 'y' must be numeric vectors.")
    }
    if (any(is.na(x)) || any(is.na(y))) {
        stop(
            "Inputs contain NA values. Please remove or impute them before running xi_ccf()."
        )
    }
    n <- length(x)
    if (n != length(y)) {
        stop("Lengths of 'x' and 'y' must be exactly the same.")
    }
    if (n < (max_lag + 2)) {
        stop(sprintf(
            "Time series length (%d) is too short for max_lag = %d.",
            n,
            max_lag
        ))
    }
    if (stats::var(x) == 0 || stats::var(y) == 0) {
        stop(
            "One or both input vectors have zero variance. Xi coefficient cannot be computed."
        )
    }

    num_tests <- 2 * max_lag + 1
    check_surrogate_count <- function(n_surr, sig_level, num_tests) {
        min_required <- ceiling(1 / sig_level) - 1
        if (n_surr < min_required) {
            stop(sprintf(
                "Error: n_surr = %d is too small to calculate the %d%% threshold. Minimum required is %d.",
                n_surr,
                as.integer((1 - sig_level) * 100),
                min_required
            ))
        }
        recommended <- ceiling(num_tests / sig_level)
        if (n_surr < recommended) {
            warning(sprintf(
                "Warning: For %d simultaneous tests at sig_level = %g, the empirical distribution of the max-statistic may be unstable with n_surr = %d. Recommended n_surr is at least %d.",
                num_tests,
                sig_level,
                n_surr,
                recommended
            ))
        }
    }
    check_surrogate_count(n_surr, sig_level, num_tests)

    ccf_obj <- stats::ccf(x, y, lag.max = max_lag, plot = FALSE)
    lags <- as.vector(ccf_obj$lag)
    ccf_vals <- as.vector(ccf_obj$acf)

    cpp_res <- compute_xi_ccf_maxstat_cpp(
        x = as.numeric(x),
        y = as.numeric(y),
        max_lag = as.integer(max_lag),
        n_surr = as.integer(n_surr),
        max_iter = as.integer(max_iter)
    )

    global_threshold <- stats::quantile(
        cpp_res$max_statistic_dist,
        probs = 1 - sig_level,
        names = FALSE
    )

    ccf_ci <- stats::qnorm((1 + (1 - sig_level)) / 2) / sqrt(n)

    df_res <- data.frame(
        Lag = seq(-max_lag, max_lag, by = 1),
        CCF = ccf_vals,
        Xi = cpp_res$xi_empirical,
        Global_Threshold = rep(global_threshold, num_tests),
        CCF_CI = rep(ccf_ci, num_tests)
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
        class = "xi_ccf"
    )
}

#' @rdname xi_ccf
#' @importFrom utils head
#' @export
print.xi_ccf <- function(x, ...) {
    cat("\n=== Bivariate Xi-Cross-Correlation Function ===\n")
    cat(sprintf("Time series length: %d\n", x$n))
    cat(sprintf("Max Lag Range: [-%d, %d]\n", x$max_lag, x$max_lag))
    cat(sprintf("Surrogates (MIAAFT): %d\n", x$n_surr))
    cat(sprintf("Significance Level: %g (FWER controlled)\n", x$sig_level))
    cat("===============================================\n")

    sig_data <- x$data[
        x$data$Xi_Excess > 0,
        c("Lag", "Xi", "Global_Threshold", "Xi_Excess")
    ]

    if (nrow(sig_data) == 0) {
        cat(
            "No significant cross-correlations found above the global threshold.\n"
        )
    } else {
        cat("Top Significant Lead-Lag Relationships:\n")
        print(
            utils::head(sig_data[order(-sig_data$Xi_Excess), ], 5),
            row.names = FALSE
        )
    }
    cat("\n")
    invisible(x)
}
