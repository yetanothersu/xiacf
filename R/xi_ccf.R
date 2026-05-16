#' Directional Xi-CCF Test for Bivariate Time Series
#'
#' Computes the empirical Cross-Correlation Function (CCF) based on Chatterjee's Xi,
#' and evaluates its statistical significance using MIAAFT surrogates.
#' To strictly control the Family-Wise Error Rate (FWER) and prevent data snooping,
#' this function systematically evaluates all lead and lag relationships
#' simultaneously up to the specified maximum lag.
#'
#' @param x A numeric vector. Must not contain missing values (NA) or be a constant.
#' @param y A numeric vector of the same length as x. Must not contain missing values or be a constant.
#' @param max_lag An integer specifying the maximum lag (evaluates both forward and backward). Default is 10.
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
    # 1. Strict Input Validation
    if (!is.numeric(x) || !is.numeric(y)) {
        stop("Both 'x' and 'y' must be numeric vectors.")
    }
    if (length(x) != length(y)) {
        stop("Lengths of 'x' and 'y' must be exactly the same.")
    }
    if (any(is.na(x)) || any(is.na(y))) {
        stop(
            "Inputs contain NA values. Please remove them before running xi_ccf()."
        )
    }
    if (stats::var(x) == 0 || stats::var(y) == 0) {
        stop("One or both inputs have zero variance (constant vector).")
    }
    if (sig_level <= 0 || sig_level >= 1) {
        stop("'sig_level' must be strictly between 0 and 1.")
    }

    n <- length(x)
    if (n < max_lag + 2) {
        stop(sprintf(
            "Time series length (%d) is too short for max_lag = %d.",
            n,
            max_lag
        ))
    }

    # 2. Check surrogate count (Total tests = 2 * max_lag + 1 for CCF)
    num_tests <- 2 * max_lag + 1
    check_surrogate_count(n_surr, sig_level, num_tests)

    # 3. Call the C++ Engine
    cpp_res <- compute_xi_ccf_maxstat_cpp(
        x = as.numeric(x),
        y = as.numeric(y),
        max_lag = as.integer(max_lag),
        n_surr = as.integer(n_surr),
        max_iter = as.integer(max_iter)
    )

    # 4. Process Thresholds
    pointwise_threshold <- apply(cpp_res$pointwise_dist, 1, function(row) {
        stats::quantile(row, probs = 1 - sig_level, names = FALSE)
    })
    global_threshold <- stats::quantile(
        cpp_res$max_statistic_dist,
        probs = 1 - sig_level,
        names = FALSE
    )

    # Standard linear CCF calculation
    linear_ccf <- stats::ccf(
        x,
        y,
        lag.max = max_lag,
        plot = FALSE,
        na.action = stats::na.pass
    )
    ccf_vals <- as.numeric(linear_ccf$acf)
    ccf_ci <- stats::qnorm(1 - sig_level / 2) / sqrt(n)

    # Prepare output
    df_res <- data.frame(
        Lag = seq(-max_lag, max_lag, by = 1),
        CCF = ccf_vals,
        Xi = cpp_res$xi_empirical,
        Pointwise_Threshold = pointwise_threshold,
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

    # Extract significant relationships
    sig_data <- x$data[
        x$data$Xi_Excess > 0,
        c("Lag", "Xi", "Global_Threshold", "Xi_Excess")
    ]

    if (nrow(sig_data) == 0) {
        cat(
            "No significant cross-correlations found above the global threshold.\n"
        )
    } else {
        # Sort in descending order of impact (Xi_Excess) and show the top 5
        cat("Top Significant Lead-Lag Relationships:\n")
        top_sig <- sig_data[order(sig_data$Xi_Excess, decreasing = TRUE), ]
        print(head(top_sig, 5), row.names = FALSE)

        if (nrow(sig_data) > 5) {
            cat(sprintf(
                "... and %d other significant lags.\n",
                nrow(sig_data) - 5
            ))
        }
    }
    cat("\n")
    invisible(x)
}
