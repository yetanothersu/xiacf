#' Bivariate Xi-Cross-Correlation Function
#'
#' Computes the directional Chatterjee's Xi coefficient between two time series
#' across multiple lags, with Family-Wise Error Rate (FWER) control.
#'
#' @param x A numeric vector representing the first time series.
#' @param y A numeric vector representing the second time series.
#' @param max_lag An integer specifying the maximum lag to compute. Default is 20.
#' @param n_surr An integer specifying the number of MIAAFT surrogate datasets. Default is 399.
#' @param sig_level A numeric value specifying the significance level (FWER). Default is 0.05.
#' @param max_iter An integer specifying the maximum iterations for the MIAAFT algorithm. Default is 100.
#' @param direction A character string specifying the testing direction.
#'        "both" computes X->Y and Y->X. "x_leads" computes only X->Y. Default is "both".
#' @param ... Additional arguments.
#'
#' @return An S3 object of class \code{xi_ccf}.
#' @export
xi_ccf <- function(
    x,
    y,
    max_lag = 20,
    n_surr = 399,
    sig_level = 0.05,
    max_iter = 100,
    direction = c("both", "x_leads"),
    ...
) {
    # Parse direction argument
    direction <- match.arg(direction)
    both_directions <- (direction == "both")

    # Extract variable names for plotting
    x_name <- deparse(substitute(x))
    y_name <- deparse(substitute(y))

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

    # Input validation
    if (length(x) != length(y)) {
        stop("Time series 'x' and 'y' must have the same length.")
    }
    if (any(is.na(x)) || any(is.na(y))) {
        stop(
            "Input time series contain NA values. Please handle them before running xi_ccf()."
        )
    }
    if (stats::var(x) == 0 || stats::var(y) == 0) {
        stop("Time series 'x' and 'y' must have non-zero variance.")
    }

    n <- length(x)
    if (n < max_lag + 2) {
        stop(sprintf(
            "Time series length (%d) is too short for max_lag = %d.",
            n,
            max_lag
        ))
    }

    # Dynamic calculation of family size for Dual-Family FWER control
    num_tests_lag0 <- if (both_directions) 2 else 1
    num_tests_lagged <- if (both_directions) 2 * max_lag else max_lag

    # Use the maximum family size between Lag 0 and Lag > 0 for surrogate stability check
    num_tests <- max(num_tests_lag0, num_tests_lagged)

    check_surrogate_count(
        n_surr = n_surr,
        sig_level = sig_level,
        num_tests = num_tests
    )

    # Execute C++ engine
    cpp_res <- compute_xi_ccf_maxstat_cpp(
        x = as.numeric(x),
        y = as.numeric(y),
        max_lag = as.integer(max_lag),
        n_surr = as.integer(n_surr),
        max_iter = as.integer(max_iter),
        both_directions = both_directions
    )

    # Calculate two independent global thresholds
    threshold_lag0 <- stats::quantile(
        cpp_res$max_dist_lag0,
        probs = 1 - sig_level,
        names = FALSE,
        na.rm = TRUE
    )
    threshold_lagged <- stats::quantile(
        cpp_res$max_dist_lagged,
        probs = 1 - sig_level,
        names = FALSE,
        na.rm = TRUE
    )

    # =========================================================================
    # Compute linear CCF and Confidence Interval (CI)
    # Note: In standard stats::ccf(), negative lags represent X->Y,
    # and positive lags represent Y->X.
    # =========================================================================
    ccf_obj <- stats::ccf(
        as.numeric(x),
        as.numeric(y),
        lag.max = max_lag,
        plot = FALSE
    )

    # X leads Y (lag 0 to max_lag): extract and reverse the negative lag side of CCF
    ccf_x_leads <- rev(as.vector(ccf_obj$acf[1:(max_lag + 1)]))
    # Y leads X (lag 0 to max_lag): extract the positive lag side of CCF
    ccf_y_leads <- as.vector(ccf_obj$acf[(max_lag + 1):(2 * max_lag + 1)])

    # Asymptotic confidence interval for CCF
    ccf_ci <- stats::qnorm(1 - sig_level / 2) / sqrt(n)

    # Format output into a tidy data frame
    lag_vec <- 0:max_lag

    if (both_directions) {
        # Bind both directions vertically (Long format)
        res_df <- data.frame(
            Lead_Var = c(
                rep(x_name, length(lag_vec)),
                rep(y_name, length(lag_vec))
            ),
            Lag_Var = c(
                rep(y_name, length(lag_vec)),
                rep(x_name, length(lag_vec))
            ),
            Lag = c(lag_vec, lag_vec),
            Xi = c(cpp_res$xi_emp_x_leads, cpp_res$xi_emp_y_leads),
            CCF = c(ccf_x_leads, ccf_y_leads),
            CCF_CI = rep(ccf_ci, 2 * length(lag_vec)),
            stringsAsFactors = FALSE
        )
    } else {
        # X leads Y only
        res_df <- data.frame(
            Lead_Var = rep(x_name, length(lag_vec)),
            Lag_Var = rep(y_name, length(lag_vec)),
            Lag = lag_vec,
            Xi = cpp_res$xi_emp_x_leads,
            CCF = ccf_x_leads,
            CCF_CI = rep(ccf_ci, length(lag_vec)),
            stringsAsFactors = FALSE
        )
    }

    # Dynamic threshold assignment based on Lag
    is_lag_zero <- res_df$Lag == 0
    res_df$Global_Threshold <- NA_real_
    res_df$Global_Threshold[is_lag_zero] <- threshold_lag0
    res_df$Global_Threshold[!is_lag_zero] <- threshold_lagged

    # Calculate Xi_Excess
    res_df$Xi_Excess <- pmax(0, res_df$Xi - res_df$Global_Threshold)

    # Construct output object
    out <- list(
        data = res_df,
        n = n,
        max_lag = max_lag,
        n_surr = n_surr,
        sig_level = sig_level,
        direction = direction,
        x_name = x_name,
        y_name = y_name
    )

    class(out) <- "xi_ccf"
    return(out)
}

#' Print method for xi_ccf
#' @param x An object of class \code{xi_ccf}.
#' @param ... Additional arguments passed to print.
#' @return The original object \code{x} invisibly.
#' @importFrom utils head
#' @export
print.xi_ccf <- function(x, ...) {
    cat("\n=== Bivariate Xi-Cross-Correlation (CCF) ===\n")
    cat(sprintf("Variables: %s, %s\n", x$x_name, x$y_name))
    cat(sprintf("Time series length: %d\n", x$n))
    cat(sprintf("Max Lag: %d\n", x$max_lag))
    cat(sprintf("Direction: %s\n", x$direction))
    cat(sprintf("Surrogates (MIAAFT): %d\n", x$n_surr))
    cat(sprintf("Significance Level: %g (FWER controlled)\n", x$sig_level))
    cat("============================================\n")

    # Safely extract significant paths
    sig_data <- x$data[
        which(x$data$Xi_Excess > 0),
        c(
            "Lead_Var",
            "Lag_Var",
            "Lag",
            "Xi",
            "CCF",
            "Global_Threshold",
            "Xi_Excess"
        )
    ]

    if (nrow(sig_data) == 0) {
        cat(
            "No significant directional dependencies found above the global threshold.\n"
        )
    } else {
        cat("Top 5 Strongest Causal Pathways:\n")
        print(
            utils::head(sig_data[order(-sig_data$Xi_Excess), ], 5),
            row.names = FALSE
        )
    }
    cat("\n")
    invisible(x)
}
