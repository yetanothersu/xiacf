#' Multivariate Xi-Correlogram Matrix
#'
#' Computes the pairwise directional Chatterjee's Xi coefficient for a multivariate
#' time series dataset. It evaluates both "Lead -> Lag" and "Lag -> Lead" relationships
#' across all variable pairs, as well as the Xi-ACF (autocorrelation) for individual variables.
#'
#' @param x A numeric matrix or data.frame containing the multivariate time series (columns = variables).
#' @param max_lag An integer specifying the maximum positive lag to compute. Default is 10.
#' @param n_surr An integer specifying the number of MIAAFT surrogate datasets for hypothesis testing. Default is 399.
#' @param sig_level A numeric value between 0 and 1 specifying the significance level. Default is 0.05.
#' @param max_iter An integer specifying the maximum iterations for the MIAAFT algorithm. Default is 100.
#' @param ... Additional arguments.
#'
#' @return An S3 object of class \code{xi_matrix} containing a tidy data frame of pairwise results.
#' @importFrom stats quantile var
#' @export
xi_matrix <- function(
    x,
    max_lag = 10,
    n_surr = 399,
    sig_level = 0.05,
    max_iter = 100,
    ...
) {
    # --- 1. Strict Input Validation ---
    if (!is.matrix(x) && !is.data.frame(x)) {
        stop("Input 'x' must be a numeric matrix or data.frame.")
    }

    x_mat <- as.matrix(x)

    if (!is.numeric(x_mat)) {
        stop("All columns in 'x' must be numeric.")
    }
    if (any(is.na(x_mat))) {
        stop(
            "Input contains NA values. Please handle missing values before running xi_matrix()."
        )
    }

    n <- nrow(x_mat)
    p <- ncol(x_mat)

    if (p < 2) {
        stop(
            "'x' must have at least 2 columns to compute pairwise correlations."
        )
    }
    if (n < max_lag + 2) {
        stop(sprintf(
            "Time series length (%d) is too short for max_lag = %d.",
            n,
            max_lag
        ))
    }

    # Check zero variance for each column
    vars <- apply(x_mat, 2, stats::var)
    if (any(vars == 0)) {
        bad_cols <- paste(which(vars == 0), collapse = ", ")
        stop(sprintf(
            "The following columns have zero variance (constant values): %s",
            bad_cols
        ))
    }

    if (sig_level <= 0 || sig_level >= 1) {
        stop("'sig_level' must be strictly between 0 and 1 (e.g., 0.05).")
    }

    # --- 2. Check surrogate count for stable Max-Statistic ---
    # Total tests = p * p * max_lag (including ACF on the diagonal)
    num_tests <- p * p * max_lag
    check_surrogate_count(n_surr, sig_level, num_tests)

    # Handle variable names for clean output
    var_names <- colnames(x_mat)
    if (is.null(var_names)) {
        var_names <- paste0("V", 1:p)
    }

    # --- 3. Call the highly optimized C++ engine ---
    # C++ returns: var_lead, var_lag, lag, xi_original, and max_statistic_dist
    res_cpp <- compute_xi_matrix_maxstat_cpp(
        X = x_mat,
        max_lag = as.integer(max_lag),
        n_surr = as.integer(n_surr),
        max_iter = as.integer(max_iter)
    )

    # --- 4. Process Thresholds ---
    global_threshold <- stats::quantile(
        res_cpp$max_statistic_dist,
        probs = 1 - sig_level,
        names = FALSE
    )

    # --- 5. Format Output as a Tidy Data Frame ---
    res_df <- data.frame(
        Lead_Var = var_names[res_cpp$var_lead],
        Lag_Var = var_names[res_cpp$var_lag],
        Lag = res_cpp$lag,
        Xi = res_cpp$xi_original,
        Global_Threshold = rep(global_threshold, length(res_cpp$lag))
    )

    # Calculate excess Xi (clamped at 0 for non-significant values)
    res_df$Xi_Excess <- pmax(0, res_df$Xi - res_df$Global_Threshold)

    # --- 6. Return S3 Object ---
    out <- list(
        data = res_df,
        n = n,
        p = p,
        max_lag = max_lag,
        n_surr = n_surr,
        sig_level = sig_level,
        var_names = var_names
    )
    class(out) <- "xi_matrix"

    return(out)
}

#' Print method for xi_matrix
#' @param x An object of class \code{xi_matrix}.
#' @param ... Additional arguments passed to print.
#' @return The original object \code{x} invisibly.
#' @export
print.xi_matrix <- function(x, ...) {
    cat("\n=== Multivariate Xi-Correlogram Matrix ===\n")
    cat(sprintf(
        "Variables: %d (%s)\n",
        x$p,
        paste(x$var_names, collapse = ", ")
    ))
    cat(sprintf("Time series length: %d\n", x$n))
    cat(sprintf("Max Lag: %d\n", x$max_lag))
    cat(sprintf("Surrogates (MIAAFT): %d\n", x$n_surr))
    cat(sprintf("Significance Level: %g (FWER controlled)\n", x$sig_level))
    cat("==========================================\n")

    # Print top 5 strongest dependencies (based on Excess Xi)
    cat("Top 5 Significant Dependencies:\n")
    top_deps <- x$data[order(x$data$Xi_Excess, decreasing = TRUE), ]
    top_deps <- head(top_deps[top_deps$Xi_Excess > 0, ], 5)

    if (nrow(top_deps) == 0) {
        cat("No significant dependencies found above the global threshold.\n")
    } else {
        print(top_deps, row.names = FALSE)
    }
    cat("\n")
    invisible(x)
}
