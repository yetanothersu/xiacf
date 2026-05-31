#' Multivariate Xi-Correlogram Matrix
#'
#' Computes the pairwise directional Chatterjee's Xi coefficient for a multivariate time series.
#'
#' @param x A numeric matrix or data.frame containing the multivariate time series (columns = variables).
#' @param max_lag An integer specifying the maximum positive lag to compute. Default is 10.
#' @param n_surr An integer specifying the number of MIAAFT surrogate datasets. Default is 399.
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

    if (!is.matrix(x) && !is.data.frame(x)) {
        stop("Input 'x' must be a numeric matrix or data.frame.")
    }
    x_mat <- as.matrix(x)

    if (!is.numeric(x_mat)) {
        stop("Input 'x' must be completely numeric.")
    }
    if (any(is.na(x_mat))) {
        stop(
            "Input 'x' contains NA values. Please handle them before running xi_matrix()."
        )
    }

    n <- nrow(x_mat)
    p <- ncol(x_mat)

    if (p < 2) {
        stop(
            "Input 'x' must have at least 2 columns to compute pairwise correlations."
        )
    }
    if (n < (max_lag + 2)) {
        stop(sprintf(
            "Time series length (%d) is too short for max_lag = %d.",
            n,
            max_lag
        ))
    }
    if (sig_level <= 0 || sig_level >= 1) {
        stop("Parameter 'sig_level' must be strictly between 0 and 1.")
    }

    var_names <- colnames(x)
    if (is.null(var_names)) {
        var_names <- paste0("V", 1:p)
    }

    for (j in 1:p) {
        if (stats::var(x_mat[, j]) == 0) {
            stop(sprintf(
                "Variable '%s' has zero variance. Matrix calculations cannot proceed.",
                var_names[j]
            ))
        }
    }

    num_pairs <- p * (p - 1)

    # Family A (Lag 0) has 'num_pairs' tests.
    # Family B (Lag > 0) has 'num_pairs * max_lag' tests.
    # We use the maximum family size to evaluate the required number of surrogates.
    num_tests <- num_pairs * max(1, max_lag)

    check_surrogate_count(
        n_surr = n_surr,
        sig_level = sig_level,
        num_tests = num_tests
    )

    # Execute C++ engine for the entire matrix
    cpp_res <- compute_xi_matrix_maxstat_cpp(
        X = x_mat,
        max_lag = as.integer(max_lag),
        n_surr = as.integer(n_surr),
        max_iter = as.integer(max_iter)
    )

    # Extract results and populate dataframe
    res_list <- list()
    idx <- 1

    for (i in 1:p) {
        for (j in 1:p) {
            for (k in 0:max_lag) {
                res_list[[idx]] <- data.frame(
                    Lead_Var = var_names[i],
                    Lag_Var = var_names[j],
                    Lag = k,
                    Xi = cpp_res$xi_emp[i, j, k + 1],
                    stringsAsFactors = FALSE
                )
                idx <- idx + 1
            }
        }
    }

    res_df <- dplyr::bind_rows(res_list)

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

    # Identify row attributes
    is_self <- res_df$Lead_Var == res_df$Lag_Var
    is_lag_zero <- res_df$Lag == 0

    # Dynamic threshold assignment (Disable for autocorrelation by setting to NA)
    res_df$Global_Threshold <- NA_real_
    res_df$Global_Threshold[!is_self & is_lag_zero] <- threshold_lag0 # Strict threshold for Lag 0 confounding
    res_df$Global_Threshold[!is_self & !is_lag_zero] <- threshold_lagged # Normal threshold for Lag > 0 propagation

    # Calculate Xi_Excess
    res_df$Xi_Excess <- NA_real_
    res_df$Xi_Excess[!is_self] <- pmax(
        0,
        res_df$Xi[!is_self] - res_df$Global_Threshold[!is_self]
    )

    # Store raw data for extractors
    df_raw <- as.data.frame(x_mat)
    colnames(df_raw) <- var_names

    out <- list(
        data = res_df,
        n = n,
        p = p,
        max_lag = max_lag,
        n_surr = n_surr,
        sig_level = sig_level,
        data_raw = df_raw
    )

    class(out) <- "xi_matrix"
    return(out)
}

#' Print method for xi_matrix
#' @param x An object of class \code{xi_matrix}.
#' @param ... Additional arguments passed to print.
#' @return The original object \code{x} invisibly.
#' @importFrom utils head
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

    sig_data <- x$data[
        which(x$data$Xi_Excess > 0),
        c("Lead_Var", "Lag_Var", "Lag", "Xi", "Global_Threshold", "Xi_Excess")
    ]

    if (nrow(sig_data) == 0) {
        cat(
            "No significant multivariate dependencies found above the global threshold.\n"
        )
    } else {
        cat("Top 5 Strongest Non-linear Causal Pathways:\n")
        print(
            utils::head(sig_data[order(-sig_data$Xi_Excess), ], 5),
            row.names = FALSE
        )
    }
    cat("\n")
    invisible(x)
}
