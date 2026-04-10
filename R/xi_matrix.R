#' Multivariate Xi-Correlogram Matrix
#'
#' Computes the pairwise directional Chatterjee's Xi coefficient for a multivariate
#' time series dataset. It evaluates both "Lead -> Lag" and "Lag -> Lead" relationships
#' across all variable pairs, as well as the Xi-ACF (autocorrelation) for individual variables.
#'
#' @param x A numeric matrix or data.frame containing the multivariate time series (columns = variables).
#' @param max_lag An integer specifying the maximum positive lag to compute.
#' @param n_surr An integer specifying the number of MIAAFT surrogate datasets for hypothesis testing.
#' @param sig_level A numeric value between 0 and 1 specifying the significance level for the surrogate threshold. Default is 0.95.
#'
#' @return An S3 object of class \code{xi_matrix} containing a tidy data frame of pairwise results.
#' @importFrom stats quantile na.omit
#' @export
xi_matrix <- function(x, max_lag = 20, n_surr = 100, sig_level = 0.95) {
    # --- 1. Input Validation ---
    if (!is.matrix(x) && !is.data.frame(x)) {
        stop("'x' must be a numeric matrix or data.frame.")
    }

    x_mat <- as.matrix(x)

    if (!is.numeric(x_mat)) {
        stop("All columns in 'x' must be numeric.")
    }
    if (any(is.na(x_mat))) {
        stop(
            "Input contains NA values. Please handle missing values before running."
        )
    }
    if (sig_level <= 0 || sig_level >= 1) {
        stop("'sig_level' must be strictly between 0 and 1.")
    }

    n <- nrow(x_mat)
    M <- ncol(x_mat)

    if (M < 2) {
        stop("'x' must contain at least two variables (columns).")
    }
    if (n <= max_lag) {
        stop("Time series length must be greater than 'max_lag'.")
    }

    # Handle column names
    var_names <- colnames(x_mat)
    if (is.null(var_names)) {
        var_names <- paste0("V", seq_len(M))
    }

    # --- 2. Call the C++ Matrix Engine ---
    res_cpp <- compute_xi_matrix_miaaft(x_mat, max_lag, n_surr)

    # --- 3. Process Surrogate Thresholds ---
    # Calculate the dynamically specified quantile for each row
    calc_threshold <- function(surr_matrix) {
        if (n_surr == 0) {
            return(rep(NA_real_, nrow(surr_matrix)))
        }
        apply(surr_matrix, 1, function(row) {
            stats::quantile(row, probs = sig_level, na.rm = TRUE)
        })
    }

    xi_threshold <- calc_threshold(res_cpp$xi_surrogates)

    # --- 4. Build Tidy DataFrame ---
    res_df <- data.frame(
        Lead_Var = var_names[res_cpp$var_lead],
        Lag_Var = var_names[res_cpp$var_lag],
        Lag = res_cpp$lag,
        Xi = res_cpp$xi_original,
        Xi_Threshold = xi_threshold # 汎用的な名前に変更
    )

    # Calculate excess Xi (clamped at 0 for non-significant values)
    res_df$Xi_Excess <- pmax(0, res_df$Xi - res_df$Xi_Threshold)

    # --- 5. Return S3 Object ---
    out <- list(
        data = res_df,
        n = n,
        M = M,
        max_lag = max_lag,
        n_surr = n_surr,
        sig_level = sig_level, # S3オブジェクトに保存
        var_names = var_names
    )
    class(out) <- "xi_matrix"

    return(out)
}

#' Print method for xi_matrix
#' @param x An object of class \code{xi_matrix}.
#' @param ... Additional arguments passed to print.
#' @return The original object \code{x} invisibly. Called primarily for its side effect of printing the matrix to the console.
#' @export
print.xi_matrix <- function(x, ...) {
    cat("\n=== Multivariate Xi-Correlogram Matrix ===\n")
    cat(sprintf(
        "Variables: %d (%s)\n",
        x$M,
        paste(x$var_names, collapse = ", ")
    ))
    cat(sprintf("Data length: %d\n", x$n))
    cat(sprintf("Max lag: %d\n", x$max_lag))
    cat(sprintf(
        "Surrogate threshold: %g%% (MIAAFT, n_surr = %d)\n\n",
        x$sig_level * 100,
        x$n_surr
    ))

    print(head(x$data))
    if (nrow(x$data) > 6) {
        cat(sprintf("... and %d more rows\n", nrow(x$data) - 6))
    }
    invisible(x)
}
