#' Directional Xi-CCF Test for Multivariate Time Series
#'
#' Calculates Chatterjee's Xi cross-correlation and the standard Cross-Correlation Function (CCF)
#' across positive lags to evaluate directional lead-lag relationships.
#'
#' @useDynLib xiacf, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats na.omit sd ccf quantile qnorm
#'
#' @param x A numeric vector representing the first time series.
#' @param y A numeric vector representing the second time series.
#' @param max_lag An integer specifying the maximum positive lag.
#' @param n_surr An integer specifying the number of MIAAFT surrogate datasets.
#' @param bidirectional Logical. If TRUE (default), computes both directions.
#' @param sig_level A numeric value between 0 and 1 specifying the significance level. Default is 0.95.
#' @return An object of class \code{"xi_ccf"}.
#' @rdname xi_ccf
#' @export
xi_ccf <- function(
    x,
    y,
    max_lag = 20,
    n_surr = 100,
    bidirectional = TRUE,
    sig_level = 0.95
) {
    if (!is.numeric(x) || !is.numeric(y)) {
        stop("Inputs 'x' and 'y' must be numeric vectors.")
    }
    if (length(x) != length(y)) {
        stop("Time series 'x' and 'y' must have the exact same length.")
    }
    if (any(is.na(x)) || any(is.na(y))) {
        stop("Inputs contain NA values.")
    }
    if (stats::var(x) == 0 || stats::var(y) == 0) {
        stop("Time series has zero variance.")
    }
    if (sig_level <= 0 || sig_level >= 1) {
        stop("'sig_level' must be strictly between 0 and 1.")
    }

    n <- length(x)
    xi_res <- compute_xi_ccf_miaaft(x, y, max_lag, n_surr)

    calc_threshold <- function(surr_matrix) {
        apply(surr_matrix, 1, function(row) {
            stats::quantile(row, probs = sig_level, na.rm = TRUE)
        })
    }

    # Dynamically link CCF confidence interval to sig_level
    ccf_ci <- stats::qnorm((1 + sig_level) / 2) / sqrt(n)

    ccf_fwd_full <- stats::ccf(
        y,
        x,
        lag.max = max_lag,
        plot = FALSE,
        na.action = stats::na.pass
    )
    df_fwd <- data.frame(
        Direction = "X leads Y",
        Lag = as.numeric(xi_res$lags),
        CCF = as.numeric(ccf_fwd_full$acf[(max_lag + 1):(2 * max_lag + 1)]),
        Xi = as.numeric(xi_res$xi_original_forward),
        Xi_Threshold = calc_threshold(xi_res$xi_surrogates_forward),
        CCF_CI = ccf_ci
    )
    df_fwd$Xi_Excess <- pmax(0, df_fwd$Xi - df_fwd$Xi_Threshold)

    if (bidirectional) {
        ccf_bwd_full <- stats::ccf(
            x,
            y,
            lag.max = max_lag,
            plot = FALSE,
            na.action = stats::na.pass
        )
        df_bwd <- data.frame(
            Direction = "Y leads X",
            Lag = as.numeric(xi_res$lags),
            CCF = as.numeric(ccf_bwd_full$acf[(max_lag + 1):(2 * max_lag + 1)]),
            Xi = as.numeric(xi_res$xi_original_backward),
            Xi_Threshold = calc_threshold(xi_res$xi_surrogates_backward),
            CCF_CI = ccf_ci
        )
        df_bwd$Xi_Excess <- pmax(0, df_bwd$Xi - df_bwd$Xi_Threshold)
        df_res <- rbind(df_fwd, df_bwd)
    } else {
        df_res <- df_fwd
    }

    structure(
        list(
            data = df_res,
            n = n,
            max_lag = max_lag,
            n_surr = n_surr,
            bidirectional = bidirectional,
            sig_level = sig_level
        ),
        class = "xi_ccf"
    )
}

#' Print method for xi_ccf objects
#'
#' @param x An object of class \code{"xi_ccf"}.
#' @param ... Additional arguments.
#' @return Invisibly returns the original object.
#' @importFrom utils head
#' @export
print.xi_ccf <- function(x, ...) {
    cat("\n=== Directional Chatterjee's Xi Cross-Correlation ===\n")
    cat("Time series length (n):", x$n, "\n")
    cat("Maximum positive lag:  ", x$max_lag, "\n")
    cat(sprintf(
        "Significance level:    %g%% (MIAAFT, n_surr = %d)\n",
        x$sig_level * 100,
        x$n_surr
    ))
    cat("Bidirectional evaluation:", x$bidirectional, "\n\n")
    print(utils::head(x$data, 10), row.names = FALSE)
    if (nrow(x$data) > 10) {
        cat("... (Showing first 10 rows)\n\n")
    }
    invisible(x)
}
