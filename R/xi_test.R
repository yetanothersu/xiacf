#' Xi-ACF Test for Time Series
#'
#' Calculates Chatterjee's Xi and Standard ACF with significance thresholds.
#'
#' @param x Numeric vector (time series).
#' @param max_lag Integer. Maximum lag to compute.
#' @param n_surr Integer. Number of surrogates for IAAFT test.
#' @return An object of class "xi_test".
#' @export
xi_test <- function(x, max_lag = 20, n_surr = 100) {
    # 1. Validate input
    if (!is.numeric(x) || length(x) < 5) {
        stop("x must be a numeric vector with length >= 5")
    }

    # 2. Calculate Standard ACF (Linear)
    # plot=FALSE で数値だけ取得
    acf_res <- stats::acf(x, lag.max = max_lag, plot = FALSE)
    acf_vals <- as.numeric(acf_res$acf)[-1] # lag 0 (常に1) を除く

    # Calculate ACF Confidence Interval (95% i.i.d.)
    n <- length(x)
    acf_ci <- qnorm(0.975) / sqrt(n)

    # 3. Calculate Xi (Non-linear) using C++ backend
    # パッケージ内のC++関数を呼び出す
    xi_res <- run_xi_test_cpp(x, max_lag, n_surr)

    # Calculate Xi Threshold from Surrogates
    xi_threshold <- apply(xi_res$xi_surrogates, 1, function(row) {
        stats::quantile(row, 0.95, na.rm = TRUE)
    })

    # 4. Construct S3 Object
    # データフレームにまとめておくとggplotで使いやすい
    df_summary <- data.frame(
        Lag = 1:max_lag,
        ACF = acf_vals,
        Xi = as.numeric(xi_res$xi_original),
        Xi_Threshold_95 = xi_threshold,
        ACF_CI = acf_ci
    )

    structure(
        list(
            data = df_summary,
            n = n,
            max_lag = max_lag,
            n_surr = n_surr,
            call = match.call()
        ),
        class = "xi_test"
    )
}

#' @export
print.xi_test <- function(x, ...) {
    cat("\n=== Xi-ACF Test Results ===\n")
    cat("Sample size:", x$n, "\n")
    cat("Max Lag:    ", x$max_lag, "\n")
    cat("Surrogates: ", x$n_surr, "\n\n")
    print(head(x$data, 5))
    cat("... (showing first 5 lags)\n")
}
