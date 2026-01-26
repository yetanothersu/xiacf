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
    acf_res <- stats::acf(x, lag.max = max_lag, plot = FALSE)
    acf_vals <- as.numeric(acf_res$acf)[-1] # lag 0を除外

    # Calculate ACF Confidence Interval (95% i.i.d.)
    n <- length(x)
    acf_ci <- qnorm(0.975) / sqrt(n)

    # 3. Calculate Xi (Non-linear) using C++ backend
    xi_res <- run_xi_test_cpp(x, max_lag, n_surr)

    # Calculate Xi Threshold from Surrogates (Safe for n_surr = 0)
    if (n_surr > 0) {
        xi_threshold <- apply(xi_res$xi_surrogates, 1, function(row) {
            stats::quantile(row, 0.95, na.rm = TRUE)
        })
    } else {
        xi_threshold <- rep(NA_real_, max_lag)
    }

    # 4. Construct S3 Object
    # 注: カラム名は 'Xi' です (Xi_Original ではありません)
    df_summary <- data.frame(
        Lag = 1:max_lag,
        ACF = acf_vals,
        Xi = as.numeric(xi_res$xi_original),
        Xi_Threshold_95 = xi_threshold,
        ACF_CI = rep(acf_ci, max_lag)
    )

    structure(
        list(
            summary = df_summary,
            params = list(max_lag = max_lag, n_surr = n_surr),
            data = x
        ),
        class = "xi_test"
    )
}

# printメソッドなども必要ならここに定義
