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
    acf_vals <- as.numeric(acf_res$acf)[-1]

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

#' @export
print.xi_test <- function(x, ...) {
    cat("Xi-ACF Test Results\n")
    cat("-------------------\n")
    cat(
        "Parameters: max_lag =",
        x$params$max_lag,
        ", n_surr =",
        x$params$n_surr,
        "\n\n"
    )
    print(head(x$summary, 10))
    if (nrow(x$summary) > 10) {
        cat("... (", nrow(x$summary) - 10, " more lags)\n")
    }
}

#' @import ggplot2
#' @export
autoplot.xi_test <- function(object, ...) {
    df <- object$summary

    ggplot(df, aes(x = Lag)) +
        # 1. Linear ACF (Gray Bars: 白黒印刷でも白飛びしないグレー)
        geom_bar(
            aes(y = ACF, fill = "Linear ACF"),
            stat = "identity",
            alpha = 0.5,
            width = 0.6
        ) +
        geom_hline(
            aes(yintercept = ACF_CI, color = "ACF 95% CI"),
            linetype = "dotted"
        ) +
        geom_hline(
            aes(yintercept = -ACF_CI, color = "ACF 95% CI"),
            linetype = "dotted"
        ) +
        # 2. Xi Threshold (Red Dashed) - 線の下に描画した方が綺麗
        {
            if (!all(is.na(df$Xi_Threshold_95))) {
                geom_line(
                    aes(y = Xi_Threshold_95, color = "Xi 95% Threshold"),
                    linetype = "dashed"
                )
            }
        } +
        # 3. Xi Coefficient (Dark Red Line & Points: 白黒印刷でも黒く出る)
        geom_line(aes(y = Xi, color = "Xi Coefficient"), linewidth = 1) +
        geom_point(aes(y = Xi, color = "Xi Coefficient"), size = 2) +

        # 配色の設定
        scale_fill_manual(
            name = "",
            values = c("Linear ACF" = "gray60") # 青からグレーに変更
        ) +
        scale_color_manual(
            name = "",
            values = c(
                "Xi Coefficient" = "firebrick",
                "Xi 95% Threshold" = "firebrick",
                "ACF 95% CI" = "gray50"
            )
        ) +
        # 前回決めた最強の軸ラベルをデフォルトに設定
        labs(
            x = "Lag",
            y = expression(hat(xi)[X](k))
        ) +
        theme_minimal() +
        theme(
            legend.position = "bottom",
            plot.title = element_text(face = "bold")
        )
}
