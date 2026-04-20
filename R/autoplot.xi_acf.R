#' Plot Xi-ACF Comparison
#'
#' Visualizes the comparison between the standard linear Autocorrelation Function (ACF)
#' and the non-linear Chatterjee's Xi coefficient, including their respective significance thresholds.
#'
#' @param object An object of class \code{"xi_acf"}.
#' @param ... Additional arguments passed to other methods.
#' @return A \code{ggplot} object representing the correlogram.
#' @importFrom ggplot2 autoplot ggplot aes geom_hline geom_ribbon geom_line geom_point guides guide_legend
#' @importFrom ggplot2 scale_color_manual scale_fill_manual scale_linetype_manual scale_x_continuous
#' @importFrom ggplot2 labs theme_minimal theme element_text coord_cartesian
#' @importFrom latex2exp TeX
#' @importFrom stats setNames
#' @method autoplot xi_acf
#' @export
autoplot.xi_acf <- function(object, ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }
    if (!requireNamespace("latex2exp", quietly = TRUE)) {
        stop("Package 'latex2exp' is required.")
    }

    df <- object$data
    acf_ci <- df$ACF_CI[1]
    sig_pct <- object$sig_level * 100

    acf_ci_label <- paste0("ACF ", sig_pct, "% CI")
    xi_thresh_label <- paste0("Xi ", sig_pct, "% Threshold")

    p <- ggplot(df, aes(x = Lag)) +
        geom_hline(yintercept = 0, color = "gray50", linewidth = 0.3) +
        geom_hline(
            aes(yintercept = acf_ci, linetype = acf_ci_label),
            color = "blue",
            linewidth = 0.4,
            alpha = 0.6,
            na.rm = TRUE
        ) +
        geom_hline(
            yintercept = -acf_ci,
            color = "blue",
            linetype = "dotted",
            linewidth = 0.4,
            alpha = 0.6,
            na.rm = TRUE
        )

    if (!all(is.na(df$Xi_Threshold))) {
        p <- p +
            geom_ribbon(
                aes(ymin = 0, ymax = Xi_Threshold, fill = xi_thresh_label),
                alpha = 0.2
            )
    }

    p <- p +
        geom_line(
            aes(y = ACF, color = "Standard ACF (Linear)"),
            linewidth = 0.6,
            linetype = "dashed",
            na.rm = TRUE
        ) +
        geom_point(
            aes(y = ACF, color = "Standard ACF (Linear)"),
            shape = 16,
            size = 3,
            na.rm = TRUE
        ) +
        geom_line(
            aes(y = Xi, color = "Chatterjee's Xi (Non-linear)"),
            linewidth = 0.8
        ) +

        # Highlight significant points (shape = 24)
        geom_point(
            aes(y = Xi, color = "Chatterjee's Xi (Non-linear)"),
            fill = ifelse(df$Xi_Excess > 0, "firebrick", "white"),
            shape = 24,
            size = ifelse(df$Xi_Excess > 0, 3, 1.5)
        ) +

        scale_color_manual(
            name = "Correlation Measure",
            breaks = c("Standard ACF (Linear)", "Chatterjee's Xi (Non-linear)"),
            values = c(
                "Standard ACF (Linear)" = "steelblue",
                "Chatterjee's Xi (Non-linear)" = "firebrick"
            ),
            labels = c(
                "Standard ACF (Linear)",
                TeX(r"($\xi$-ACF (Non-linear))")
            )
        ) +
        scale_fill_manual(
            name = "Significance",
            values = setNames("gray50", xi_thresh_label),
            labels = TeX(c(paste0(r"($\xi$-ACF )", sig_pct, r"(% Threshold)")))
        ) +
        scale_linetype_manual(
            name = "Significance",
            values = setNames("dotted", acf_ci_label)
        ) +
        # Override aesthetics to display legend icons correctly
        guides(
            color = guide_legend(
                override.aes = list(
                    shape = c(16, 24),
                    fill = c(NA, "firebrick"),
                    size = 3
                )
            )
        ) +
        labs(
            title = TeX(r"($\xi$-ACF Correlogram)"),
            subtitle = paste0(
                "Linear vs Non-linear Dependence (n = ",
                object$n,
                ")"
            ),
            y = "Correlation Coefficient",
            x = "Lag"
        ) +
        theme_minimal() +
        theme(
            legend.position = "bottom",
            legend.box = "vertical",
            legend.title = element_text(size = 9, face = "bold"),
            plot.title = element_text(face = "bold")
        )

    max_val <- max(c(df$ACF, df$Xi, df$Xi_Threshold, acf_ci), na.rm = TRUE)
    min_val <- min(c(df$ACF, df$Xi, -acf_ci, 0), na.rm = TRUE)
    y_margin <- (max_val - min_val) * 0.1
    p <- p +
        coord_cartesian(
            ylim = c(
                max(min_val - y_margin, -1.05),
                min(max_val + y_margin, 1.05)
            )
        ) +
        scale_x_continuous(breaks = function(x) {
            seq(
                ceiling(x[1]),
                floor(x[2]),
                by = max(1, floor((x[2] - x[1]) / 10))
            )
        })
    return(p)
}
