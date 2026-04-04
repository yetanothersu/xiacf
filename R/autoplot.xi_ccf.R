#' Plot Directional Xi-CCF Comparison
#'
#' Visualizes the comparison between the standard linear Cross-Correlation Function (CCF)
#' and the non-linear Chatterjee's Xi cross-correlation across two directions (X leads Y, Y leads X).
#'
#' @param object An object of class \code{"xi_ccf"}.
#' @param ... Additional arguments passed to other methods.
#' @return A \code{ggplot} object representing the directional cross-correlogram.
#' @importFrom ggplot2 autoplot ggplot aes geom_hline geom_ribbon geom_line geom_point guides guide_legend
#' @importFrom ggplot2 scale_color_manual scale_fill_manual scale_linetype_manual coord_cartesian facet_wrap
#' @importFrom ggplot2 labs theme_minimal theme element_text element_rect
#' @importFrom latex2exp TeX
#' @importFrom stats setNames
#' @method autoplot xi_ccf
#' @export
autoplot.xi_ccf <- function(object, ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }
    if (!requireNamespace("latex2exp", quietly = TRUE)) {
        stop("Package 'latex2exp' is required.")
    }

    df <- object$data
    ccf_ci <- df$CCF_CI[1]
    sig_pct <- object$sig_level * 100

    ccf_ci_label <- paste0("CCF ", sig_pct, "% CI")
    xi_thresh_label <- paste0("Xi ", sig_pct, "% Threshold")

    p <- ggplot(df, aes(x = Lag)) +
        geom_hline(yintercept = 0, color = "gray50", linewidth = 0.3) +
        geom_hline(
            data = data.frame(y_vals = c(ccf_ci, -ccf_ci)),
            aes(yintercept = y_vals, linetype = ccf_ci_label),
            color = "blue",
            linewidth = 0.4,
            alpha = 0.6
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
            aes(y = CCF, color = "Standard CCF (Linear)"),
            linewidth = 0.6,
            linetype = "dashed"
        ) +
        geom_point(
            aes(y = CCF, color = "Standard CCF (Linear)"),
            shape = 16,
            size = 3
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
        facet_wrap(~Direction, ncol = 1) +

        scale_fill_manual(
            name = "Significance",
            values = setNames("gray50", xi_thresh_label),
            labels = TeX(c(paste0(r"($\xi$-CCF )", sig_pct, r"(% Threshold)")))
        ) +
        scale_linetype_manual(
            name = "Significance",
            values = setNames("dotted", ccf_ci_label)
        ) +
        scale_color_manual(
            name = "Correlation Measure",
            breaks = c("Standard CCF (Linear)", "Chatterjee's Xi (Non-linear)"),
            values = c(
                "Standard CCF (Linear)" = "steelblue",
                "Chatterjee's Xi (Non-linear)" = "firebrick"
            ),
            labels = c(
                "Standard CCF (Linear)",
                TeX(r"($\xi$-CCF (Non-linear))")
            )
        ) +
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
            title = latex2exp::TeX(r"($\xi$-CCF Correlogram)"),
            subtitle = paste0(
                "Linear vs Non-linear Cross-Dependence (n = ",
                object$n,
                ", surrogates = ",
                object$n_surr,
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
            plot.title = element_text(face = "bold"),
            strip.text = element_text(
                size = 11,
                face = "bold",
                color = "white"
            ),
            strip.background = element_rect(fill = "gray40", color = NA)
        )

    max_val <- max(c(df$CCF, df$Xi, df$Xi_Threshold, ccf_ci), na.rm = TRUE)
    min_val <- min(c(df$CCF, df$Xi, -ccf_ci, 0), na.rm = TRUE)
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
