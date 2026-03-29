#' Plot Xi-CCF Comparison
#'
#' Visualizes the comparison between the standard linear Cross-Correlation Function (CCF)
#' and the non-linear Chatterjee's Xi cross-correlation, including their significance thresholds.
#'
#' @param object An object of class \code{"xi_ccf"}.
#' @param ... Additional arguments passed to other methods.
#' @return A \code{ggplot} object representing the cross-correlogram.
#' @importFrom ggplot2 autoplot ggplot aes geom_hline geom_vline geom_ribbon geom_line geom_point
#' @importFrom ggplot2 scale_color_manual scale_fill_manual scale_linetype_manual scale_x_continuous coord_cartesian
#' @importFrom ggplot2 labs theme_minimal theme element_text
#' @importFrom latex2exp TeX
#' @method autoplot xi_ccf
#' @export
autoplot.xi_ccf <- function(object, ...) {
    # Check for required packages
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting. Please install it.")
    }
    if (!requireNamespace("latex2exp", quietly = TRUE)) {
        stop(
            "Package 'latex2exp' is required for LaTeX rendering in plots. Please install it."
        )
    }
    # Extract data from the object
    df <- object$data
    ccf_ci <- df$CCF_CI[1]

    # Initialize the base plot
    p <- ggplot(df, aes(x = Lag)) +
        # Zero lines
        geom_hline(yintercept = 0, color = "gray50", linewidth = 0.3) +
        geom_vline(
            xintercept = 0,
            color = "gray50",
            linewidth = 0.3,
            linetype = "dashed"
        ) +

        # 1. Standard CCF 95% Confidence Interval (blue dotted line)
        geom_hline(
            aes(yintercept = ccf_ci, linetype = "CCF 95% CI"),
            color = "blue",
            linewidth = 0.4,
            alpha = 0.6
        ) +
        geom_hline(
            yintercept = -ccf_ci,
            color = "blue",
            linewidth = 0.4,
            alpha = 0.6,
            linetype = "dotted"
        )

    # 2. Xi Surrogate Threshold (gray ribbon)
    # Draw only if surrogate data exists (i.e., n_surr > 0, no NA in threshold)
    if (!all(is.na(df$Xi_Threshold_95))) {
        p <- p +
            geom_ribbon(
                aes(
                    ymin = 0,
                    ymax = Xi_Threshold_95,
                    fill = "Xi 95% Threshold"
                ),
                alpha = 0.2
            )
    }

    # 3. Main lines and points
    p <- p +
        # CCF (Linear)
        geom_line(
            aes(y = CCF, color = "Standard CCF (Linear)"),
            linewidth = 0.6,
            linetype = "dashed"
        ) +
        geom_point(
            aes(y = CCF, color = "Standard CCF (Linear)"),
            size = 3,
            shape = 16
        ) +

        # Xi (Non-linear)
        geom_line(
            aes(y = Xi, color = "Chatterjee's Xi (Non-linear)"),
            linewidth = 0.8
        ) +
        geom_point(
            aes(y = Xi, color = "Chatterjee's Xi (Non-linear)"),
            size = 3,
            shape = 17
        ) +

        # --- Scales (Color and Label definitions) ---
        scale_color_manual(
            name = "Correlation Measure",
            # Prevent automatic alphabetical sorting to ensure consistent ordering
            breaks = c(
                "Standard CCF (Linear)",
                "Chatterjee's Xi (Non-linear)"
            ),
            values = c(
                "Standard CCF (Linear)" = "steelblue",
                "Chatterjee's Xi (Non-linear)" = "firebrick"
            ),
            # Convert display names to TeX format for mathematical symbols
            labels = c(
                "Standard CCF (Linear)",
                TeX(r"($\xi$-CCF (Non-linear))")
            )
        ) +
        scale_fill_manual(
            name = "Significance",
            values = c("Xi 95% Threshold" = "gray50"),
            labels = TeX(c(r"($\xi$-CCF 95% Threshold)"))
        ) +
        scale_linetype_manual(
            name = "Significance",
            values = c("CCF 95% CI" = "dotted")
        ) +

        # --- Theme & Labels ---
        labs(
            title = TeX(r"(Multivariate $\xi$-CCF Correlogram)"),
            subtitle = paste0(
                "Linear vs Non-linear Cross-Dependence (n = ",
                object$n,
                ")"
            ),
            y = "Correlation Coefficient",
            x = "Lag (Positive lag means X leads Y)"
        ) +
        theme_minimal() +
        theme(
            legend.position = "bottom",
            legend.box = "vertical",
            legend.title = element_text(size = 9, face = "bold"),
            plot.title = element_text(face = "bold")
        )

    # --- Dynamic Y-axis Zoom ---
    # Retrieve the maximum and minimum values among all plotted elements
    max_val <- max(c(df$CCF, df$Xi, df$Xi_Threshold_95, ccf_ci), na.rm = TRUE)
    min_val <- min(c(df$CCF, df$Xi, -ccf_ci, 0), na.rm = TRUE)

    # Add a 10% margin to the top and bottom
    y_margin <- (max_val - min_val) * 0.1
    min_y <- min_val - y_margin
    max_y <- max_val + y_margin

    # Cap the limits so they do not exceed the theoretical bounds of correlation (-1.0 to 1.0)
    # Using 1.05 to give a slight visual buffer even at maximum correlation
    max_y <- min(max_y, 1.05)
    min_y <- max(min_y, -1.05)

    p <- p + coord_cartesian(ylim = c(min_y, max_y))

    return(p)
}
