#' Plot Xi-ACF Comparison
#'
#' Visualizes the comparison between the standard linear Autocorrelation Function (ACF)
#' and the non-linear Chatterjee's Xi coefficient, including their respective significance thresholds.
#'
#' @param object An object of class \code{"xi_test"}.
#' @param ... Additional arguments passed to other methods (currently ignored).
#' @return A \code{ggplot} object representing the correlogram.
#' @importFrom ggplot2 autoplot ggplot aes geom_hline geom_ribbon geom_line geom_point
#' @importFrom ggplot2 scale_color_manual scale_fill_manual scale_linetype_manual
#' @importFrom ggplot2 labs theme_minimal theme element_text coord_cartesian
#' @importFrom latex2exp TeX
#' @export
autoplot.xi_test <- function(object, ...) {
    # Extract data from the object
    df <- object$data

    # ACF confidence interval (constant across lags)
    acf_ci <- object$data$ACF_CI[1]

    # Initialize the base plot
    p <- ggplot(df, aes(x = Lag)) +
        # Zero line
        geom_hline(yintercept = 0, color = "gray50", linewidth = 0.3) +

        # 1. Standard ACF 95% Confidence Interval (blue dotted line)
        geom_hline(
            aes(yintercept = acf_ci, linetype = "ACF 95% CI"),
            color = "blue",
            linewidth = 0.4,
            alpha = 0.6
        ) +
        geom_hline(
            yintercept = -acf_ci,
            color = "blue",
            linetype = "dotted",
            linewidth = 0.4,
            alpha = 0.6
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
        # ACF (Linear)
        geom_line(
            aes(y = ACF, color = "Standard ACF (Linear)"),
            linewidth = 0.6,
            linetype = "dashed"
        ) +
        geom_point(
            aes(y = ACF, color = "Standard ACF (Linear)"),
            shape = 16,
            size = 3
        ) +

        # Xi (Non-linear)
        geom_line(
            aes(y = Xi, color = "Chatterjee's Xi (Non-linear)"),
            linewidth = 0.8
        ) +
        geom_point(
            aes(y = Xi, color = "Chatterjee's Xi (Non-linear)"),
            shape = 17,
            size = 3
        ) +

        # --- Scales (Color and Label definitions) ---
        scale_color_manual(
            name = "Correlation Measure",
            # Prevent automatic alphabetical sorting to ensure consistent ordering
            breaks = c(
                "Standard ACF (Linear)",
                "Chatterjee's Xi (Non-linear)"
            ),
            values = c(
                "Standard ACF (Linear)" = "steelblue",
                "Chatterjee's Xi (Non-linear)" = "firebrick"
            ),
            # Convert display names to TeX format for mathematical symbols
            labels = c(
                "Standard ACF (Linear)",
                TeX(r"($\xi$-ACF (Non-linear))")
            )
        ) +
        scale_fill_manual(
            name = "Significance",
            values = c("Xi 95% Threshold" = "gray50"),
            labels = TeX(c(r"($\xi$-ACF 95% Threshold)"))
        ) +
        scale_linetype_manual(
            name = "Significance",
            values = c("ACF 95% CI" = "dotted")
        ) +

        # --- Theme & Labels ---
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

    # --- Dynamic Y-axis Zoom ---
    # Retrieve the maximum and minimum values among all plotted elements
    max_val <- max(c(df$ACF, df$Xi, df$Xi_Threshold_95, acf_ci), na.rm = TRUE)
    min_val <- min(c(df$ACF, df$Xi, -acf_ci, 0), na.rm = TRUE)

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
