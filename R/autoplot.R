#' @importFrom ggplot2 autoplot ggplot aes geom_hline geom_vline geom_ribbon geom_line geom_point scale_color_manual scale_fill_manual scale_size_manual scale_x_continuous labs theme_minimal theme element_text element_blank element_rect coord_cartesian guides guide_legend
#' @importFrom patchwork wrap_plots plot_annotation
#' @export
ggplot2::autoplot

#' Plot method for xi_acf objects
#'
#' @param object An object of class \code{xi_acf}.
#' @param ... Additional arguments passed to other methods.
#' @method autoplot xi_acf
#' @export
autoplot.xi_acf <- function(object, ...) {
    df <- object$data
    sig_pct <- object$sig_level * 100

    df$Significant <- df$Xi_Excess > 0

    p <- ggplot2::ggplot(df, ggplot2::aes(x = Lag)) +
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = 0, ymax = Global_Threshold),
            fill = "gray50",
            alpha = 0.2
        ) +
        ggplot2::geom_hline(yintercept = 0, color = "gray50", linewidth = 0.3) +
        ggplot2::geom_hline(
            yintercept = c(df$ACF_CI[1], -df$ACF_CI[1]),
            color = "blue",
            linetype = "dotted",
            linewidth = 0.4,
            alpha = 0.6
        ) +
        ggplot2::geom_line(
            ggplot2::aes(y = ACF, color = "Pearson (Linear)"),
            linetype = "dashed",
            linewidth = 0.6
        ) +
        ggplot2::geom_point(
            ggplot2::aes(y = ACF, color = "Pearson (Linear)"),
            shape = 16,
            size = 3
        ) +
        ggplot2::geom_line(
            ggplot2::aes(y = Xi, color = "Xi (Non-linear)"),
            linewidth = 0.8
        ) +
        ggplot2::geom_point(
            ggplot2::aes(
                y = Xi,
                color = "Xi (Non-linear)",
                fill = Significant,
                size = Significant
            ),
            shape = 24
        ) +
        ggplot2::scale_fill_manual(
            values = c("FALSE" = "white", "TRUE" = "firebrick"),
            guide = "none"
        ) +
        ggplot2::scale_size_manual(
            values = c("FALSE" = 1.5, "TRUE" = 3.0),
            guide = "none"
        ) +
        ggplot2::scale_color_manual(
            name = "Correlation Measure",
            values = c(
                "Xi (Non-linear)" = "firebrick",
                "Pearson (Linear)" = "steelblue"
            ),
            labels = c(
                "Xi (Non-linear)" = bquote(xi * " (Non-linear)"),
                "Pearson (Linear)" = "Pearson (Linear)"
            )
        ) +
        ggplot2::guides(
            color = ggplot2::guide_legend(
                override.aes = list(
                    shape = c(16, 24),
                    fill = c(NA, "firebrick"),
                    size = 3,
                    linetype = c("dashed", "solid")
                )
            )
        ) +
        ggplot2::scale_x_continuous(breaks = function(x) {
            seq(floor(x[1]), ceiling(x[2]), by = 1)
        }) +
        ggplot2::labs(
            title = expression(xi * "-Autocorrelation Function"),
            subtitle = bquote(
                "Gray Ribbon: FWER " * .(sig_pct) * "% Threshold"
            ),
            y = "Correlation Coefficient",
            x = "Lag"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "bottom")

    return(p)
}

#' Plot method for xi_ccf objects
#'
#' @param object An object of class \code{xi_ccf}.
#' @param ... Additional arguments passed to other methods.
#' @method autoplot xi_ccf
#' @export
autoplot.xi_ccf <- function(object, ...) {
    df <- object$data
    sig_pct <- object$sig_level * 100

    # Flag for significant pathways and create labels for faceting
    df$Significant <- df$Xi_Excess > 0
    df$Pathway <- paste(df$Lead_Var, " leads ", df$Lag_Var)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = Lag)) +
        # Max-Statistic Global Threshold (Shared Ribbon)
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = 0, ymax = Global_Threshold),
            fill = "gray50",
            alpha = 0.2
        ) +
        # Zero baseline
        ggplot2::geom_hline(yintercept = 0, color = "gray50", linewidth = 0.3) +
        # Linear CCF Confidence Intervals
        ggplot2::geom_hline(
            yintercept = c(df$CCF_CI[1], -df$CCF_CI[1]),
            color = "blue",
            linetype = "dotted",
            linewidth = 0.4,
            alpha = 0.6
        ) +
        # Pearson CCF (Linear) - Line
        ggplot2::geom_line(
            ggplot2::aes(y = CCF, color = "Pearson (Linear)"),
            linetype = "dashed",
            linewidth = 0.6
        ) +
        # Pearson CCF (Linear) - Points
        ggplot2::geom_point(
            ggplot2::aes(y = CCF, color = "Pearson (Linear)"),
            shape = 16,
            size = 3
        ) +
        # Xi CCF (Non-linear) - Line
        ggplot2::geom_line(
            ggplot2::aes(y = Xi, color = "Xi (Non-linear)"),
            linewidth = 0.8
        ) +
        # Xi CCF (Non-linear) - Points (Triangles)
        ggplot2::geom_point(
            ggplot2::aes(
                y = Xi,
                color = "Xi (Non-linear)",
                fill = Significant,
                size = Significant
            ),
            shape = 24
        ) +
        # Manual fill and size scales for significance highlighting
        ggplot2::scale_fill_manual(
            values = c("FALSE" = "white", "TRUE" = "firebrick"),
            guide = "none"
        ) +
        ggplot2::scale_size_manual(
            values = c("FALSE" = 1.5, "TRUE" = 3.0),
            guide = "none"
        ) +
        # Color mapping and custom legend labels
        ggplot2::scale_color_manual(
            name = "Correlation Measure",
            values = c(
                "Xi (Non-linear)" = "firebrick",
                "Pearson (Linear)" = "steelblue"
            ),
            labels = c(
                "Xi (Non-linear)" = bquote(xi * " (Non-linear)"),
                "Pearson (Linear)" = "Pearson (Linear)"
            )
        ) +
        # Custom legend aesthetic overrides
        ggplot2::guides(
            color = ggplot2::guide_legend(
                override.aes = list(
                    shape = c(16, 24),
                    fill = c(NA, "firebrick"),
                    size = 3,
                    linetype = c("dashed", "solid")
                )
            )
        ) +
        # Integer breaks for X-axis
        ggplot2::scale_x_continuous(breaks = function(x) {
            seq(floor(x[1]), ceiling(x[2]), by = 1)
        }) +
        # Facet by causal direction (stacked vertically)
        ggplot2::facet_wrap(~Pathway, ncol = 1) +
        # Labels and titles
        ggplot2::labs(
            title = expression(xi * "-Cross-Correlation Function"),
            subtitle = bquote(
                "Gray Ribbon: FWER " * .(sig_pct) * "% Threshold"
            ),
            y = "Correlation Coefficient",
            x = "Lag"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            legend.position = "bottom",
            strip.text = ggplot2::element_text(face = "bold", size = 11)
        )

    return(p)
}

#' Plot method for xi_matrix objects
#'
#' @param object An object of class \code{xi_matrix}.
#' @param ... Additional arguments passed to other methods.
#' @method autoplot xi_matrix
#' @export
autoplot.xi_matrix <- function(object, ...) {
    df <- object$data
    sig_pct <- object$sig_level * 100

    # Identify diagonal (auto-correlation) and off-diagonal (cross-correlation) elements
    df$is_diag <- df$Lead_Var == df$Lag_Var
    df$Significant <- df$Xi_Excess > 0

    # Separate data for plotting to completely avoid NA warnings
    df_off_diag <- df[!df$is_diag, ]

    # Add a dummy column for yintercept to strictly bind hline to off-diagonal panels
    df_off_diag$zero_line <- 0

    # Calculate exact mathematical center of the Lag axis for precise text placement
    center_lag <- (max(df$Lag, na.rm = TRUE) + min(df$Lag, na.rm = TRUE)) / 2

    # Create a minimal data frame for diagonal text labels
    vars <- unique(df$Lead_Var)
    df_diag <- data.frame(
        Lead_Var = vars,
        Lag_Var = vars,
        Lag = center_lag, # Placed at the exact horizontal center
        Xi = 0.5, # Placed at the exact vertical center
        Label = vars,
        stringsAsFactors = FALSE
    )

    p <- ggplot2::ggplot(mapping = ggplot2::aes(x = Lag, y = Xi)) +
        # Off-diagonal: Max-Statistic Global Threshold (Shared Ribbon)
        ggplot2::geom_ribbon(
            data = df_off_diag,
            ggplot2::aes(ymin = 0, ymax = Global_Threshold),
            fill = "gray50",
            alpha = 0.2
        ) +
        # Off-diagonal: Zero baseline (Mapped inside aes to strictly avoid panel bleeding and warnings)
        ggplot2::geom_hline(
            data = df_off_diag,
            ggplot2::aes(yintercept = zero_line),
            color = "gray50",
            linewidth = 0.3
        ) +
        # Off-diagonal: Directional dependency lines
        ggplot2::geom_line(
            data = df_off_diag,
            color = "firebrick",
            linewidth = 0.8
        ) +
        # Off-diagonal: Highlight significant points
        ggplot2::geom_point(
            data = df_off_diag,
            ggplot2::aes(fill = Significant, size = Significant),
            color = "firebrick",
            shape = 24
        ) +
        # Manual scales for significance highlighting
        ggplot2::scale_fill_manual(
            values = c("FALSE" = "white", "TRUE" = "firebrick"),
            guide = "none"
        ) +
        ggplot2::scale_size_manual(
            values = c("FALSE" = 1.5, "TRUE" = 3.0),
            guide = "none"
        ) +
        # Diagonal: Variable names as large text
        ggplot2::geom_text(
            data = df_diag,
            ggplot2::aes(label = Label),
            size = 7,
            fontface = "bold",
            color = "black",
            hjust = 0.5, # Force strict horizontal centering
            vjust = 0.5 # Force strict vertical centering
        ) +
        # Grid layout: Lead_Var (Rows) -> Lag_Var (Columns)
        ggplot2::facet_grid(Lead_Var ~ Lag_Var) +
        # Fixed Y-axis limits
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        # Integer breaks for X-axis
        ggplot2::scale_x_continuous(breaks = function(x) {
            seq(floor(x[1]), ceiling(x[2]), by = 1)
        }) +
        ggplot2::labs(
            title = "Multivariate Xi-Correlogram Matrix",
            subtitle = bquote(
                "Gray Ribbon: Max-Statistic FWER " * .(sig_pct) * "% Threshold"
            ),
            y = expression(xi ~ ~Coefficient),
            x = "Lag"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            strip.text = ggplot2::element_text(face = "bold", size = 11),
            # Add subtle borders to make the grid structure clear
            panel.border = ggplot2::element_rect(
                fill = NA,
                color = "gray80",
                linewidth = 0.5
            )
        )

    return(p)
}
