#' @importFrom ggplot2 autoplot ggplot aes geom_hline geom_vline geom_ribbon geom_line geom_point scale_color_manual scale_fill_manual scale_size_manual scale_x_continuous labs theme_minimal theme element_text element_blank element_rect coord_cartesian guides guide_legend
#' @importFrom patchwork wrap_plots plot_annotation
#' @export
ggplot2::autoplot

#' Plot method for xi_acf objects
#' @method autoplot xi_acf
#' @export
autoplot.xi_acf <- function(object, ...) {
    df <- object$data
    sig_pct <- object$sig_level * 100

    # 有意性のフラグを作成
    df$Significant <- df$Xi_Excess > 0

    p <- ggplot2::ggplot(df, ggplot2::aes(x = Lag)) +
        # FWER基準のグレーリボン (旧バージョンに合わせて薄めの gray50, alpha = 0.2)
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = 0, ymax = Global_Threshold),
            fill = "gray50",
            alpha = 0.2
        ) +
        ggplot2::geom_hline(yintercept = 0, color = "gray50", linewidth = 0.3) +
        # Pearsonの線形CI (blue, dotted, 0.4)
        ggplot2::geom_hline(
            yintercept = c(df$ACF_CI[1], -df$ACF_CI[1]),
            color = "blue",
            linetype = "dotted",
            linewidth = 0.4,
            alpha = 0.6
        ) +
        # Pearson (Linear): steelblue, dashed, linewidth 0.6, shape 16, size 3
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
        # Xi (Non-linear): firebrick, solid, linewidth 0.8
        ggplot2::geom_line(
            ggplot2::aes(y = Xi, color = "Xi (Non-linear)"),
            linewidth = 0.8
        ) +
        # Xiのポイント (有意なら大きく塗りつぶし、無意なら小さく白抜き)
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
            )
        ) +
        # 凡例のアイコンを正しく表示させるためのオーバーライド
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
            title = "Xi-Autocorrelation Function",
            subtitle = sprintf("Gray Ribbon: FWER %g%% Threshold", sig_pct),
            y = "Correlation Coefficient",
            x = "Lag"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "bottom")

    return(p)
}

#' Plot method for xi_ccf objects
#' @method autoplot xi_ccf
#' @export
autoplot.xi_ccf <- function(object, ...) {
    df <- object$data
    sig_pct <- object$sig_level * 100

    # 有意性のフラグを作成
    df$Significant <- df$Xi_Excess > 0

    base_plot <- function(data_sub, title_sub) {
        ggplot2::ggplot(data_sub, ggplot2::aes(x = Lag)) +
            ggplot2::geom_ribbon(
                ggplot2::aes(ymin = 0, ymax = Global_Threshold),
                fill = "gray50",
                alpha = 0.2
            ) +
            ggplot2::geom_hline(
                yintercept = 0,
                color = "gray50",
                linewidth = 0.3
            ) +
            ggplot2::geom_hline(
                yintercept = c(data_sub$CCF_CI[1], -data_sub$CCF_CI[1]),
                color = "blue",
                linetype = "dotted",
                linewidth = 0.4,
                alpha = 0.6
            ) +
            ggplot2::geom_line(
                ggplot2::aes(y = CCF, color = "Pearson (Linear)"),
                linetype = "dashed",
                linewidth = 0.6
            ) +
            ggplot2::geom_point(
                ggplot2::aes(y = CCF, color = "Pearson (Linear)"),
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
                )
            ) +
            ggplot2::scale_x_continuous(breaks = function(x) {
                seq(floor(x[1]), ceiling(x[2]), by = 1)
            }) +
            ggplot2::theme_minimal() +
            ggplot2::labs(
                subtitle = title_sub,
                y = "Correlation Coefficient",
                x = "Lag"
            ) +
            ggplot2::coord_cartesian(
                ylim = c(min(data_sub$CCF, -0.2, na.rm = TRUE), 1.05)
            )
    }

    df_fwd <- df[df$Lag >= 0, ]
    df_bwd <- df[df$Lag <= 0, ]
    df_bwd$Lag <- abs(df_bwd$Lag)
    df_bwd <- df_bwd[order(df_bwd$Lag), ]

    p1 <- base_plot(df_fwd, "Direction: X leads Y")
    p2 <- base_plot(df_bwd, "Direction: Y leads X")

    combined_plot <- patchwork::wrap_plots(
        p1,
        p2,
        ncol = 1,
        guides = "collect"
    ) &
        ggplot2::theme(
            legend.position = "bottom",
            legend.box = "vertical"
        ) &
        ggplot2::guides(
            color = ggplot2::guide_legend(
                override.aes = list(
                    shape = c(16, 24),
                    fill = c(NA, "firebrick"),
                    size = 3,
                    linetype = c("dashed", "solid")
                )
            )
        )

    return(
        combined_plot +
            patchwork::plot_annotation(
                title = "Xi-Cross-Correlation Function",
                subtitle = sprintf("MIAAFT FWER %g%% Control", sig_pct)
            )
    )
}

#' Plot method for xi_matrix objects
#' @method autoplot xi_matrix
#' @export
autoplot.xi_matrix <- function(object, ...) {
    df <- object$data
    sig_pct <- object$sig_level * 100

    # 有意性のフラグを作成
    df$Significant <- df$Xi_Excess > 0

    p <- ggplot2::ggplot(df, ggplot2::aes(x = Lag, y = Xi)) +
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = 0, ymax = Global_Threshold),
            fill = "gray50",
            alpha = 0.2
        ) +
        ggplot2::geom_hline(yintercept = 0, color = "gray50", linewidth = 0.3) +
        ggplot2::geom_line(color = "firebrick", linewidth = 0.8) +
        ggplot2::geom_point(
            ggplot2::aes(fill = Significant, size = Significant),
            color = "firebrick",
            shape = 24
        ) +
        ggplot2::scale_fill_manual(
            values = c("FALSE" = "white", "TRUE" = "firebrick")
        ) +
        ggplot2::scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 3.0)) +
        ggplot2::facet_grid(Lead_Var ~ Lag_Var) +
        ggplot2::scale_x_continuous(breaks = function(x) {
            seq(floor(x[1]), ceiling(x[2]), by = 1)
        }) +
        ggplot2::labs(
            title = "Multivariate Xi-Correlogram Matrix",
            subtitle = sprintf(
                "Rows: Lead (Predictor)  |  Columns: Lag (Response)\nGray ribbon: %g%% MIAAFT Threshold",
                sig_pct
            ),
            y = "Chatterjee's Xi",
            x = "Lag"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            strip.background = ggplot2::element_rect(
                fill = "gray40",
                color = NA
            ),
            strip.text = ggplot2::element_text(
                face = "bold",
                color = "white",
                size = 11
            ),
            axis.text = ggplot2::element_text(size = 7),
            panel.border = ggplot2::element_rect(color = "gray80", fill = NA),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position = "none"
        )

    return(p)
}
