#' Plot Multivariate Xi-Correlogram Matrix
#'
#' Visualizes the result of a multivariate Xi-matrix analysis using a facet grid.
#' Rows represent the leading (predictor) variable, and columns represent the lagging (response) variable.
#'
#' @param object An object of class \code{xi_matrix}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A \code{ggplot} object.
#' @import ggplot2
#' @export
autoplot.xi_matrix <- function(object, ...) {
    df <- object$data
    sig_pct <- object$sig_level * 100

    # 有意性のフラグを作成（Xi_Excess が 0 より大きいか）
    df$Significant <- df$Xi_Excess > 0

    p <- ggplot2::ggplot(df, ggplot2::aes(x = Lag, y = Xi)) +
        # Zero line
        ggplot2::geom_hline(yintercept = 0, color = "gray50", linewidth = 0.3) +

        # 1. Xi Surrogate Threshold (Gray ribbon)
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = 0, ymax = Xi_Threshold),
            fill = "gray50",
            alpha = 0.2
        ) +

        # 2. Chatterjee's Xi (Line)
        ggplot2::geom_line(
            color = "firebrick",
            linewidth = 0.8
        ) +

        # 3. Points (Highlight significant ones)
        # shape = 24 は枠線(color)と中身(fill)の色を個別に設定できる三角です
        ggplot2::geom_point(
            ggplot2::aes(fill = Significant, size = Significant),
            color = "firebrick",
            shape = 24
        ) +

        # 有意かどうかに応じて、塗りつぶしの色とサイズを動的に変更
        ggplot2::scale_fill_manual(
            values = c("FALSE" = "white", "TRUE" = "firebrick")
        ) +
        ggplot2::scale_size_manual(
            values = c("FALSE" = 1.5, "TRUE" = 3.0) # 有意な点は2倍の大きさに！
        ) +

        # Create the N x N Matrix Grid
        ggplot2::facet_grid(Lead_Var ~ Lag_Var) +

        # Labels and Theme
        ggplot2::labs(
            title = "Multivariate Xi-Correlogram Matrix",
            subtitle = sprintf(
                "Rows: Lead (Predictor)  |  Columns: Lag (Response)\nGray ribbon: %g%% MIAAFT Threshold",
                sig_pct
            ),
            x = "Lag",
            y = "Chatterjee's Xi"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            strip.background = ggplot2::element_rect(
                fill = "gray40",
                color = NA
            ),
            strip.text = ggplot2::element_text(
                face = "bold",
                size = 11,
                color = "white"
            ),
            panel.border = ggplot2::element_rect(color = "gray80", fill = NA),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position = "none" # 凡例は不要なので消す
        )

    # --- Dynamic Y-axis Zoom ---
    max_val <- max(c(df$Xi, df$Xi_Threshold), na.rm = TRUE)
    min_val <- min(c(df$Xi, 0), na.rm = TRUE)

    y_margin <- (max_val - min_val) * 0.1
    min_y <- min_val - y_margin
    max_y <- max_val + y_margin

    max_y <- min(max_y, 1.05)
    min_y <- max(min_y, -0.1)

    p <- p +
        ggplot2::coord_cartesian(ylim = c(min_y, max_y)) +
        ggplot2::scale_x_continuous(breaks = function(x) {
            seq(
                ceiling(x[1]),
                floor(x[2]),
                by = max(1, floor((x[2] - x[1]) / 10))
            )
        })

    return(p)
}
