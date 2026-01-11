#' Plot Xi-ACF Comparison
#'
#' @param object An object of class "xi_test".
#' @param ... Additional arguments (ignored).
#' @import ggplot2
#' @export
autoplot.xi_test <- function(object, ...) {
    df <- object$data
    acf_ci <- object$data$ACF_CI[1] # 定数として取得

    ggplot(df, aes(x = Lag)) +
        # ゼロライン
        geom_hline(yintercept = 0, color = "gray50", linewidth = 0.3) +

        # 1. ACF信頼区間 (青点線)
        geom_hline(
            yintercept = c(-acf_ci, acf_ci),
            color = "blue",
            linetype = "dotted",
            linewidth = 0.5,
            alpha = 0.6
        ) +

        # 2. Xiサロゲート閾値 (グレー帯)
        geom_ribbon(
            aes(
                ymin = 0,
                ymax = Xi_Threshold_95,
                fill = "95% Null Threshold (IAAFT)"
            ),
            alpha = 0.3
        ) +

        # 3. ACF Line & Point
        geom_line(
            aes(
                y = ACF,
                color = "Standard ACF (Linear)",
                linetype = "Standard ACF (Linear)"
            ),
            linewidth = 0.8
        ) +
        geom_point(
            aes(
                y = ACF,
                color = "Standard ACF (Linear)",
                shape = "Standard ACF (Linear)"
            ),
            size = 2.5
        ) +

        # 4. Xi Line & Point
        geom_line(
            aes(
                y = Xi,
                color = "Chatterjee's Xi (Non-linear)",
                linetype = "Chatterjee's Xi (Non-linear)"
            ),
            linewidth = 0.8
        ) +
        geom_point(
            aes(
                y = Xi,
                color = "Chatterjee's Xi (Non-linear)",
                shape = "Chatterjee's Xi (Non-linear)"
            ),
            size = 3.0
        ) +

        # --- Scales ---
        scale_color_manual(
            name = NULL,
            values = c(
                "Standard ACF (Linear)" = "blue",
                "Chatterjee's Xi (Non-linear)" = "red"
            )
        ) +
        scale_linetype_manual(
            name = NULL,
            values = c(
                "Standard ACF (Linear)" = "dashed",
                "Chatterjee's Xi (Non-linear)" = "solid"
            )
        ) +
        scale_shape_manual(
            name = NULL,
            values = c(
                "Standard ACF (Linear)" = 16,
                "Chatterjee's Xi (Non-linear)" = 17
            )
        ) +
        scale_fill_manual(
            name = NULL,
            values = c("95% Null Threshold (IAAFT)" = "gray80")
        ) +

        # --- Design ---
        labs(
            title = "Xi-ACF Correlogram",
            subtitle = paste0(
                "Comparison of Linear (ACF) vs Non-linear (Xi) Dependence (N=",
                object$n,
                ")"
            ),
            y = "Correlation Coefficient",
            x = "Lag"
        ) +
        scale_y_continuous(limits = c(-0.2, 1.05)) +
        theme_minimal() +
        theme(
            legend.position = c(0.8, 0.8),
            legend.background = element_rect(fill = "white", color = "gray90"),
            legend.box = "vertical",
            plot.title = element_text(face = "bold"),
            axis.title = element_text(face = "bold")
        )
}
