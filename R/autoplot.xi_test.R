#' Plot Xi-ACF Comparison
#'
#' Visualizes the comparison between standard linear ACF and non-linear Chatterjee's Xi
#' coefficient, including significance thresholds.
#'
#' @param object An object of class "xi_test".
#' @param ... Additional arguments (ignored).
#' @import ggplot2
#' @import latex2exp
#' @export
autoplot.xi_test <- function(object, ...) {
    # データを取り出し
    df <- object$data

    # ACFの信頼区間（定数）
    acf_ci <- object$data$ACF_CI[1]

    # ベースのプロット作成
    p <- ggplot(df, aes(x = Lag)) +
        # 0ライン
        geom_hline(yintercept = 0, color = "gray50", linewidth = 0.3) +

        # 1. ACF 95% 信頼区間 (青の点線)
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

    # 2. Xi サロゲート閾値 (グレーの帯)
    # n_surr > 0 の場合のみ描画 (NAチェック)
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

    # 3. メインのライン描画
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

        # --- Scales (色とラベルの定義) ---
        scale_color_manual(
            name = "Correlation Measure",
            # ★追加：アルファベット順の自動ソートを防ぎ、順番を固定する
            breaks = c(
                "Standard ACF (Linear)",
                "Chatterjee's Xi (Non-linear)"
            ),
            values = c(
                "Standard ACF (Linear)" = "steelblue",
                "Chatterjee's Xi (Non-linear)" = "firebrick"
            ),
            # ★追加：表示名だけTeXに変換
            labels = c(
                "Standard ACF (Linear)",
                TeX(r"($\xi$-ACF (Non-linear))")
            )
        ) +
        scale_fill_manual(
            name = "Significance",
            values = c("Xi 95% Threshold" = "gray50"),
            # ★追加：表示名だけTeXに変換
            labels = TeX(c(r"($\xi$-ACF 95% Threshold)"))
        ) +
        scale_linetype_manual(
            name = "Significance",
            values = c("ACF 95% CI" = "dotted")
        ) +

        # --- Theme & Labs ---
        labs(
            title = TeX(r"($\xi$-ACF Correlogram)"), # ★タイトルもXiを数式化
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

    # Y軸の範囲設定 (データに合わせて動的にズームイン)
    # グラフ内に描画される全要素の最大値と最小値を取得
    max_val <- max(c(df$ACF, df$Xi, df$Xi_Threshold_95, acf_ci), na.rm = TRUE)
    min_val <- min(c(df$ACF, df$Xi, -acf_ci, 0), na.rm = TRUE)

    # 上下に10%の余白 (マージン) を持たせる
    y_margin <- (max_val - min_val) * 0.1
    min_y <- min_val - y_margin
    max_y <- max_val + y_margin

    # 相関係数の理論上の上下限 (-1.0 ~ 1.0) を超えないようにキャップする
    max_y <- min(max_y, 1.05)
    min_y <- max(min_y, -1.05)

    p <- p + coord_cartesian(ylim = c(min_y, max_y))

    return(p)
}
