#' Xi-ACF Test for Time Series
#'
#' Calculates Chatterjee's Xi and Standard ACF with significance thresholds.
#'
#' @useDynLib xiacf, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats na.omit sd acf quantile
#'
#' @param x Numeric vector (time series).
#' @param max_lag Integer. Maximum lag to compute.
#' @param n_surr Integer. Number of surrogates for IAAFT test.
#' @return An object of class "xi_test".
#' @export
xi_test <- function(x, max_lag = 20, n_surr = 100) {
    # --- 1. Robust Input Validation (監査修正ポイント) ---

    # 数値ベクトルか確認
    if (!is.numeric(x)) {
        stop("Input 'x' must be a numeric vector.")
    }

    # NAチェック: 削除して警告を出す、またはエラーにする
    if (any(is.na(x))) {
        warning("'x' contains NA values. Removing them before analysis.")
        x <- stats::na.omit(x)
        # na.omit後の属性(attributes)を消して純粋なベクトルにする
        x <- as.numeric(x)
    }

    n <- length(x)

    # データ長チェック
    if (n < 5) {
        stop("Time series length is too short (n < 5).")
    }

    # 定数ベクトルチェック (分散ゼロ回避)
    if (stats::sd(x) == 0) {
        stop(
            "Input 'x' is constant (zero variance). Correlation cannot be computed."
        )
    }

    # ラグの長さチェック
    if (max_lag >= n) {
        warning("max_lag is >= series length. Reducing max_lag to n - 1.")
        max_lag <- n - 1
    }

    # --- 2. Calculate Standard ACF (Linear) ---
    # stats::acf は plot=FALSE でも計算してくれる
    # na.action = na.pass は不要(上で除去済み)
    acf_res <- stats::acf(x, lag.max = max_lag, plot = FALSE)

    # lag=0 (必ず1.0) を除外
    acf_vals <- as.numeric(acf_res$acf)[-1]

    # Calculate ACF Confidence Interval (95% i.i.d. assumption: +/- 1.96/sqrt(n))
    acf_ci <- stats::qnorm(0.975) / sqrt(n)

    # --- 3. Calculate Xi (Non-linear) using C++ backend ---
    # C++側で初期化修正済みなので安全
    xi_res <- compute_xi_lags(x, max_lag, n_surr)

    # --- 4. Calculate Threshold from Surrogates ---
    # n_surr = 0 の場合のハンドリング
    xi_threshold <- rep(NA, max_lag)

    if (n_surr > 0) {
        # 行ごと(各ラグごと)に95%点検定値を計算
        # C++側で計算不能なラグは NaN で埋められているため、na.rm = TRUE が必須
        xi_threshold <- apply(xi_res$xi_surrogates, 1, function(row) {
            stats::quantile(row, 0.95, na.rm = TRUE)
        })
    }

    # --- 5. Construct Result Object ---
    # データフレームにまとめることで、ggplot2での描画が楽になる
    df_res <- data.frame(
        Lag = 1:max_lag,
        ACF = acf_vals[1:max_lag], # acfの結果がmax_lagより短い場合に備えて安全策
        Xi = as.numeric(xi_res$xi_original),
        Xi_Threshold_95 = xi_threshold,
        ACF_CI = acf_ci # 定数だがプロット用に列として持たせても良いし、属性でも良い
    )

    # クラス定義
    structure(
        list(
            data = df_res,
            n = n,
            max_lag = max_lag,
            n_surr = n_surr
        ),
        class = "xi_test"
    )
}

#' Print method for xi_test objects
#'
#' @param x An object of class "xi_test".
#' @param ... Additional arguments (ignored).
#' @export
print.xi_test <- function(x, ...) {
    cat("\n\tChatterjee's Xi-ACF Test\n\n")
    cat("Data length:  ", x$n, "\n")
    cat("Max lag:      ", x$max_lag, "\n")
    cat("Surrogates:   ", x$n_surr, " (IAAFT)\n\n")

    # 表示したい列だけをピックアップ
    # ACF_CI (定数) は表示せず、メインの指標に絞る
    show_cols <- c("Lag", "ACF", "Xi", "Xi_Threshold_95")

    # 念のため、データフレームに存在する列だけを選ぶ (安全策)
    cols_exist <- intersect(show_cols, names(x$data))

    # 最初の5行だけ表示
    print(utils::head(x$data[, cols_exist], 5), row.names = FALSE)

    # 残りがある場合はメッセージを表示
    if (nrow(x$data) > 5) {
        cat("... with", nrow(x$data) - 5, "more lags\n")
    }

    invisible(x)
}
