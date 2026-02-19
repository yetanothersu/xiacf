#' Rolling Xi-ACF Analysis
#'
#' Performs a rolling window analysis using Chatterjee's Xi coefficient.
#'
#' @param x A numeric vector of the time series (e.g., log-returns).
#' @param window_size Integer. The size of the rolling window.
#' @param step_size Integer. The step size for moving the window.
#' @param max_lag Integer. Maximum lag to calculate Xi for.
#' @param n_surr Integer. Number of surrogates for the null hypothesis test.
#' @param n_cores Integer (optional). Number of cores.
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom future plan multisession sequential
#' @importFrom progressr progressor
#' @importFrom dplyr bind_rows
#' @export
run_rolling_xi_analysis <- function(
    x, # ts_vec -> x に変更
    window_size,
    step_size = 1,
    max_lag = 20,
    n_surr = 100,
    n_cores = NULL
) {
    # --- 1. Input Validation ---
    n_total <- length(x)
    if (window_size > n_total) {
        stop("window_size cannot be larger than the time series length.")
    }

    # --- 2. Safe Parallel Setup (Polite Programming) ---
    # ユーザーが明示的にコア数を指定した場合のみプランを変更し、終わったら戻す
    if (!is.null(n_cores)) {
        # CRANチェック対策: コア数制限がある環境では最大2コアに抑える
        chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
        if (nzchar(chk) && chk == "TRUE") {
            n_cores <- min(n_cores, 2)
        }

        # 現在のプランを保存
        old_plan <- future::plan()

        # 新しいプランを設定
        future::plan(future::multisession, workers = n_cores)

        # 関数終了時(エラー時含む)に必ず元のプランに戻す
        on.exit(future::plan(old_plan), add = TRUE)
    }

    # --- 3. Prepare Windows ---
    starts <- seq(1, n_total - window_size + 1, by = step_size)
    n_windows <- length(starts)

    p <- progressr::progressor(steps = n_windows)

    # --- 4. Execution ---
    results_df <- foreach::foreach(
        i = seq_along(starts),
        .combine = dplyr::bind_rows,
        .packages = c("xiacf", "stats"),
        .options.future = list(seed = TRUE),
        .errorhandling = "remove"
    ) %dopar%
        {
            p()

            tryCatch(
                {
                    idx_start <- starts[i]
                    idx_end <- idx_start + window_size - 1

                    # ts_vec -> x
                    y_sub <- x[idx_start:idx_end]

                    # 定数チェック
                    if (sd(y_sub, na.rm = TRUE) == 0) {
                        return(NULL)
                    }

                    # C++エンジンの呼び出し (後でリネームするならここも変更)
                    res <- compute_xi_lags(y_sub, max_lag, n_surr)

                    # 閾値計算 (NA除去必須)
                    xi_threshold <- rep(NA, max_lag)
                    if (n_surr > 0) {
                        xi_threshold <- apply(
                            res$xi_surrogates,
                            1,
                            function(r) {
                                stats::quantile(r, 0.95, na.rm = TRUE)
                            }
                        )
                    }

                    # データフレーム構築
                    data.frame(
                        Window_ID = i, # 何番目のウィンドウか
                        Window_Start_Idx = idx_start,
                        Window_End_Idx = idx_end,
                        Lag = 1:max_lag,
                        Xi_Original = as.numeric(res$xi_original),
                        Xi_Threshold_95 = xi_threshold,
                        # 余剰量 (有意でなければ0にする処理を入れても良いが、生の値を入れておく)
                        Xi_Excess = as.numeric(res$xi_original) - xi_threshold
                    )
                },
                error = function(e) {
                    # エラー時は警告を出して NULL を返す (bind_rowsで無視される)
                    warning(paste("Error in window", i, ":", e$message))
                    return(NULL)
                }
            )
        }

    return(results_df)
}
