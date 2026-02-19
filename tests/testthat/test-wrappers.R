test_that("xi_test function works and returns correct class", {
    # データ生成
    set.seed(42)
    x <- rnorm(100)

    # 実行
    res <- xiacf::xi_test(x, max_lag = 10, n_surr = 10)

    # クラスの確認
    expect_s3_class(res, "xi_test")

    # データフレームの構造確認
    expect_true("data.frame" %in% class(res$data))
    expect_equal(
        names(res$data),
        c("Lag", "ACF", "Xi", "Xi_Threshold_95", "ACF_CI")
    )
    expect_equal(nrow(res$data), 10)
})

test_that("xi_test handles input errors gracefully", {
    # 短すぎるデータ
    expect_error(xiacf::xi_test(c(1, 2, 3)), "too short")

    # 定数（分散ゼロ）
    expect_error(xiacf::xi_test(rep(1, 100)), "zero variance")

    # NA混入（警告が出て計算されるか）
    x_na <- c(rnorm(10), NA)
    expect_warning(res <- xiacf::xi_test(x_na, max_lag = 5), "contains NA")
    expect_s3_class(res, "xi_test")
})

test_that("print method produces output", {
    x <- rnorm(50)
    res <- xiacf::xi_test(x, max_lag = 5, n_surr = 5)

    # print実行時の出力をキャプチャして確認
    out <- capture.output(print(res))

    # 出力に特定の単語が含まれているか
    expect_true(any(grepl("Chatterjee's Xi-ACF Test", out)))
    expect_true(any(grepl("Data length:", out)))
})

test_that("run_rolling_xi_analysis works with new 'x' argument", {
    # 小さなデータでテスト
    x <- rnorm(50)

    # 逐次処理で実行 (n_cores = 1 or NULL)
    res_df <- xiacf::run_rolling_xi_analysis(
        x = x, # 引数名 x が効いているかチェック
        window_size = 20,
        step_size = 10,
        max_lag = 5,
        n_surr = 10,
        n_cores = 1 # テストなので並列化せず軽く済ませる
    )

    # 結果の型確認
    expect_true(is.data.frame(res_df))

    # ウィンドウ数の確認 (50 - 20)/10 + 1 = 4 ウィンドウ
    # 各ウィンドウで 5ラグ分あるので、合計 20行のはず
    expected_rows <- 4 * 5
    expect_equal(nrow(res_df), expected_rows)

    # 列名の確認
    expect_true("Xi_Excess" %in% names(res_df))
})
