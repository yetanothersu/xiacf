# tests/testthat/test-xi_matrix.R

test_that("xi_matrix input validation works", {
    # 不正な入力（ベクトルや文字列）
    expect_error(xi_matrix(c(1, 2, 3)), "numeric matrix or data.frame")
    expect_error(
        xi_matrix(data.frame(A = letters[1:10], B = 1:10)),
        "must be numeric"
    )

    # NAのチェック
    df_na <- data.frame(A = c(1, NA, 3, 4, 5), B = 1:5)
    expect_error(xi_matrix(df_na), "NA values")

    # 変数数のチェック
    expect_error(
        xi_matrix(matrix(rnorm(10), ncol = 1)),
        "at least two variables"
    )

    # 有意水準のチェック
    expect_error(
        xi_matrix(data.frame(A = 1:10, B = 1:10), sig_level = 1.5),
        "strictly between 0 and 1"
    )
})

test_that("xi_matrix computes correctly and returns S3 object", {
    set.seed(42)
    # テスト用に小さなデータセットを作成
    df <- data.frame(A = rnorm(30), B = rnorm(30), C = rnorm(30))

    # 計算量削減のため max_lag と n_surr を小さめに設定
    res <- xi_matrix(df, max_lag = 3, n_surr = 5)

    # クラスと基本情報の確認
    expect_s3_class(res, "xi_matrix")
    expect_equal(res$M, 3)
    expect_equal(res$max_lag, 3)
    expect_equal(res$n_surr, 5)
    expect_equal(res$sig_level, 0.95) # デフォルト値

    # データフレームの行数確認: 3変数 x 3変数 x (3ラグ + 1) = 36行
    expect_equal(nrow(res$data), 36)

    # 出力カラムの確認
    expect_true(all(
        c("Lead_Var", "Lag_Var", "Lag", "Xi", "Xi_Threshold", "Xi_Excess") %in%
            colnames(res$data)
    ))
})

test_that("autoplot.xi_matrix works", {
    set.seed(42)
    df <- data.frame(A = rnorm(20), B = rnorm(20))
    res <- xi_matrix(df, max_lag = 2, n_surr = 2)

    p <- autoplot(res)
    expect_s3_class(p, "ggplot")
})
