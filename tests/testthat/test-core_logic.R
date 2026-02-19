test_that("Xi coefficient captures functional dependencies correctly", {
    # Case 1: 完全な線形関係 (y = x)
    # 理論上、nが増えれば 1.0 に収束するはず (厳密には 1 - 3/(n+1) 付近)
    n <- 100
    x <- 1:n
    y <- x

    # xi_coefficient は RcppExports で export されているはず
    xi_val <- xiacf::xi_coefficient(x, y)

    # 完全に1にはならない仕様なので、0.96以上なら合格とする
    expect_gt(xi_val, 0.96)

    # Case 2: 完全なU字型 (y = x^2)
    # ピアソン相関なら0になるが、Xiなら高い値が出るはず
    x_para <- seq(-10, 10, length.out = 100)
    y_para <- x_para^2

    xi_val_para <- xiacf::xi_coefficient(x_para, y_para)

    # これが「非線形検知」の証拠。0.5以上なら合格（実際はかなり高くなる）
    expect_gt(xi_val_para, 0.5)
})

test_that("Xi coefficient is close to 0 for independent noise", {
    # Case 3: ランダムノイズ (独立)
    set.seed(123)
    n <- 500
    x <- rnorm(n)
    y <- rnorm(n)

    xi_val_noise <- xiacf::xi_coefficient(x, y)

    # 0付近であること（絶対値が0.1以下なら合格）
    expect_lt(abs(xi_val_noise), 0.1)
})

test_that("compute_xi_lags handles initialization correctly (NA check)", {
    # Case 4: 初期化漏れバグの回帰テスト
    # データ長(n=10) より ラグ(max_lag=20) が大きい意地悪なケース
    x_short <- rnorm(10)
    max_lag <- 20
    n_surr <- 5

    # エラー落ちせずに結果が返ってくるか？
    res <- xiacf::compute_xi_lags(x_short, max_lag, n_surr)

    # 計算できないラグ（11以降）は NaN (または NA) になっているか確認
    # Rcppの datum::nan は R では NaN として扱われる
    expect_true(is.nan(res$xi_original[15]))
    expect_true(all(is.nan(res$xi_surrogates[15, ])))

    # 計算できるラグ（1〜9）は数値が入っているか
    expect_false(is.nan(res$xi_original[5]))
})
