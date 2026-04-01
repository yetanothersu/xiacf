# tests/testthat/test-xi_ccf.R

test_that("xi_ccf computes correctly (bidirectional = TRUE) and returns proper structure", {
    set.seed(42)
    x <- rnorm(100)
    y <- rnorm(100)

    # Run with small parameters for testing
    max_lag <- 3
    res <- xi_ccf(x, y, max_lag = max_lag, n_surr = 10, bidirectional = TRUE)

    # Check S3 class and metadata
    expect_s3_class(res, "xi_ccf")
    expect_equal(res$n, 100)
    expect_equal(res$max_lag, max_lag)
    expect_true(res$bidirectional)

    # Check data frame dimensions
    # 0 to max_lag (max_lag + 1 rows) * 2 directions (X leads Y, Y leads X)
    expect_equal(nrow(res$data), 2 * (max_lag + 1))

    # Check new long-format columns
    expect_true(all(
        c("Direction", "Lag", "CCF", "Xi", "Xi_Threshold_95", "CCF_CI") %in%
            colnames(res$data)
    ))

    # Check if both directions are correctly included
    expect_setequal(unique(res$data$Direction), c("X leads Y", "Y leads X"))
})

test_that("xi_ccf computes correctly (bidirectional = FALSE)", {
    set.seed(42)
    x <- rnorm(100)
    y <- rnorm(100)

    max_lag <- 3
    res <- xi_ccf(x, y, max_lag = max_lag, n_surr = 10, bidirectional = FALSE)

    # Check data frame dimensions (unidirectional only)
    expect_equal(nrow(res$data), max_lag + 1)
    expect_setequal(unique(res$data$Direction), "X leads Y")
})

test_that("xi_ccf handles invalid inputs correctly", {
    set.seed(42)
    x <- rnorm(100)
    y <- rnorm(100)

    # When lengths differ
    expect_error(xi_ccf(x, y[1:50]), "exact same length")

    # When NA is included
    x_na <- x
    x_na[10] <- NA
    expect_error(xi_ccf(x_na, y), "NA values")

    # When variance is zero (ensure internal stats::ccf throws an error)
    expect_error(xi_ccf(rep(1, 100), y))
})

test_that("autoplot.xi_ccf returns a valid ggplot object", {
    set.seed(42)
    x <- rnorm(50)
    y <- rnorm(50)
    res <- xi_ccf(x, y, max_lag = 2, n_surr = 2)

    p <- ggplot2::autoplot(res)
    expect_s3_class(p, "ggplot")
})

# TODO: Temporarily commented out until Phase 3 adaptation (long-format adaptation)
# for xi_rolling_ccf.R is complete.

# test_that("run_rolling_xi_ccf works sequentially", {
#     set.seed(42)
#     x <- rnorm(100)
#     y <- rnorm(100)
#
#     window_size <- 50
#     step_size <- 10
#     max_lag <- 2
#
#     res <- run_rolling_xi_ccf(
#         x,
#         y,
#         window_size = window_size,
#         step_size = step_size,
#         max_lag = max_lag,
#         n_surr = 5,
#         n_cores = NULL
#     )
#
#     # クラスと列の確認
#     expect_s3_class(res, "data.frame")
#     expect_true(all(
#         c("Window_ID", "Window_Start_Idx", "Lag", "Xi_Excess") %in%
#             colnames(res)
#     ))
#
#     # 行数の確認 (ウィンドウ数 × ラグ数)
#     expected_windows <- length(seq(1, 100 - window_size + 1, by = step_size))
#     expected_lags <- 2 * max_lag + 1
#     expect_equal(nrow(res), expected_windows * expected_lags)
# })
