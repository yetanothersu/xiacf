test_that("xi_ccf computes correctly and returns proper structure", {
    set.seed(42)
    x <- rnorm(100)
    y <- rnorm(100)
    max_lag <- 3
    res <- suppressWarnings(xi_ccf(x, y, max_lag = max_lag, n_surr = 20))

    expect_s3_class(res, "xi_ccf")
    expect_equal(nrow(res$data), 2 * max_lag + 1)
    expect_true(all(
        c("Lag", "CCF", "Xi", "Global_Threshold", "CCF_CI", "Xi_Excess") %in%
            colnames(res$data)
    ))
})

test_that("xi_ccf handles invalid inputs correctly", {
    set.seed(42)
    x <- rnorm(100)
    y <- rnorm(100)

    expect_error(xi_ccf(x, y[1:50]), "exactly the same")
    expect_error(xi_ccf(c(x[1:99], NA), y), "NA values")
    expect_error(xi_ccf(rep(1, 100), y))
})

test_that("run_rolling_xi_ccf works sequentially", {
    set.seed(42)
    x <- rnorm(100)
    y <- rnorm(100)
    window_size <- 50
    step_size <- 10
    max_lag <- 2

    res <- suppressWarnings(run_rolling_xi_ccf(
        x,
        y,
        window_size = window_size,
        step_size = step_size,
        max_lag = max_lag,
        n_surr = 20
    ))

    expect_s3_class(res, "data.frame")
    expected_windows <- length(seq(1, 100 - window_size + 1, by = step_size))
    expected_rows_per_window <- 2 * max_lag + 1
    expect_equal(nrow(res), expected_windows * expected_rows_per_window)
})
