test_that("run_rolling_xi_acf works and saves intermediate results", {
    set.seed(42)
    x_test <- rnorm(100)
    tmp_dir <- tempfile("xi_rolling_test")
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE))

    res_initial <- suppressWarnings(run_rolling_xi_acf(
        x = x_test,
        window_size = 80,
        step_size = 10,
        max_lag = 2,
        n_surr = 20,
        save_dir = tmp_dir
    ))
    expect_length(list.files(tmp_dir, pattern = "\\.rds$"), 3)
    expect_s3_class(res_initial, "data.frame")
    expect_true("Global_Threshold" %in% names(res_initial))
})

test_that("run_rolling_xi_ccf works sequentially", {
    set.seed(42)
    x <- rnorm(100)
    y <- rnorm(100)
    window_size <- 50
    step_size <- 10
    max_lag <- 2

    res <- suppressWarnings(run_rolling_xi_ccf(
        x = x,
        y = y,
        window_size = window_size,
        step_size = step_size,
        max_lag = max_lag,
        n_surr = 20
    ))

    expect_s3_class(res, "data.frame")
    expected_windows <- length(seq(1, 100 - window_size + 1, by = step_size))
    expected_rows_per_window <- 2 * (max_lag + 1) # Lag 0起点の双方向
    expect_equal(nrow(res), expected_windows * expected_rows_per_window)
})

test_that("Rolling functions strictly preserve Window_ID for traceability", {
    set.seed(42)
    x <- rnorm(100)
    y <- rnorm(100)

    # Check univariate rolling
    res_acf <- suppressWarnings(xiacf::run_rolling_xi_acf(
        x,
        window_size = 50,
        step_size = 10,
        max_lag = 2,
        n_surr = 50
    ))
    expect_true("Window_ID" %in% colnames(res_acf))
    expect_true(is.numeric(res_acf$Window_ID))
    # Starting at 1, step 10 -> indices 1, 11, 21, 31, 41, 51 (Total 6 windows)
    expect_equal(length(unique(res_acf$Window_ID)), 6)

    # Check bivariate rolling
    res_ccf <- suppressWarnings(xiacf::run_rolling_xi_ccf(
        x,
        y,
        window_size = 50,
        step_size = 10,
        max_lag = 2,
        n_surr = 50
    ))
    expect_true("Window_ID" %in% colnames(res_ccf))
    expect_true(is.numeric(res_ccf$Window_ID))
})
