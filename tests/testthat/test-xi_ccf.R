# tests/testthat/test-xi_ccf.R

# tests/testthat/test-xi_ccf.R

test_that("xi_ccf computes correctly and returns proper structure", {
    set.seed(42)
    x <- rnorm(100)
    y <- rnorm(100)
    max_lag <- 3

    # Test for direction = "both"
    res_both <- suppressWarnings(xi_ccf(
        x,
        y,
        max_lag = max_lag,
        n_surr = 20,
        direction = "both"
    ))
    expect_s3_class(res_both, "xi_ccf")
    expect_equal(nrow(res_both$data), 2 * (max_lag + 1)) # 2*(max_lag+1) rows because it includes Lag 0

    # Check for mandatory columns in tidy format
    expected_cols <- c(
        "Lead_Var",
        "Lag_Var",
        "Lag",
        "Xi",
        "CCF",
        "Global_Threshold",
        "Xi_Excess",
        "CCF_CI"
    )
    expect_true(all(expected_cols %in% colnames(res_both$data)))

    # Test for direction = "x_leads"
    res_x <- suppressWarnings(xi_ccf(
        x,
        y,
        max_lag = max_lag,
        n_surr = 20,
        direction = "x_leads"
    ))
    expect_equal(nrow(res_x$data), max_lag + 1)
    expect_true(all(res_x$data$Lead_Var == "x" & res_x$data$Lag_Var == "y"))
}) # ← ★先ほど誤って消去されてしまったのはこの行です

test_that("xi_ccf handles invalid inputs correctly", {
    set.seed(42)
    x <- rnorm(100)
    y <- rnorm(100)

    # Different lengths
    expect_error(xi_ccf(x[1:90], y), "same length")
    # Contains NA
    expect_error(xi_ccf(c(x[1:99], NA), y), "NA values")
    # Zero variance
    expect_error(xi_ccf(rep(1, 100), y), "non-zero variance")
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
