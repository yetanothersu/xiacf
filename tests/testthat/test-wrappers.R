test_that("xi_test function works and returns correct class", {
    # Generate synthetic data
    set.seed(42)
    x <- rnorm(100)

    # Execute the main function
    res <- xiacf::xi_test(x, max_lag = 10, n_surr = 10)

    # Verify the returned S3 class
    expect_s3_class(res, "xi_test")

    # Verify the structure and contents of the data frame
    expect_true("data.frame" %in% class(res$data))
    expect_equal(
        names(res$data),
        c("Lag", "ACF", "Xi", "Xi_Threshold_95", "ACF_CI")
    )
    expect_equal(nrow(res$data), 10)
})

test_that("xi_test handles input errors gracefully", {
    # Edge case 1: Time series is too short
    expect_error(xiacf::xi_test(c(1, 2, 3)), "too short")

    # Edge case 2: Constant series (zero variance)
    expect_error(xiacf::xi_test(rep(1, 100)), "zero variance")

    # Edge case 3: Handling of NA values (should warn, remove NA, and compute)
    x_na <- c(rnorm(10), NA)
    expect_warning(res <- xiacf::xi_test(x_na, max_lag = 5), "contains NA")
    expect_s3_class(res, "xi_test")
})

test_that("print method produces output", {
    x <- rnorm(50)
    res <- xiacf::xi_test(x, max_lag = 5, n_surr = 5)

    # Capture the output of the print method
    out <- capture.output(print(res))

    # Verify that specific keywords are present in the console output
    expect_true(any(grepl("Chatterjee's Xi-ACF Test", out)))
    expect_true(any(grepl("Data length:", out)))
})

test_that("run_rolling_xi_analysis works with new 'x' argument", {
    # Test with a small synthetic dataset
    x <- rnorm(50)

    # Execute sequentially to keep the automated test lightweight
    res_df <- xiacf::run_rolling_xi_analysis(
        x = x, # Ensure the argument 'x' is correctly mapped
        window_size = 20,
        step_size = 10,
        max_lag = 5,
        n_surr = 10,
        n_cores = 1 # Avoid parallel overhead in CRAN checks
    )

    # Verify the type of the result
    expect_true(is.data.frame(res_df))

    # Check the total number of windows: (50 - 20) / 10 + 1 = 4 windows
    # Each window computes 5 lags, resulting in exactly 4 * 5 = 20 rows
    expected_rows <- 4 * 5
    expect_equal(nrow(res_df), expected_rows)

    # Verify the presence of specific computed columns
    expect_true("Xi_Excess" %in% names(res_df))
})
