test_that("xi_acf function works and returns correct structure", {
  # Generate synthetic data
  set.seed(42)
  x <- rnorm(100)
  max_lag <- 10
  n_surr <- 10

  # Execute the main function
  res <- xiacf::xi_acf(x, max_lag = max_lag, n_surr = n_surr)

  # Verify the returned S3 class
  expect_s3_class(res, "xi_acf")
  expect_equal(res$max_lag, max_lag)
  expect_equal(res$n_surr, n_surr)
  expect_equal(res$n, length(x))

  # Verify the structure and contents of the data frame
  expect_true("data.frame" %in% class(res$data))
  expected_cols <- c("Lag", "ACF", "Xi", "Xi_Threshold_95", "ACF_CI")
  expect_true(all(expected_cols %in% names(res$data)))
  expect_equal(nrow(res$data), max_lag)
})

test_that("xi_acf handles input errors gracefully", {
  # Edge case 1: Time series is too short
  expect_error(xiacf::xi_acf(c(1, 2, 3)), "too short")

  # Edge case 2: Constant series (zero variance)
  expect_error(xiacf::xi_acf(rep(1, 100)), "zero variance")

  # Edge case 3: Handling of NA values (should warn, remove NA, and compute)
  x_na <- c(rnorm(10), NA)
  expect_warning(
    res <- xiacf::xi_acf(x_na, max_lag = 5, n_surr = 5),
    "contains NA"
  )
  expect_s3_class(res, "xi_acf")
})

test_that("xi_test wrapper throws deprecation warning but returns xi_acf object", {
  set.seed(42)
  x <- rnorm(100)

  # Verify that the wrapper function throws a deprecation warning and returns the correct object type
  expect_warning(
    res <- xiacf::xi_test(x, max_lag = 5, n_surr = 5),
    "deprecated"
  )
  expect_s3_class(res, "xi_acf")
})

test_that("print method produces correct output", {
  x <- rnorm(50)
  res <- xiacf::xi_acf(x, max_lag = 5, n_surr = 5)

  # Capture the output of the print method
  out <- capture.output(print(res))

  # Verify that specific keywords are present in the console output
  expect_true(any(grepl("Chatterjee's Xi-ACF Test", out)))
  expect_true(any(grepl("Data length:", out)))
})

test_that("autoplot.xi_acf returns a ggplot object", {
  set.seed(42)
  x <- rnorm(100)
  res <- xiacf::xi_acf(x, max_lag = 5, n_surr = 5)

  # Verify that autoplot.xi_acf returns a ggplot object
  p <- ggplot2::autoplot(res)
  expect_s3_class(p, "ggplot")
})

test_that("run_rolling_xi_analysis works with sequential execution", {
  # Test with a small synthetic dataset
  set.seed(42)
  x <- rnorm(50)

  # Execute sequentially to keep the automated test lightweight
  res_df <- xiacf::run_rolling_xi_analysis(
    x = x,
    window_size = 20,
    step_size = 10,
    max_lag = 5,
    n_surr = 10,
    n_cores = 1 # Avoid parallel overhead in CRAN checks
  )

  # Verify the type of the result
  expect_s3_class(res_df, "data.frame")
  expect_true("Xi_Excess" %in% names(res_df))
  expect_true(nrow(res_df) > 0)
})
