test_that("xi_acf function works and returns correct structure", {
  set.seed(42)
  x <- rnorm(100)
  max_lag <- 10
  res <- suppressWarnings(xiacf::xi_acf(x, max_lag = max_lag, n_surr = 20))

  expect_s3_class(res, "xi_acf")
  expected_cols <- c(
    "Lag",
    "ACF",
    "Xi",
    "Global_Threshold",
    "ACF_CI",
    "Xi_Excess"
  )
  expect_true(all(expected_cols %in% names(res$data)))
})

test_that("xi_acf handles input errors gracefully", {
  expect_error(xiacf::xi_acf(c(1, 2, 3)), "too short")
  expect_error(xiacf::xi_acf(rep(1, 100)), "zero variance")
  expect_error(xiacf::xi_acf(c(rnorm(10), NA), max_lag = 5), "NA values")
})

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

test_that("sig_level dynamically adjusts thresholds and CIs", {
  set.seed(42)
  x <- rnorm(100)

  res_05 <- suppressWarnings(xiacf::xi_acf(
    x,
    max_lag = 5,
    n_surr = 100,
    sig_level = 0.05
  ))
  res_01 <- suppressWarnings(xiacf::xi_acf(
    x,
    max_lag = 5,
    n_surr = 100,
    sig_level = 0.01
  ))

  expect_gt(res_01$data$ACF_CI[1], res_05$data$ACF_CI[1])
  expect_gt(
    mean(res_01$data$Global_Threshold),
    mean(res_05$data$Global_Threshold)
  )
})
