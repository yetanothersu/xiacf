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
})

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

test_that("Dual-Family FWER Control separates Lag 0 and Lag > 0 thresholds", {
    set.seed(42)
    x <- rnorm(100)
    y <- rnorm(100)

    # Suppress warnings for low surrogate count to speed up tests
    res <- suppressWarnings(xiacf::xi_ccf(
        x,
        y,
        max_lag = 3,
        n_surr = 100,
        sig_level = 0.05
    ))

    df <- res$data

    thresh_lag0 <- df$Global_Threshold[df$Lag == 0][1]
    thresh_lagged <- df$Global_Threshold[df$Lag > 0][1]

    expect_false(is.na(thresh_lag0))
    expect_false(is.na(thresh_lagged))

    # Because MIAAFT preserves Lag 0 structures, the threshold for Lag 0
    # should be calculated independently and generally differs from Lag > 0.
    expect_true(thresh_lag0 != thresh_lagged)
})
