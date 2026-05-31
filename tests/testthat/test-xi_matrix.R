# tests/testthat/test-xi_matrix.R

test_that("xi_matrix input validation works", {
    expect_error(xi_matrix(c(1, 2, 3)), "numeric matrix or data.frame")
    expect_error(
        xi_matrix(data.frame(A = letters[1:10], B = 1:10)),
        "completely numeric"
    )
    expect_error(
        xi_matrix(data.frame(A = c(1, NA, 3, 4, 5), B = 1:5)),
        "NA values"
    )
    expect_error(xi_matrix(matrix(rnorm(10), ncol = 1)), "at least 2 columns")
})

test_that("xi_matrix computes correctly and returns proper structure", {
    set.seed(42)
    df <- data.frame(A = rnorm(30), B = rnorm(30), C = rnorm(30))
    p <- ncol(df)
    max_lag <- 2

    res <- suppressWarnings(xi_matrix(df, max_lag = max_lag, n_surr = 20))

    expect_s3_class(res, "xi_matrix")

    # Because it includes Lag 0, the total number of rows is p * p * (max_lag + 1)
    expect_equal(nrow(res$data), p * p * (max_lag + 1))

    # Check if data_raw is saved correctly
    expect_identical(res$data_raw, df)

    # Check if threshold and Excess for diagonal elements (autocorrelation) are NA (confirming self-loop exclusion)
    diags <- res$data[res$data$Lead_Var == res$data$Lag_Var, ]
    expect_true(all(is.na(diags$Global_Threshold)))
    expect_true(all(is.na(diags$Xi_Excess)))

    # Check if threshold is calculated for off-diagonal elements (cross-correlation)
    off_diags <- res$data[res$data$Lead_Var != res$data$Lag_Var, ]
    expect_true(all(!is.na(off_diags$Global_Threshold)))
    expect_true(all(!is.na(off_diags$Xi_Excess)))
})

test_that("xi_matrix properly handles Dual-Family FWER and completely excludes autocorrelations", {
    set.seed(42)
    df_input <- data.frame(A = rnorm(50), B = rnorm(50), C = rnorm(50))

    res <- suppressWarnings(xiacf::xi_matrix(
        df_input,
        max_lag = 2,
        n_surr = 100,
        sig_level = 0.05
    ))

    df <- res$data

    # 1. Autocorrelation Exclusion (Must be NA to prevent FWER inflation)
    auto_corr <- df[df$Lead_Var == df$Lag_Var, ]
    expect_true(all(is.na(auto_corr$Global_Threshold)))
    expect_true(all(is.na(auto_corr$Xi_Excess)))

    # 2. Dual-Family Thresholds for Cross-correlation
    cross_corr <- df[df$Lead_Var != df$Lag_Var, ]
    thresh_lag0 <- cross_corr$Global_Threshold[cross_corr$Lag == 0][1]
    thresh_lagged <- cross_corr$Global_Threshold[cross_corr$Lag > 0][1]

    expect_false(is.na(thresh_lag0))
    expect_false(is.na(thresh_lagged))

    # Thresholds must be distinct due to the Dual-Family separation policy
    expect_true(thresh_lag0 != thresh_lagged)
})
