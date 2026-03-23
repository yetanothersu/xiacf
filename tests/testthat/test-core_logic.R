test_that("xi_coefficient detects linear and non-linear relationships", {
    # Case 1: Perfect linear relationship (y = x)
    # Theoretically, it converges to 1.0 as n increases (strictly around 1 - 3/(n+1)).
    n <- 100
    x <- 1:n
    y <- x

    # xi_coefficient is exported via RcppExports
    xi_val <- xiacf::xi_coefficient(x, y)

    # Since it does not become exactly 1.0 by design, we accept values > 0.96.
    expect_gt(xi_val, 0.96)

    # Case 2: Perfect U-shaped relationship (y = x^2)
    # Pearson correlation would be 0, but Xi should yield a high value.
    x_para <- seq(-10, 10, length.out = 100)
    y_para <- x_para^2

    xi_val_para <- xiacf::xi_coefficient(x_para, y_para)

    # This is evidence of "non-linear detection".
    # We accept values > 0.5 (it will actually be much higher).
    expect_gt(xi_val_para, 0.5)
})

test_that("Xi coefficient is close to 0 for independent noise", {
    # Case 3: Random noise (independent variables)
    set.seed(123)
    n <- 500
    x <- rnorm(n)
    y <- rnorm(n)

    xi_val_noise <- xiacf::xi_coefficient(x, y)

    # Should be close to 0 (absolute value < 0.1 is acceptable)
    expect_lt(abs(xi_val_noise), 0.1)
})

test_that("compute_xi_lags handles initialization correctly (NA check)", {
    # Case 4: Regression test for initialization bug
    # Edge case where max_lag (20) is greater than the data length (10).
    x_short <- rnorm(10)
    max_lag <- 20
    n_surr <- 5

    # Ensure it returns a result without throwing an error
    res <- xiacf::compute_xi_lags(x_short, max_lag, n_surr)

    # Check if uncomputable lags (>= 11) are filled with NaN.
    # Note: Rcpp's datum::nan is treated as NaN in R.
    expect_true(is.nan(res$xi_original[15]))
    expect_true(all(is.nan(res$xi_surrogates[15, ])))

    # Check if computable lags (1 to 9) contain valid numeric values.
    expect_false(is.nan(res$xi_original[5]))
})
