test_that("xi_coefficient detects linear and non-linear relationships", {
    n <- 100
    x <- 1:n
    y <- x
    xi_val <- xiacf::xi_coefficient(x, y)
    expect_gt(xi_val, 0.96)

    x_para <- seq(-10, 10, length.out = 100)
    y_para <- x_para^2
    xi_val_para <- xiacf::xi_coefficient(x_para, y_para)
    expect_gt(xi_val_para, 0.5)
})

test_that("Xi coefficient is close to 0 for independent noise", {
    set.seed(123)
    n <- 500
    x <- rnorm(n)
    y <- rnorm(n)
    xi_val_noise <- xiacf::xi_coefficient(x, y)
    expect_lt(abs(xi_val_noise), 0.35)
})

test_that("compute_xi_acf_maxstat_cpp throws error when max_lag is too large", {
    x_short <- rnorm(10)
    max_lag <- 20
    n_surr <- 20

    expect_error(
        xiacf:::compute_xi_acf_maxstat_cpp(x_short, max_lag, n_surr, 10L),
        "strictly greater"
    )
})
