test_that("xi_matrix input validation works", {
    expect_error(xi_matrix(c(1, 2, 3)), "numeric matrix or data.frame")
    expect_error(
        xi_matrix(data.frame(A = letters[1:10], B = 1:10)),
        "must be numeric"
    )
    df_na <- data.frame(A = c(1, NA, 3, 4, 5), B = 1:5)
    expect_error(xi_matrix(df_na), "NA values")
    expect_error(
        xi_matrix(matrix(rnorm(10), ncol = 1)),
        "at least 2 columns"
    )

    expect_error(
        xi_matrix(data.frame(A = 1:10, B = 1:10), max_lag = 2, sig_level = 1.5),
        "strictly between 0 and 1"
    )
})

test_that("xi_matrix computes correctly and returns S3 object", {
    set.seed(42)
    df <- data.frame(A = rnorm(30), B = rnorm(30), C = rnorm(30))
    res <- suppressWarnings(xi_matrix(df, max_lag = 3, n_surr = 20))

    expect_s3_class(res, "xi_matrix")
    expect_equal(res$sig_level, 0.05)

    expect_equal(nrow(res$data), 27)
    expect_true(all(
        c(
            "Lead_Var",
            "Lag_Var",
            "Lag",
            "Xi",
            "Global_Threshold",
            "Xi_Excess"
        ) %in%
            colnames(res$data)
    ))
})

test_that("autoplot.xi_matrix works", {
    set.seed(42)
    df <- data.frame(A = rnorm(20), B = rnorm(20))
    res <- suppressWarnings(xi_matrix(df, max_lag = 2, n_surr = 20))
    expect_s3_class(autoplot(res), "ggplot")
})
