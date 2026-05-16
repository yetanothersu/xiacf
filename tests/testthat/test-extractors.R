set.seed(42)
n_test <- 50
df_test <- data.frame(A = rnorm(n_test), B = rnorm(n_test), C = rnorm(n_test))

mat_obj <- suppressWarnings(xi_matrix(df_test, max_lag = 3, n_surr = 20))

test_that("extract_xi_acf works correctly with and without x_raw", {
    acf_only <- extract_xi_acf(mat_obj, var = "A")
    expect_s3_class(acf_only, "xi_acf")
    expect_equal(acf_only$max_lag, 3)
    expect_true(all(is.na(acf_only$data$ACF)))

    acf_full <- extract_xi_acf(mat_obj, var = "A", x_raw = df_test)
    expect_false(all(is.na(acf_full$data$ACF)))
})

test_that("extract_xi_ccf works correctly with and without x_raw", {
    ccf_only <- extract_xi_ccf(mat_obj, var_x = "A", var_y = "B")
    expect_s3_class(ccf_only, "xi_ccf")
    expect_true(any(ccf_only$data$Lag < 0))
    expect_true(any(ccf_only$data$Lag > 0))
    expect_true(all(is.na(ccf_only$data$CCF)))

    ccf_full <- extract_xi_ccf(
        mat_obj,
        var_x = "A",
        var_y = "B",
        x_raw = df_test
    )
    expect_false(all(is.na(ccf_full$data$CCF)))
})

test_that("extractor functions handle errors appropriately", {
    expect_error(extract_xi_acf(df_test, "A"), "must be a 'xi_matrix' object")
    expect_error(extract_xi_acf(mat_obj, var = "Z"), "not found")
})

test_that("extracted objects can be plotted without errors", {
    acf_full <- extract_xi_acf(mat_obj, var = "A", x_raw = df_test)
    ccf_only <- extract_xi_ccf(mat_obj, var_x = "A", var_y = "B")
    expect_s3_class(ggplot2::autoplot(acf_full), "ggplot")
    expect_s3_class(ggplot2::autoplot(ccf_only), "ggplot")
})
