# Generate lightweight test data and compute the matrix
set.seed(42)
n_test <- 50
df_test <- data.frame(
    A = rnorm(n_test),
    B = rnorm(n_test),
    C = rnorm(n_test)
)
mat_obj <- xi_matrix(df_test, max_lag = 3, n_surr = 10)

test_that("extract_xi_acf works correctly with and without x_raw", {
    # 1. Without raw data (Xi metrics only)
    acf_only <- extract_xi_acf(mat_obj, var = "A")

    expect_s3_class(acf_only, "xi_acf")
    expect_equal(acf_only$max_lag, 3)
    expect_equal(acf_only$n, n_test)
    # Ensure all standard linear metrics are NA
    expect_true(all(is.na(acf_only$data$ACF)))
    expect_true(all(is.na(acf_only$data$ACF_CI)))

    # 2. With raw data (Recalculate linear ACF)
    acf_full <- extract_xi_acf(mat_obj, var = "A", x_raw = df_test)

    expect_s3_class(acf_full, "xi_acf")
    # Ensure standard ACF metrics are properly populated (not NA)
    expect_false(all(is.na(acf_full$data$ACF)))
    expect_false(all(is.na(acf_full$data$ACF_CI)))
})

test_that("extract_xi_ccf works correctly with and without x_raw", {
    # 1. Without raw data (Xi metrics only)
    ccf_only <- extract_xi_ccf(mat_obj, var_x = "A", var_y = "B")

    expect_s3_class(ccf_only, "xi_ccf")
    expect_true(ccf_only$bidirectional)
    # Ensure both forward and backward directions are combined
    expect_true(all(c("X leads Y", "Y leads X") %in% ccf_only$data$Direction))
    expect_true(all(is.na(ccf_only$data$CCF)))

    # 2. With raw data (Recalculate linear CCF)
    ccf_full <- extract_xi_ccf(
        mat_obj,
        var_x = "A",
        var_y = "B",
        x_raw = df_test
    )

    expect_s3_class(ccf_full, "xi_ccf")
    expect_false(all(is.na(ccf_full$data$CCF)))
})

test_that("extractor functions handle errors appropriately", {
    # Invalid object class
    expect_error(extract_xi_acf(df_test, "A"), "must be a 'xi_matrix' object")
    expect_error(
        extract_xi_ccf(df_test, "A", "B"),
        "must be a 'xi_matrix' object"
    )

    # Non-existent variables
    expect_error(extract_xi_acf(mat_obj, var = "Z"), "not found")
    expect_error(extract_xi_ccf(mat_obj, var_x = "A", var_y = "Z"), "not found")
})

test_that("extracted objects can be plotted without errors", {
    # Check if autoplot runs successfully and returns a ggplot object.
    # This also implicitly verifies that no unexpected warnings (e.g., from NAs) are thrown.

    acf_full <- extract_xi_acf(mat_obj, var = "A", x_raw = df_test)
    ccf_only <- extract_xi_ccf(mat_obj, var_x = "A", var_y = "B")

    p_acf <- ggplot2::autoplot(acf_full)
    p_ccf <- ggplot2::autoplot(ccf_only)

    expect_s3_class(p_acf, "ggplot")
    expect_s3_class(p_ccf, "ggplot")
})
