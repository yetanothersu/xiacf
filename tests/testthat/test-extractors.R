# tests/testthat/test-extractors.R

test_that("extract_xi_acf and extract_xi_ccf work correctly and clean variable names", {
    set.seed(42)
    n_test <- 50
    df_test <- data.frame(
        X = rnorm(n_test),
        Y = rnorm(n_test),
        Z = rnorm(n_test)
    )

    # Calculation of the parent xi_matrix
    res_mat <- suppressWarnings(xi_matrix(df_test, max_lag = 2, n_surr = 20))

    # 1. Test for extract_xi_acf
    ext_acf <- suppressWarnings(extract_xi_acf(res_mat, var = "X"))
    expect_s3_class(ext_acf, "xi_acf")
    # Check if variable names are cleaned up correctly
    expect_equal(ext_acf$x_name, "X")

    # 2. Test for extract_xi_ccf
    ext_ccf <- suppressWarnings(extract_xi_ccf(
        res_mat,
        var_x = "X",
        var_y = "Y",
        direction = "both"
    ))
    expect_s3_class(ext_ccf, "xi_ccf")

    # Check cleanup of variable names (metadata)
    expect_equal(ext_ccf$x_name, "X")
    expect_equal(ext_ccf$y_name, "Y")

    # Check cleanup of variable names (Lead_Var, Lag_Var) in the tidy data frame
    expect_true(all(ext_ccf$data$Lead_Var %in% c("X", "Y")))
    expect_true(all(ext_ccf$data$Lag_Var %in% c("X", "Y")))
})

test_that("extractors throw errors for invalid inputs", {
    set.seed(42)
    df_test <- data.frame(X = rnorm(30), Y = rnorm(30))
    res_mat <- suppressWarnings(xi_matrix(df_test, max_lag = 2, n_surr = 20))

    # Specify a non-existent variable name
    expect_error(extract_xi_acf(res_mat, var = "Missing"), "not found")
    expect_error(
        extract_xi_ccf(res_mat, var_x = "X", var_y = "Missing"),
        "not found"
    )

    # Pass an invalid object
    expect_error(extract_xi_acf(list(), var = "X"), "must be of class")
})
