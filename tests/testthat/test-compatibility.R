# tests/testthat/test-compatibility.R

test_that("sig_level is correctly inverted for backward compatibility", {
    set.seed(42)
    x <- rnorm(100)
    y <- rnorm(100)
    df <- data.frame(A = x, B = y)

    # Specify n_surr = 400 to avoid the "too few surrogates" warning
    # 1. xi_acf check
    expect_warning(
        res_acf <- xiacf::xi_acf(
            x,
            max_lag = 2,
            n_surr = 400,
            sig_level = 0.95
        ),
        "Confidence Level to Significance Level"
    )
    expect_equal(res_acf$sig_level, 0.05)

    # 2. xi_ccf check
    expect_warning(
        res_ccf <- xiacf::xi_ccf(
            x,
            y,
            max_lag = 2,
            n_surr = 400,
            sig_level = 0.95
        ),
        "Confidence Level to Significance Level"
    )
    expect_equal(res_ccf$sig_level, 0.05)

    # 3. xi_matrix check
    expect_warning(
        res_mat <- xiacf::xi_matrix(
            df,
            max_lag = 2,
            n_surr = 400,
            sig_level = 0.95
        ),
        "Confidence Level to Significance Level"
    )
    expect_equal(res_mat$sig_level, 0.05)
})

test_that("deprecated functions issue warnings and function correctly", {
    set.seed(42)
    x <- rnorm(100)

    # xi_test issues two warnings simultaneously: "deprecated" and
    # "Confidence Level" (inverted significance level).
    # Therefore, we use a double expect_warning to capture both.
    expect_warning(
        expect_warning(
            res_test <- xi_test(x, max_lag = 2, n_surr = 400, sig_level = 0.95),
            "Confidence Level"
        ),
        "deprecated"
    )
    expect_s3_class(res_test, "xi_acf")
    expect_equal(res_test$sig_level, 0.05)

    tmp_dir <- tempfile("xi_compat_test")
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE))

    # Remove suppressWarnings to cleanly capture the two warnings here as well
    expect_warning(
        expect_warning(
            res_roll <- run_rolling_xi_analysis(
                x = x,
                window_size = 80,
                step_size = 10,
                max_lag = 2,
                n_surr = 400,
                sig_level = 0.95,
                save_dir = tmp_dir
            ),
            "Confidence Level"
        ),
        "deprecated"
    )
    expect_s3_class(res_roll, "data.frame")
})
