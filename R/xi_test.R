#' Xi-ACF Test for Time Series
#'
#' Calculates Chatterjee's Xi and the standard Autocorrelation Function (ACF)
#' along with their respective significance thresholds.
#'
#' @useDynLib xiacf, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats na.omit sd acf quantile qnorm
#'
#' @param x A numeric vector representing the time series data.
#' @param max_lag An integer specifying the maximum number of lags to compute.
#' @param n_surr An integer specifying the number of surrogate datasets to generate for the IAAFT test.
#' @return An object of class \code{"xi_test"} containing the computed statistics and metadata.
#' @export
xi_test <- function(x, max_lag = 20, n_surr = 100) {
    # --- 1. Robust Input Validation ---

    # Check if input is a numeric vector
    if (!is.numeric(x)) {
        stop("Input 'x' must be a numeric vector.")
    }

    # Check for NA values: remove them and issue a warning
    if (any(is.na(x))) {
        warning("'x' contains NA values. Removing them before analysis.")
        x <- stats::na.omit(x)
        # Remove attributes added by na.omit to keep it as a pure numeric vector
        x <- as.numeric(x)
    }

    n <- length(x)

    # Check data length
    if (n < 5) {
        stop("Time series length is too short (n < 5).")
    }

    # Check for constant vector (zero variance)
    if (stats::sd(x) == 0) {
        stop(
            "Input 'x' is constant (zero variance). Correlation cannot be computed."
        )
    }

    # Check lag length against series length
    if (max_lag >= n) {
        warning("max_lag is >= series length. Reducing max_lag to n - 1.")
        max_lag <- n - 1
    }

    # --- 2. Calculate Standard ACF (Linear) ---
    # Compute standard ACF without plotting
    acf_res <- stats::acf(x, lag.max = max_lag, plot = FALSE)

    # Exclude lag 0 (which is always identically 1.0)
    acf_vals <- as.numeric(acf_res$acf)[-1]

    # Calculate ACF Confidence Interval (95% i.i.d. assumption: +/- 1.96/sqrt(n))
    acf_ci <- stats::qnorm(0.975) / sqrt(n)

    # --- 3. Calculate Xi (Non-linear) using C++ backend ---
    xi_res <- compute_xi_lags(x, max_lag, n_surr)

    # --- 4. Calculate Threshold from Surrogates ---
    # Handle the case where n_surr = 0
    xi_threshold <- rep(NA, max_lag)

    if (n_surr > 0) {
        # Calculate the 95th percentile threshold for each lag
        # na.rm = TRUE is required as uncomputable lags are filled with NaN in the C++ backend
        xi_threshold <- apply(xi_res$xi_surrogates, 1, function(row) {
            stats::quantile(row, 0.95, na.rm = TRUE)
        })
    }

    # --- 5. Construct Result Object ---
    # Compile results into a data frame for easier plotting with ggplot2
    df_res <- data.frame(
        Lag = 1:max_lag,
        ACF = acf_vals[1:max_lag], # Safeguard in case ACF results are shorter than max_lag
        Xi = as.numeric(xi_res$xi_original),
        Xi_Threshold_95 = xi_threshold,
        ACF_CI = acf_ci # Store the constant ACF CI as a column for easy plotting
    )

    # Define and return the S3 object
    structure(
        list(
            data = df_res,
            n = n,
            max_lag = max_lag,
            n_surr = n_surr
        ),
        class = "xi_test"
    )
}

#' Print method for xi_test objects
#'
#' Prints a concise summary of the \code{xi_test} object, including metadata
#' and the first few rows of the computed statistics.
#'
#' @param x An object of class \code{"xi_test"}.
#' @param ... Additional arguments passed to other methods (currently ignored).
#' @return Invisibly returns the original object.
#' @export
print.xi_test <- function(x, ...) {
    cat("\n\tChatterjee's Xi-ACF Test\n\n")
    cat("Data length:  ", x$n, "\n")
    cat("Max lag:      ", x$max_lag, "\n")
    cat("Surrogates:   ", x$n_surr, " (IAAFT)\n\n")

    # Select specific columns to display
    # Exclude the constant ACF_CI column to focus on the main metrics
    show_cols <- c("Lag", "ACF", "Xi", "Xi_Threshold_95")

    # Safely intersect with existing column names
    cols_exist <- intersect(show_cols, names(x$data))

    # Display only the first 5 rows
    print(utils::head(x$data[, cols_exist], 5), row.names = FALSE)

    # Display a message if there are remaining rows
    if (nrow(x$data) > 5) {
        cat("... with", nrow(x$data) - 5, "more lags\n")
    }

    invisible(x)
}
