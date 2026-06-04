# R/deprecated.R
# Backward compatibility mappings for deprecated functions

#' @title Deprecated functions in xiacf
#' @description These functions are provided for backward compatibility with
#' older versions of xiacf and will be removed in future releases.
#' @param x A numeric vector or matrix depending on the function.
#' @param max_lag An integer specifying the maximum lag.
#' @param n_surr An integer specifying the number of surrogate datasets.
#' @param sig_level A numeric value specifying the significance or confidence level.
#' @param max_iter An integer specifying the maximum iterations.
#' @param ... Additional arguments passed to the updated functions.
#' @name xiacf-deprecated
NULL

#' @rdname xiacf-deprecated
#' @export
xi_test <- function(
    x,
    max_lag = 10,
    n_surr = 399,
    sig_level = 0.95,
    max_iter = 100,
    ...
) {
    .Deprecated("xi_acf", package = "xiacf")

    # Automatically detect and convert old confidence level (e.g., 0.95) to new significance level (e.g., 0.05)
    if (sig_level > 0.5) {
        warning(
            sprintf(
                "In legacy xi_test(), sig_level = %g is interpreted as a Confidence Level and has been automatically converted to a Significance Level (sig_level = %g) for xi_acf(). Please migrate to xi_acf(sig_level = 0.05) in the future.",
                sig_level,
                1 - sig_level
            ),
            call. = FALSE
        )
        sig_level <- 1 - sig_level
    }

    xi_acf(
        x = x,
        max_lag = max_lag,
        n_surr = n_surr,
        sig_level = sig_level,
        max_iter = max_iter,
        ...
    )
}

#' @rdname xiacf-deprecated
#' @export
generate_iaaft_surrogate <- function(...) {
    .Deprecated("surrogate_iaaft_cpp", package = "xiacf")
    surrogate_iaaft_cpp(...)
}

#' @rdname xiacf-deprecated
#' @export
generate_miaaft_surrogates <- function(...) {
    .Deprecated("surrogate_miaaft_cpp", package = "xiacf")
    surrogate_miaaft_cpp(...)
}

#' @rdname xiacf-deprecated
#' @export
generate_miaaft_surrogate_cpp <- function(...) {
    .Deprecated("surrogate_miaaft_cpp", package = "xiacf")
    surrogate_miaaft_cpp(...)
}

#' @rdname xiacf-deprecated
#' @export
run_rolling_xi_analysis <- function(...) {
    .Deprecated("run_rolling_xi_acf", package = "xiacf")
    run_rolling_xi_acf(...)
}

#' @rdname xiacf-deprecated
#' @export
compute_xi_acf_iaaft <- function(...) {
    .Deprecated("compute_xi_acf_maxstat_cpp", package = "xiacf")
    compute_xi_acf_maxstat_cpp(...)
}

#' @rdname xiacf-deprecated
#' @export
compute_xi_ccf_miaaft <- function(...) {
    .Deprecated("compute_xi_ccf_maxstat_cpp", package = "xiacf")
    compute_xi_ccf_maxstat_cpp(...)
}

#' @rdname xiacf-deprecated
#' @export
compute_xi_matrix_miaaft <- function(...) {
    .Deprecated("compute_xi_matrix_maxstat_cpp", package = "xiacf")
    compute_xi_matrix_maxstat_cpp(...)
}

#' @rdname xiacf-deprecated
#' @export
surrogate_miaaft_cpp <- function(...) {
    .Deprecated(
        msg = "surrogate_miaaft_cpp() is deprecated and its C++ backend has been removed to prevent OOM (Out-Of-Memory) crashes. The package now uses in-place C++ memory management for surrogate generation."
    )
    stop("Execution aborted: Function removed for memory safety.")
}
