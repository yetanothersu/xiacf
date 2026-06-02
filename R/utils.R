#' Validate Number of Surrogates for Max-Statistic FWER Control
#'
#' This is an internal helper function to ensure that the user has provided
#' a sufficient number of surrogate datasets to estimate the specified
#' significance level (e.g., 0.05) stably in the tail of the empirical distribution.
#'
#' @param n_surr An integer specifying the number of surrogates requested.
#' @param sig_level A numeric value for the significance level (e.g., 0.05).
#' @param num_tests An integer specifying the number of simultaneous tests (e.g., max_lag or p*(p-1)/2).
#' @return NULL. Throws an error or warning if n_surr is insufficient.
#' @noRd
check_surrogate_count <- function(n_surr, sig_level, num_tests) {
    # 1. Minimum required just to calculate the quantile
    # Example: If sig_level = 0.05, we need at least 1/0.05 - 1 = 19 surrogates.
    min_required <- ceiling(1 / sig_level) - 1
    if (n_surr < min_required) {
        stop(sprintf(
            "n_surr = %d is too small to calculate the %g%% threshold. Minimum required is %d.",
            n_surr,
            (1 - sig_level) * 100,
            min_required
        ))
    }

    # 2. Recommended number for stable max-statistic estimation in the right tail
    # The tail becomes thicker with more simultaneous tests.
    recommended <- ceiling(num_tests / sig_level)
    recommended <- max(399, recommended) # Enforce typical non-parametric baseline (e.g., 399 or 999)

    if (n_surr < recommended) {
        warning(sprintf(
            "For %d simultaneous tests at sig_level = %g, the empirical distribution of the max-statistic may be unstable with n_surr = %d. Recommended n_surr is at least %d.",
            num_tests,
            sig_level,
            n_surr,
            recommended
        ))
    }

    invisible(NULL)
}
