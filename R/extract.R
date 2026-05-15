#' Extract Individual Xi-ACF from a Multivariate Xi-Matrix
#'
#' @param object An object of class \code{xi_matrix}.
#' @param var A character string specifying the variable name to extract.
#' @param x_raw Optional. The original multivariate data matrix/data.frame.
#'   If provided, standard linear ACF will be calculated and included.
#'
#' @return An object of class \code{xi_acf}.
#' @export
extract_xi_acf <- function(object, var, x_raw = NULL) {
    if (!inherits(object, "xi_matrix")) {
        stop("Input must be a 'xi_matrix' object.")
    }

    # Extract diagonal components (autocorrelation)
    sub_df <- object$data[
        object$data$Lead_Var == var & object$data$Lag_Var == var,
    ]
    if (nrow(sub_df) == 0) {
        stop(sprintf("Variable '%s' not found in the matrix.", var))
    }

    sub_df <- sub_df[order(sub_df$Lag), ]

    # Construct the exact same data structure as xi_acf
    res_df <- data.frame(
        Lag = sub_df$Lag,
        Xi = sub_df$Xi,
        Pointwise_Threshold = NA_real_, # NA because it cannot be obtained during matrix calculation
        Global_Threshold = sub_df$Global_Threshold,
        ACF_CI = NA_real_
    )
    res_df$Xi_Excess <- pmax(0, res_df$Xi - res_df$Global_Threshold)

    # Calculate standard linear ACF as well if raw data is provided
    if (!is.null(x_raw)) {
        x_vec <- as.numeric(x_raw[, var])
        acf_res <- stats::acf(
            x_vec,
            lag.max = object$max_lag,
            plot = FALSE,
            na.action = stats::na.fail
        )
        res_df$ACF_CI <- stats::qnorm(1 - object$sig_level / 2) / sqrt(object$n)
        # You may add res_df$ACF column if necessary
    }

    structure(
        list(
            data = res_df,
            n = object$n,
            max_lag = object$max_lag,
            n_surr = object$n_surr,
            sig_level = object$sig_level
        ),
        class = "xi_acf"
    )
}

#' Extract Pairwise Xi-CCF from a Multivariate Xi-Matrix
#'
#' @param object An object of class \code{xi_matrix}.
#' @param var_x A character string specifying the first variable (X).
#' @param var_y A character string specifying the second variable (Y).
#' @param x_raw Optional. The original multivariate data matrix/data.frame.
#'
#' @return An object of class \code{xi_ccf}.
#' @export
extract_xi_ccf <- function(object, var_x, var_y, x_raw = NULL) {
    if (!inherits(object, "xi_matrix")) {
        stop("Input must be a 'xi_matrix' object.")
    }

    # Direction where X leads (Lag > 0)
    sub_fwd <- object$data[
        object$data$Lead_Var == var_x & object$data$Lag_Var == var_y,
    ]
    # Direction where Y leads (X is lag, Lag < 0)
    sub_bwd <- object$data[
        object$data$Lead_Var == var_y & object$data$Lag_Var == var_x,
    ]

    if (nrow(sub_fwd) == 0 && nrow(sub_bwd) == 0) {
        stop("Specified variable pair not found in the matrix.")
    }

    # Invert the sign of Lag and combine (to match xi_ccf structure)
    if (nrow(sub_bwd) > 0) {
        sub_bwd$Lag <- -sub_bwd$Lag
    }

    combined <- rbind(sub_bwd, sub_fwd)
    combined <- combined[order(combined$Lag), ]

    res_df <- data.frame(
        Lag = combined$Lag,
        Xi = combined$Xi,
        Pointwise_Threshold = NA_real_,
        Global_Threshold = combined$Global_Threshold
    )
    res_df$Xi_Excess <- pmax(0, res_df$Xi - res_df$Global_Threshold)

    structure(
        list(
            data = res_df,
            n = object$n,
            max_lag = object$max_lag,
            n_surr = object$n_surr,
            sig_level = object$sig_level
        ),
        class = "xi_ccf"
    )
}
