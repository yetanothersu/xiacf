#' Extract Individual Xi-ACF from Xi-Matrix
#'
#' Extracts the autocorrelation results for a specific variable from a \code{xi_matrix} object
#' and converts it into a \code{xi_acf} S3 object.
#'
#' @param object An object of class \code{xi_matrix}.
#' @param var A character string specifying the variable name to extract.
#' @param x_raw Optional. The original multivariate data (matrix or data.frame) used to
#'     compute the matrix. If provided, standard linear ACF will be re-calculated.
#'
#' @return An object of class \code{xi_acf}.
#' @export
extract_xi_acf <- function(object, var, x_raw = NULL) {
    if (!inherits(object, "xi_matrix")) {
        stop("Input 'object' must be a 'xi_matrix' object.")
    }

    sub_df <- object$data[
        object$data$Lead_Var == var & object$data$Lag_Var == var,
    ]

    if (nrow(sub_df) == 0) {
        stop(sprintf("Variable '%s' not found in the matrix object.", var))
    }

    res_df <- data.frame(
        Lag = sub_df$Lag,
        ACF = NA_real_,
        Xi = sub_df$Xi,
        Xi_Threshold = sub_df$Xi_Threshold,
        ACF_CI = NA_real_,
        Xi_Excess = sub_df$Xi_Excess
    )

    if (!is.null(x_raw)) {
        x_vec <- as.numeric(x_raw[, var])
        acf_res <- stats::acf(
            x_vec,
            lag.max = object$max_lag,
            plot = FALSE,
            na.action = stats::na.pass
        )
        # Safely map ACF values using Lag as the key
        acf_map <- stats::setNames(
            as.numeric(acf_res$acf),
            as.numeric(acf_res$lag)
        )
        res_df$ACF <- unname(acf_map[as.character(res_df$Lag)])
        res_df$ACF_CI <- stats::qnorm((1 + object$sig_level) / 2) /
            sqrt(object$n)
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

#' Extract Individual Xi-CCF from Xi-Matrix
#'
#' Extracts the cross-correlation results for a specific pair of variables from
#' a \code{xi_matrix} object and converts it into a \code{xi_ccf} S3 object.
#'
#' @param object An object of class \code{xi_matrix}.
#' @param var_x A character string for the first variable (X).
#' @param var_y A character string for the second variable (Y).
#' @param x_raw Optional. The original multivariate data used to compute the matrix.
#'
#' @return An object of class \code{xi_ccf}.
#' @export
extract_xi_ccf <- function(object, var_x, var_y, x_raw = NULL) {
    if (!inherits(object, "xi_matrix")) {
        stop("Input 'object' must be a 'xi_matrix' object.")
    }

    df_fwd <- object$data[
        object$data$Lead_Var == var_x & object$data$Lag_Var == var_y,
    ]
    df_bwd <- object$data[
        object$data$Lead_Var == var_y & object$data$Lag_Var == var_x,
    ]

    if (nrow(df_fwd) == 0 || nrow(df_bwd) == 0) {
        stop(sprintf(
            "Variable pair '%s' and '%s' not found in the matrix object.",
            var_x,
            var_y
        ))
    }

    res_fwd <- data.frame(
        Direction = "X leads Y",
        Lag = df_fwd$Lag,
        CCF = NA_real_,
        Xi = df_fwd$Xi,
        Xi_Threshold = df_fwd$Xi_Threshold,
        CCF_CI = NA_real_,
        Xi_Excess = df_fwd$Xi_Excess
    )

    res_bwd <- data.frame(
        Direction = "Y leads X",
        Lag = df_bwd$Lag,
        CCF = NA_real_,
        Xi = df_bwd$Xi,
        Xi_Threshold = df_bwd$Xi_Threshold,
        CCF_CI = NA_real_,
        Xi_Excess = df_bwd$Xi_Excess
    )

    if (!is.null(x_raw)) {
        x_vec <- as.numeric(x_raw[, var_x])
        y_vec <- as.numeric(x_raw[, var_y])
        ccf_ci <- stats::qnorm((1 + object$sig_level) / 2) / sqrt(object$n)

        # Forward CCF
        ccf_fwd_full <- stats::ccf(
            y_vec,
            x_vec,
            lag.max = object$max_lag,
            plot = FALSE,
            na.action = stats::na.pass
        )
        ccf_fwd_map <- stats::setNames(
            as.numeric(ccf_fwd_full$acf),
            as.numeric(ccf_fwd_full$lag)
        )
        res_fwd$CCF <- unname(ccf_fwd_map[as.character(res_fwd$Lag)])
        res_fwd$CCF_CI <- ccf_ci

        # Backward CCF
        ccf_bwd_full <- stats::ccf(
            x_vec,
            y_vec,
            lag.max = object$max_lag,
            plot = FALSE,
            na.action = stats::na.pass
        )
        ccf_bwd_map <- stats::setNames(
            as.numeric(ccf_bwd_full$acf),
            as.numeric(ccf_bwd_full$lag)
        )
        res_bwd$CCF <- unname(ccf_bwd_map[as.character(res_bwd$Lag)])
        res_bwd$CCF_CI <- ccf_ci
    }

    structure(
        list(
            data = rbind(res_fwd, res_bwd),
            n = object$n,
            max_lag = object$max_lag,
            n_surr = object$n_surr,
            bidirectional = TRUE,
            sig_level = object$sig_level
        ),
        class = "xi_ccf"
    )
}
