#' Extract Individual Xi-ACF from a Multivariate Xi-Matrix
#' @param object An object of class \code{xi_matrix}.
#' @param var A character string specifying the variable name.
#' @param x_raw Optional. The original data to calculate linear ACF.
#' @return An object of class \code{xi_acf}.
#' @export
extract_xi_acf <- function(object, var, x_raw = NULL) {
    if (!inherits(object, "xi_matrix")) {
        stop("Input must be a 'xi_matrix' object.")
    }

    sub_df <- object$data[
        object$data$Lead_Var == var & object$data$Lag_Var == var,
    ]
    if (nrow(sub_df) == 0) {
        stop(sprintf("Variable '%s' not found.", var))
    }

    # 基本構造の構築
    res_df <- data.frame(
        Lag = sub_df$Lag,
        ACF = NA_real_,
        Xi = sub_df$Xi,
        Pointwise_Threshold = NA_real_,
        Global_Threshold = sub_df$Global_Threshold,
        ACF_CI = stats::qnorm(1 - object$sig_level / 2) / sqrt(object$n)
    )

    # 線形指標の補完
    if (!is.null(x_raw)) {
        x_vec <- as.numeric(x_raw[, var])
        lin_acf <- stats::acf(
            x_vec,
            lag.max = object$max_lag,
            plot = FALSE,
            na.action = stats::na.pass
        )
        res_df$ACF <- as.numeric(lin_acf$acf)[-1]
    }

    res_df$Xi_Excess <- pmax(0, res_df$Xi - res_df$Global_Threshold)

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
#' @param object An object of class \code{xi_matrix}.
#' @param var_x Variable X name.
#' @param var_y Variable Y name.
#' @param x_raw Optional. The original data to calculate linear CCF.
#' @return An object of class \code{xi_ccf}.
#' @export
extract_xi_ccf <- function(object, var_x, var_y, x_raw = NULL) {
    if (!inherits(object, "xi_matrix")) {
        stop("Input must be a 'xi_matrix' object.")
    }

    # X leads Y (Positive lags in CCF sense)
    sub_fwd <- object$data[
        object$data$Lead_Var == var_x & object$data$Lag_Var == var_y,
    ]
    # Y leads X (Negative lags in CCF sense)
    sub_bwd <- object$data[
        object$data$Lead_Var == var_y & object$data$Lag_Var == var_x,
    ]

    if (nrow(sub_fwd) == 0 && nrow(sub_bwd) == 0) {
        stop("Variable pair not found.")
    }

    if (nrow(sub_bwd) > 0) {
        sub_bwd$Lag <- -sub_bwd$Lag
    }
    combined <- rbind(sub_bwd, sub_fwd)
    combined <- combined[order(combined$Lag), ]

    res_df <- data.frame(
        Lag = combined$Lag,
        CCF = NA_real_,
        Xi = combined$Xi,
        Pointwise_Threshold = NA_real_,
        Global_Threshold = combined$Global_Threshold,
        CCF_CI = stats::qnorm(1 - object$sig_level / 2) / sqrt(object$n)
    )

    if (!is.null(x_raw)) {
        lin_ccf <- stats::ccf(
            as.numeric(x_raw[, var_x]),
            as.numeric(x_raw[, var_y]),
            lag.max = object$max_lag,
            plot = FALSE,
            na.action = stats::na.pass
        )
        # stats::ccfの結果とラグを正確に紐付け
        ccf_map <- stats::setNames(
            as.numeric(lin_ccf$acf),
            as.numeric(lin_ccf$lag)
        )
        res_df$CCF <- unname(ccf_map[as.character(res_df$Lag)])
    }

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
