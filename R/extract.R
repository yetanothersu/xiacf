#' Extract Univariate Xi-ACF from a Multivariate Xi-Matrix
#'
#' @param obj An object of class \code{xi_matrix}.
#' @param var A character string specifying the variable name to extract.
#' @param ... Additional arguments passed to xi_acf.
#'
#' @return An object of class \code{xi_acf}.
#' @export
extract_xi_acf <- function(obj, var, ...) {
    if (!inherits(obj, "xi_matrix")) {
        stop("Object must be of class 'xi_matrix'.")
    }

    raw_data_name <- "data_raw"
    if (is.null(obj[[raw_data_name]])) {
        stop(
            "Raw data is not preserved in the xi_matrix object. Cannot compute on-demand ACF."
        )
    }

    if (!(var %in% colnames(obj[[raw_data_name]]))) {
        stop(sprintf(
            "Variable '%s' not found in the original matrix data.",
            var
        ))
    }

    x_raw <- obj[[raw_data_name]][[var]]

    # Re-calculate xi_acf on demand
    res <- xi_acf(
        x = x_raw,
        max_lag = obj$max_lag,
        n_surr = obj$n_surr,
        sig_level = obj$sig_level,
        ...
    )

    # Clean up the messy deparsed variable name
    res$x_name <- var

    return(res)
}

#' Extract Bivariate Xi-CCF from a Multivariate Xi-Matrix
#'
#' @param obj An object of class \code{xi_matrix}.
#' @param var_x A character string specifying the lead variable.
#' @param var_y A character string specifying the lag variable.
#' @param ... Additional arguments passed to xi_ccf.
#'
#' @return An object of class \code{xi_ccf}.
#' @export
extract_xi_ccf <- function(obj, var_x, var_y, ...) {
    if (!inherits(obj, "xi_matrix")) {
        stop("Object must be of class 'xi_matrix'.")
    }

    raw_data_name <- "data_raw"
    if (is.null(obj[[raw_data_name]])) {
        stop(
            "Raw data is not preserved in the xi_matrix object. Cannot compute on-demand CCF."
        )
    }

    if (
        !(var_x %in% colnames(obj[[raw_data_name]])) ||
            !(var_y %in% colnames(obj[[raw_data_name]]))
    ) {
        stop(sprintf(
            "Variables '%s' and/or '%s' not found in the original matrix data.",
            var_x,
            var_y
        ))
    }

    # Re-calculate xi_ccf on demand
    res <- xi_ccf(
        x = obj[[raw_data_name]][[var_x]],
        y = obj[[raw_data_name]][[var_y]],
        max_lag = obj$max_lag,
        n_surr = obj$n_surr,
        sig_level = obj$sig_level,
        ...
    )

    # Capture the messy auto-generated names from substitute()
    old_x <- res$x_name
    old_y <- res$y_name

    # Override the metadata names
    res$x_name <- var_x
    res$y_name <- var_y

    # Clean up the data frame columns so autoplot facet labels are elegant
    if (!is.null(res$data$Lead_Var)) {
        res$data$Lead_Var[res$data$Lead_Var == old_x] <- var_x
        res$data$Lead_Var[res$data$Lead_Var == old_y] <- var_y
    }
    if (!is.null(res$data$Lag_Var)) {
        res$data$Lag_Var[res$data$Lag_Var == old_x] <- var_x
        res$data$Lag_Var[res$data$Lag_Var == old_y] <- var_y
    }

    return(res)
}
