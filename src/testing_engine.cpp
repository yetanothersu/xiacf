#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Forward declarations from other files
double xi_coefficient(arma::vec x, arma::vec y);
arma::mat surrogate_iaaft_cpp(const arma::vec& x, int n_surr, int max_iter);
void generate_single_miaaft(const arma::mat& X, arma::mat& X_surr, int max_iter);

//' Compute empirical Xi-ACF and Max-Statistic null distribution
//' @noRd
// [[Rcpp::export]]
List compute_xi_acf_maxstat_cpp(NumericVector x, int max_lag, int n_surr, int max_iter = 100) {
    vec y = as<vec>(x);
    int n = y.n_elem;

    if (n < max_lag + 2) {
        stop("C++ Error: Time series length 'n' must be strictly greater than 'max_lag + 1'.");
    }

    // 1. Original empirical ACF
    vec xi_emp(max_lag, fill::zeros);
    for(int k = 1; k <= max_lag; k++) {
            xi_emp(k-1) = xi_coefficient(y.subvec(0, n - k - 1), y.subvec(k, n - 1));
    }
    
    // 2. Generate Surrogates
    mat surrogates = surrogate_iaaft_cpp(y, n_surr, max_iter);
    
    // 3. Extract Max-statistic
    vec max_stat_vec(n_surr, fill::zeros);
    mat pointwise_xi(max_lag, n_surr, fill::zeros);
    
    for(int s = 0; s < n_surr; s++) {
        vec y_surr = surrogates.col(s);
        double current_max = -1.0;
        
        for(int k = 1; k <= max_lag; k++) {
                double xi_val = xi_coefficient(y_surr.subvec(0, n - k - 1), y_surr.subvec(k, n - 1));
                pointwise_xi(k-1, s) = xi_val;
                if (xi_val > current_max) current_max = xi_val; // Keep track of the maximum
        }
        max_stat_vec(s) = current_max;
    }
    
    return List::create(
        Named("xi_empirical") = xi_emp,
        Named("max_statistic_dist") = max_stat_vec,
        Named("pointwise_dist") = pointwise_xi
    );
}

//' Compute empirical Xi-CCF and Max-Statistic null distribution
//' @noRd
// [[Rcpp::export]]
List compute_xi_ccf_maxstat_cpp(NumericVector x, NumericVector y, int max_lag, int n_surr, int max_iter = 100) {
    vec vx = as<vec>(x);
    vec vy = as<vec>(y);
    int n = vx.n_elem;
    
    if (n < max_lag + 2) {
        stop("C++ Error: Time series length is too short.");
    }
    
    int num_tests = 2 * max_lag + 1; // -max_lag to +max_lag
    vec xi_emp(num_tests, fill::zeros);
    
    // 1. Empirical CCF
    int idx = 0;
    for (int k = -max_lag; k <= max_lag; k++) {
        if (k < 0) { // Y leads X
            xi_emp(idx) = xi_coefficient(vy.subvec(0, n + k - 1), vx.subvec(-k, n - 1));
        } else if (k > 0) { // X leads Y
            xi_emp(idx) = xi_coefficient(vx.subvec(0, n - k - 1), vy.subvec(k, n - 1));
        } else { // lag 0
            xi_emp(idx) = xi_coefficient(vx, vy);
        }
        idx++;
    }
    
    // 2. Generate MIAAFT Surrogates for the bivariate system [X, Y]
    mat X_mat(n, 2);
    X_mat.col(0) = vx;
    X_mat.col(1) = vy;
    mat X_surr(n, 2);
    
    vec max_stat_vec(n_surr, fill::zeros);
    mat pointwise_xi(num_tests, n_surr, fill::zeros);
    
    // 3. FWER Control Loop
    for (int s = 0; s < n_surr; s++) {
        generate_single_miaaft(X_mat, X_surr, max_iter);
        vec sx = X_surr.col(0);
        vec sy = X_surr.col(1);
        
        double current_max = -1.0;
        idx = 0;
        for (int k = -max_lag; k <= max_lag; k++) {
            double val = 0;
            if (k < 0) {
                val = xi_coefficient(sy.subvec(0, n + k - 1), sx.subvec(-k, n - 1));
            } else if (k > 0) {
                val = xi_coefficient(sx.subvec(0, n - k - 1), sy.subvec(k, n - 1));
            } else {
                val = xi_coefficient(sx, sy);
            }
            pointwise_xi(idx, s) = val;
            if (val > current_max) current_max = val; // Find max across all lags
            idx++;
        }
        max_stat_vec(s) = current_max;
    }
    
    return List::create(
        Named("xi_empirical") = xi_emp,
        Named("pointwise_dist") = pointwise_xi,
        Named("max_statistic_dist") = max_stat_vec
    );
}

//' Compute empirical Xi-Matrix and Max-Statistic null distribution
//' @noRd
// [[Rcpp::export]]
List compute_xi_matrix_maxstat_cpp(const arma::mat& X, int max_lag, int n_surr, int max_iter = 100) {
    int n = X.n_rows;
    int p = X.n_cols;
    
    if (n < max_lag + 2) {
        stop("C++ Error: Time series length is too short.");
    }
    
    // Total tests: ACF and CCF for all positive lags up to max_lag
    int num_tests = p * p * max_lag; 
    
    IntegerVector var_lead(num_tests);
    IntegerVector var_lag(num_tests);
    IntegerVector lag_vec(num_tests);
    NumericVector xi_orig(num_tests);
    
    // 1. Empirical Matrix
    int idx = 0;
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < p; j++) {
            for (int k = 1; k <= max_lag; k++) {
                var_lead[idx] = i + 1; // 1-based index for R
                var_lag[idx] = j + 1;
                lag_vec[idx] = k;
                xi_orig[idx] = xi_coefficient(X.col(i).subvec(0, n - k - 1), X.col(j).subvec(k, n - 1));
                idx++;
            }
        }
    }
    
    NumericVector max_stat_vec(n_surr);
    mat X_surr(n, p);
    
    // 2. FWER Control Loop
    for (int s = 0; s < n_surr; s++) {
        generate_single_miaaft(X, X_surr, max_iter);
        double current_max = -1.0;
        
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < p; j++) {
                if (i == j) continue; // Skip autocorrelation for max-statistic
                for (int k = 1; k <= max_lag; k++) {
                    double val = xi_coefficient(X_surr.col(i).subvec(0, n - k - 1), X_surr.col(j).subvec(k, n - 1));
                    if (val > current_max) current_max = val; // Find max across ALL pairs and lags
                }
            }
        }
        max_stat_vec[s] = current_max;
    }
    
    return List::create(
        Named("var_lead") = var_lead,
        Named("var_lag") = var_lag,
        Named("lag") = lag_vec,
        Named("xi_original") = xi_orig,
        Named("max_statistic_dist") = max_stat_vec
    );
}