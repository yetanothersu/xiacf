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
    
    // 1. Original empirical ACF
    vec xi_emp(max_lag, fill::zeros);
    for(int k = 1; k <= max_lag; k++) {
        if(n > k) {
            xi_emp(k-1) = xi_coefficient(y.subvec(0, n - k - 1), y.subvec(k, n - 1));
        }
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
            if(n > k) {
                double xi_val = xi_coefficient(y_surr.subvec(0, n - k - 1), y_surr.subvec(k, n - 1));
                pointwise_xi(k-1, s) = xi_val;
                if (xi_val > current_max) current_max = xi_val; // Keep track of the maximum
            }
        }
        max_stat_vec(s) = current_max;
    }
    
    return List::create(
        Named("xi_empirical") = xi_emp,
        Named("max_statistic_dist") = max_stat_vec,
        Named("pointwise_dist") = pointwise_xi
    );
}

// Helper: Compute full p x p xi_matrix
arma::mat compute_xi_matrix_internal(const arma::mat& X) {
    int p = X.n_cols;
    arma::mat xi_mat(p, p, fill::zeros);
    for (int i = 0; i < p; ++i) {
        for (int j = 0; j < p; ++j) {
            if (i == j) {
                xi_mat(i, j) = 1.0;
            } else {
                xi_mat(i, j) = xi_coefficient(X.col(i), X.col(j));
            }
        }
    }
    return xi_mat;
}

//' Compute empirical Xi-Matrix and Max-Statistic null distribution
//' @noRd
// [[Rcpp::export]]
List compute_xi_matrix_maxstat_cpp(const arma::mat& X, int n_surr, int max_iter = 100) {
    int n = X.n_rows;
    int p = X.n_cols;
    
    arma::mat xi_emp = compute_xi_matrix_internal(X);
    arma::vec max_stat_vec(n_surr);
    
    arma::mat X_surr(n, p);
    arma::mat xi_surr(p, p);
    
    for (int b = 0; b < n_surr; ++b) {
        // Generate surrogate in-place to save memory
        generate_single_miaaft(X, X_surr, max_iter);
        xi_surr = compute_xi_matrix_internal(X_surr);
        
        double current_max = -1.0;
        for (int i = 0; i < p; ++i) {
            for (int j = 0; j < p; ++j) {
                if (i != j) {
                    if (xi_surr(i, j) > current_max) {
                        current_max = xi_surr(i, j); // Extract global max
                    }
                }
            }
        }
        max_stat_vec(b) = current_max;
    }
    
    return List::create(
        Named("xi_empirical") = xi_emp,
        Named("max_statistic_dist") = max_stat_vec
    );
}