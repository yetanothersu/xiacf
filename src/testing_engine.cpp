#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Forward declarations from other files
double xi_coefficient(arma::vec x, arma::vec y);
void generate_single_iaaft_worker(const arma::vec& x_sorted, const arma::vec& X_amp, arma::vec& x_surr, int max_iter);
arma::mat surrogate_iaaft_cpp(const arma::vec& x, int n_surr, int max_iter);
void generate_single_miaaft(const arma::mat& X, arma::mat& X_surr, int max_iter);

//' Compute empirical Xi-ACF and Max-Statistic null distribution
//' @export
//' @keywords internal
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
    vec max_stat_vec(n_surr, fill::zeros);
    mat pointwise_xi(max_lag, n_surr, fill::zeros);
        
    // Precompute for IAAFT (Done only ONCE per time series)
    vec y_sorted = sort(y);
    cx_vec Y_f = fft(y);
    vec Y_amp = abs(Y_f);
    
    // Allocate a SINGLE vector to hold the surrogate (Memory Churn Eliminated)
    vec surr(n);
    
    for(int s = 0; s < n_surr; s++) {
        // Generate surrogate in-place
        generate_single_iaaft_worker(y_sorted, Y_amp, surr, max_iter);
        
        double max_val = -1.0;
        for(int k = 1; k <= max_lag; k++) {
            double xi_val = xi_coefficient(surr.subvec(0, n - k - 1), surr.subvec(k, n - 1));
            pointwise_xi(k-1, s) = xi_val;
            if(xi_val > max_val) max_val = xi_val;
        }
        max_stat_vec(s) = max_val;
    }
    
    return List::create(
        Named("xi_empirical") = xi_emp,
        Named("max_statistic_dist") = max_stat_vec,
        Named("pointwise_dist") = pointwise_xi
    );
}

//' Compute empirical Xi-CCF and Max-Statistic null distribution
//' @export
//' @keywords internal
// [[Rcpp::export]]
List compute_xi_ccf_maxstat_cpp(NumericVector x, NumericVector y, int max_lag, int n_surr, int max_iter = 100, bool both_directions = true) {
    vec vx = as<vec>(x);
    vec vy = as<vec>(y);
    int n = vx.n_elem;

    if (n < max_lag + 2) {
        stop("C++ Error: Time series length is too short.");
    }

    // 1. Empirical CCF (Compute directed dependencies starting from lag 0)
    vec xi_emp_x_leads(max_lag + 1, fill::zeros);
    vec xi_emp_y_leads;
    
    if (both_directions) {
        xi_emp_y_leads.zeros(max_lag + 1);
    }

    for (int k = 0; k <= max_lag; k++) {
        // X leads Y: X_t -> Y_{t+k} (k=0 is contemporaneous X -> Y)
        xi_emp_x_leads(k) = xi_coefficient(vx.subvec(0, n - k - 1), vy.subvec(k, n - 1));
        
        if (both_directions) {
            xi_emp_y_leads(k) = xi_coefficient(vy.subvec(0, n - k - 1), vx.subvec(k, n - 1));
        }
    }

    // 2. Surrogate distributions (Dual-Family FWER Control)
    arma::mat X(n, 2);
    X.col(0) = vx;
    X.col(1) = vy;

    arma::vec max_dist_lag0(n_surr, fill::zeros);
    arma::vec max_dist_lagged(n_surr, fill::zeros);

    for (int s = 0; s < n_surr; s++) {
        arma::mat X_surr(n, 2);
        generate_single_miaaft(X, X_surr, max_iter);
        vec x_surr = X_surr.col(0);
        vec y_surr = X_surr.col(1);

        double current_max_lag0 = -1.0;
        double current_max_lagged = -1.0;

        // --- Family A: Lag 0 (Evaluate both directions due to functional asymmetry) ---
        double val_0_xy = xi_coefficient(x_surr, y_surr);
        current_max_lag0 = val_0_xy; 
        
        if (both_directions) {
            double val_0_yx = xi_coefficient(y_surr, x_surr);
            if (val_0_yx > current_max_lag0) current_max_lag0 = val_0_yx;
        }

        // --- Family B: Lag > 0 (Temporal Propagation) ---
        for (int k = 1; k <= max_lag; k++) {
            double val_x_leads = xi_coefficient(x_surr.subvec(0, n - k - 1), y_surr.subvec(k, n - 1));
            if (val_x_leads > current_max_lagged) current_max_lagged = val_x_leads;

            if (both_directions) {
                double val_y_leads = xi_coefficient(y_surr.subvec(0, n - k - 1), x_surr.subvec(k, n - 1));
                if (val_y_leads > current_max_lagged) current_max_lagged = val_y_leads;
            }
        }
        max_dist_lag0(s) = current_max_lag0;
        max_dist_lagged(s) = current_max_lagged;
    }

    return List::create(
        Named("xi_emp_x_leads") = xi_emp_x_leads,
        Named("xi_emp_y_leads") = xi_emp_y_leads,
        Named("max_dist_lag0") = max_dist_lag0,
        Named("max_dist_lagged") = max_dist_lagged
    );
}

//' Compute empirical Xi-Matrix and Max-Statistic null distribution
//' @export
//' @keywords internal
// [[Rcpp::export]]
List compute_xi_matrix_maxstat_cpp(const arma::mat& X, int max_lag, int n_surr, int max_iter = 100) {
    int n = X.n_rows;
    int p = X.n_cols;

    if (n < max_lag + 2) {
        stop("C++ Error: Time series length is too short.");
    }

    // 1. Empirical Xi-Matrix
    arma::cube xi_emp(p, p, max_lag + 1, fill::zeros);
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < p; j++) {
            for (int k = 0; k <= max_lag; k++) {
                xi_emp(i, j, k) = xi_coefficient(X.col(i).subvec(0, n - k - 1), X.col(j).subvec(k, n - 1));
            }
        }
    }

    // 2. Surrogate distributions (Dual-Family FWER Control)
    arma::vec max_dist_lag0(n_surr, fill::zeros);
    arma::vec max_dist_lagged(n_surr, fill::zeros);

    for (int s = 0; s < n_surr; s++) {
        arma::mat X_surr(n, p);
        generate_single_miaaft(X, X_surr, max_iter);
        
        double current_max_lag0 = -1.0;
        double current_max_lagged = -1.0;
        
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < p; j++) {
                if (i == j) continue; // Completely exclude autocorrelation from the FWER search space
                
                // --- Family A: Lag 0 (Contemporaneous Confounding) ---
                double val_0 = xi_coefficient(X_surr.col(i), X_surr.col(j));
                if (val_0 > current_max_lag0) current_max_lag0 = val_0;
                
                // --- Family B: Lag > 0 (Temporal Causal Propagation) ---
                // * Start the loop from k = 1 to separate temporal effects
                for (int k = 1; k <= max_lag; k++) { 
                    double val_k = xi_coefficient(X_surr.col(i).subvec(0, n - k - 1), X_surr.col(j).subvec(k, n - 1));
                    if (val_k > current_max_lagged) current_max_lagged = val_k;
                }
            }
        }
        max_dist_lag0(s) = current_max_lag0;
        max_dist_lagged(s) = current_max_lagged;
    }

    // Return two independent null distributions to R
    return List::create(
        Named("xi_emp") = xi_emp,
        Named("max_dist_lag0") = max_dist_lag0,
        Named("max_dist_lagged") = max_dist_lagged
    );
}