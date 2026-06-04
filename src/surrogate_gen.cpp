#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' @noRd
void generate_single_iaaft_worker(const arma::vec& x_sorted, const arma::vec& X_amp, arma::vec& x_surr, int max_iter) {
    int n = x_sorted.n_elem;
    
    vec rand_phases(n);
    for(int k = 0; k < n; ++k) {
        rand_phases[k] = R::runif(0, 1) * 2.0 * M_PI;
    }
    
    cx_vec S_f = X_amp % exp(cx_vec(zeros<vec>(n), rand_phases));
    vec s_t = real(ifft(S_f));
    
    for (int iter = 0; iter < max_iter; ++iter) {
        uvec rank_idx = sort_index(s_t);
        vec s_t_matched(n);
        s_t_matched(rank_idx) = x_sorted;
        
        cx_vec S_matched_f = fft(s_t_matched);
        vec phases = arg(S_matched_f);
        cx_vec S_f_new = X_amp % exp(cx_vec(zeros<vec>(n), phases));
        vec s_t_new = real(ifft(S_f_new));
        
        if (abs(s_t_new - s_t).max() < 1e-6) {
            s_t = s_t_new;
            break;
        }
        s_t = s_t_new;
    }
    
    uvec rank_idx = sort_index(s_t);
    x_surr(rank_idx) = x_sorted;
}

//' Generate Multiple IAAFT Surrogates (Univariate)
//' @param x A numeric vector.
//' @param n_surr Number of surrogates to generate.
//' @param max_iter Maximum iterations for IAAFT.
//' @return A matrix of surrogates (N x n_surr).
//' @export
// [[Rcpp::export]]
arma::mat surrogate_iaaft_cpp(const arma::vec& x, int n_surr, int max_iter = 100) {
    int n = x.n_elem;
    arma::mat surrogates(n, n_surr);
    
    // Precompute FFT and sort outside the loop for extreme performance
    vec x_sorted = sort(x);
    cx_vec X_f = fft(x);
    vec X_amp = abs(X_f);
    
    vec x_surr(n);
    for (int s = 0; s < n_surr; ++s) {
        generate_single_iaaft_worker(x_sorted, X_amp, x_surr, max_iter);
        surrogates.col(s) = x_surr;
    }
    
    return surrogates;
}

// Helper: Generate shared phases for MIAAFT
vec generate_shared_random_phases(int n) {
    vec phases(n, fill::zeros);
    for(int i = 1; i <= (n - 1) / 2; i++) {
        double p = R::runif(0, 2 * M_PI);
        phases(i) = p;
        phases(n - i) = -p; 
    }
    if(n % 2 == 0) phases(n / 2) = 0;
    return phases;
}

// Helper: Single MIAAFT generation
void generate_single_miaaft(const arma::mat& X, arma::mat& X_surr, int max_iter) {
    int n = X.n_rows;
    int p = X.n_cols;
    
    mat X_sorted(n, p);
    cx_mat X_f(n, p);
    mat X_amp(n, p);
    
    for (int j = 0; j < p; ++j) {
        X_sorted.col(j) = sort(X.col(j));
        X_f.col(j) = fft(X.col(j));
        X_amp.col(j) = abs(X_f.col(j));
    }
    
    vec shared_phases = generate_shared_random_phases(n);
    mat S_t(n, p);
    
    for (int j = 0; j < p; ++j) {
        cx_vec S_f = X_amp.col(j) % exp(cx_vec(zeros<vec>(n), shared_phases));
        S_t.col(j) = real(ifft(S_f));
    }
    
    for (int iter = 0; iter < max_iter; ++iter) {
        mat S_t_matched(n, p);
        for (int j = 0; j < p; ++j) {
            uvec rank_idx = sort_index(S_t.col(j));
            vec temp(n);                           
            temp(rank_idx) = X_sorted.col(j);      
            S_t_matched.col(j) = temp;             
        }
        
        mat S_t_new(n, p);
        for (int j = 0; j < p; ++j) {
            cx_vec S_matched_f = fft(S_t_matched.col(j));
            vec phases = arg(S_matched_f);
            cx_vec S_f = X_amp.col(j) % exp(cx_vec(zeros<vec>(n), phases));
            S_t_new.col(j) = real(ifft(S_f));
        }
        
        if (abs(S_t_new - S_t).max() < 1e-6) {
            S_t = S_t_new;
            break;
        }
        S_t = S_t_new;
    }
    
    for (int j = 0; j < p; ++j) {
        uvec rank_idx = sort_index(S_t.col(j));
        vec temp(n);
        temp(rank_idx) = X_sorted.col(j);
        X_surr.col(j) = temp;
    }
}

//' Generate Multiple MIAAFT Surrogates (3D Array / Cube)
//' @param X A numeric matrix (N x p).
//' @param n_surr Number of surrogates to generate.
//' @param max_iter Maximum iterations for MIAAFT.
//' @return A 3D array (arma::cube) of dimensions N x p x n_surr.
//' @export
// [[Rcpp::export]]
arma::cube surrogate_miaaft_cpp(const arma::mat& X, int n_surr, int max_iter = 100) {
    int n = X.n_rows;
    int p = X.n_cols;
    arma::cube surr_cube(n, p, n_surr);
    
    for (int b = 0; b < n_surr; ++b) {
        arma::mat X_surr(n, p);
        generate_single_miaaft(X, X_surr, max_iter);
        surr_cube.slice(b) = X_surr;
    }
    return surr_cube;
}