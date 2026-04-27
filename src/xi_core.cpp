#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// ============================================================
// Helper 1: R-synchronized Random Shuffle Function
// ============================================================
// Shuffles indices using R's random number generator to ensure
// exact reproducibility with set.seed() from the R environment.
void shuffle_indices_with_r_seed(uvec& idx) {
    int n = idx.n_elem;
    for (int i = n - 1; i > 0; i--) {
        int j = floor(R::runif(0, 1) * (i + 1));
        uword tmp = idx[i];
        idx[i] = idx[j];
        idx[j] = tmp;
    }
}

// ============================================================
// Helper 2: Rank Calculation with Random Tie-Breaking
// ============================================================
// Calculates the rank of elements while breaking ties randomly
// to ensure Chatterjee's Xi coefficient remains well-defined.
uvec rank_random_ties_r_sync(vec x) {
    uword n = x.n_elem;
    uvec idx = sort_index(x);
    vec x_sorted = x(idx);
    
    for (uword i = 0; i < n; ) {
        uword j = i + 1;
        // Detect ties
        while (j < n && x_sorted[j] == x_sorted[i]) {
            j++;
        }
        // If ties are found, shuffle their indices randomly
        if (j > i + 1) {
            uvec sub_idx = idx.subvec(i, j - 1);
            shuffle_indices_with_r_seed(sub_idx);
            idx.subvec(i, j - 1) = sub_idx;
        }
        i = j;
    }
    
    uvec ranks(n);
    for(uword k = 0; k < n; k++) {
        ranks(idx(k)) = k; 
    }
    return ranks;
}

// ============================================================
// Main Logic: Chatterjee's Xi Coefficient
// ============================================================
//' Calculate Chatterjee's Rank Correlation Coefficient (Xi)
//' 
//' Computes Chatterjee's rank correlation coefficient (Xi) between two numeric vectors.
//' Ties are broken uniformly at random to ensure strict inequalities.
//' 
//' @param x A numeric vector.
//' @param y A numeric vector of the same length as \code{x}.
//' @return A numeric scalar representing the Chatterjee's Xi coefficient.
//' @export
// [[Rcpp::export]]
double xi_coefficient(arma::vec x, arma::vec y) {
    uword n = x.n_elem;
    if (n < 2) return 0.0;

    uvec PI = sort_index(x); 
    vec y_sorted = y(PI);
    uvec r = rank_random_ties_r_sync(y_sorted);
    
    double term1 = 0.0;
    for (uword i = 0; i < n - 1; ++i) {
        term1 += std::abs((double)r[i + 1] - (double)r[i]);
    }
    
    double n_double = (double)n;
    double xi = 1.0 - (3.0 * term1) / (n_double * n_double - 1.0);
    return xi;
}

// ============================================================
// IAAFT Surrogate Generator
// ============================================================
// Generates Iterative Amplitude Adjusted Fourier Transform (IAAFT)
// surrogate datasets to test for non-linear dependence.
// [[Rcpp::export]]
arma::mat generate_iaaft_surrogates(arma::vec x, int n_surr, int max_iter = 100) {
    int n = x.n_elem;
    mat surrogates(n, n_surr);
    
    // Explicitly define std::complex for FFT operations
    cx_vec X_fft = fft(std::complex<double>(0,0) * zeros(n) + x);
    vec amplitudes = abs(X_fft); 
    
    vec x_sorted = sort(x);

    for(int s = 0; s < n_surr; s++) {
        vec surr = shuffle(x);
        
        for(int iter = 0; iter < max_iter; iter++) {
            cx_vec S_fft = fft(std::complex<double>(0,0) * zeros(n) + surr);
            
            // abs() returns a real vector (vec)
            vec S_abs = abs(S_fft); 
            
            // Element-wise division to extract phase
            cx_vec phases = S_fft / (S_abs + 1e-10); 
            
            // Reconstruct the spectrum with original amplitudes
            cx_vec S_new_fft = phases % amplitudes; 
            vec surr_spec = real(ifft(S_new_fft));
            
            // Rank ordering to adjust amplitudes to the original distribution
            uvec r = sort_index(surr_spec);
            uvec rank_indices = sort_index(r); 
            surr = x_sorted(rank_indices);
        }
        surrogates.col(s) = surr;
    }
    
    return surrogates;
}

// ============================================================
// Core Engine: Compute Xi Statistics for Lags
// ============================================================
//' Compute Xi-ACF for Multiple Lags (Core Engine)
//' 
//' Calculates Chatterjee's Xi coefficient for multiple lags and generates
//' IAAFT (Iterative Amplitude Adjusted Fourier Transform) surrogates to establish 
//' confidence intervals for non-linear dependence.
//' 
//' @param x A numeric vector (time series).
//' @param max_lag An integer specifying the maximum number of lags to compute.
//' @param n_surr An integer specifying the number of surrogate datasets to generate.
//' @return A list containing \code{xi_original} (the Xi coefficients for the original series) 
//' and \code{xi_surrogates} (a matrix of Xi coefficients for the surrogate datasets).
//' @export
// [[Rcpp::export]]
List compute_xi_acf_iaaft(NumericVector x, int max_lag, int n_surr) {
    vec y = as<vec>(x); // Convert R vector to Armadillo vector
    int n = y.n_elem;
    
    // Initialize output vector and fill with NA
    vec xi_original(max_lag);
    xi_original.fill(datum::nan); 
    
    // Initialize surrogate matrix and fill with NA
    mat xi_surrogates(max_lag, n_surr);
    xi_surrogates.fill(datum::nan);

    // Compute Xi for the original time series across lags
    for(int k = 1; k <= max_lag; k++) {
        if(n > k) {
            vec y_lagged = y.subvec(0, n - k - 1);
            vec y_target = y.subvec(k, n - 1);
            xi_original(k-1) = xi_coefficient(y_lagged, y_target);
        }
    }
    
    // Generate IAAFT surrogates
    mat surrogates = generate_iaaft_surrogates(y, n_surr);
    
    // Compute Xi for each surrogate dataset across lags
    for(int s = 0; s < n_surr; s++) {
        vec y_surr = surrogates.col(s);
        for(int k = 1; k <= max_lag; k++) {
            if(n > k) {
                vec y_s_lagged = y_surr.subvec(0, n - k - 1);
                vec y_s_target = y_surr.subvec(k, n - 1);
                xi_surrogates(k-1, s) = xi_coefficient(y_s_lagged, y_s_target);
            }
        }
    }
    
    return List::create(
        Named("xi_original") = xi_original,
        Named("xi_surrogates") = xi_surrogates
    );
}