#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Forward declaration to use xi_coefficient from xi_core.cpp
double xi_coefficient(arma::vec x, arma::vec y);

// ============================================================
// Helper: Generate Shared Random Phases for Multivariate FT Surrogate
// ============================================================
// Generates a random phase vector that satisfies Hermitian symmetry
// so that the inverse FFT results in a real-valued time series.
vec generate_shared_random_phases(int n) {
    vec phases(n, fill::zeros);
    GetRNGstate(); // Pull RNG state from R
    
    // Generate random phases for positive frequencies
    for(int i = 1; i <= (n - 1) / 2; i++) {
        double p = R::runif(0, 2 * M_PI);
        phases(i) = p;
        phases(n - i) = -p; // Negative frequencies must be symmetric
    }
    
    // If n is even, the Nyquist frequency must be 0 (or pi) for real signals
    if(n % 2 == 0) {
        phases(n / 2) = 0;
    }
    
    PutRNGstate(); // Return RNG state to R
    return phases;
}

// ============================================================
// MIAAFT Core Algorithm (Phase Randomization Initialization + Iteration)
// ============================================================
//' Generate a Single MIAAFT Surrogate Matrix
//'
//' @param x A numeric matrix (rows = time, cols = variables).
//' @param max_iter An integer specifying the maximum number of iterations.
//' @return A numeric matrix representing the generated MIAAFT surrogate.
//' @export
// [[Rcpp::export]]
NumericMatrix generate_miaaft_surrogate_cpp(NumericMatrix x, int max_iter = 100) {
    mat Y = as<mat>(x);
    int n_rows = Y.n_rows;
    int n_cols = Y.n_cols;

    mat Y_sorted(n_rows, n_cols);
    mat Amplitudes(n_rows, n_cols);
    mat S(n_rows, n_cols);

    // --------------------------------------------------------
    // Step 1: Initialize via Shared Phase Randomization
    // --------------------------------------------------------
    vec shared_phases = generate_shared_random_phases(n_rows);

    for (int j = 0; j < n_cols; j++) {
        vec col_data = Y.col(j);
        
        // Target 1: Exact values (for Rank Replacement later)
        Y_sorted.col(j) = sort(col_data);
        
        // Transform original data to frequency domain
        cx_vec Y_freq = fft(col_data);
        
        // Target 2: Exact amplitudes
        vec amp = abs(Y_freq);
        Amplitudes.col(j) = amp;
        
        // Apply SHARED random phase to preserve cross-correlations!
        vec orig_phase = arg(Y_freq);
        vec new_phase = orig_phase + shared_phases;
        
        // Construct initial surrogate S^(0)
        cx_vec S_freq_init = cx_vec(amp % cos(new_phase), amp % sin(new_phase));
        S.col(j) = real(ifft(S_freq_init));
    }

    // --------------------------------------------------------
    // Step 2: Iteration Loop (Rank & Amplitude Replacement)
    // --------------------------------------------------------
    for (int iter = 0; iter < max_iter; iter++) {
        for (int j = 0; j < n_cols; j++) {
            // a) Transform surrogate to frequency domain
            cx_vec S_freq = fft(S.col(j));

            // b) Amplitude Replacement: Keep surrogate phases, use target amplitudes
            vec phases = arg(S_freq);
            vec target_amp = Amplitudes.col(j);
            cx_vec new_freq = cx_vec(target_amp % cos(phases), target_amp % sin(phases));

            // c) Inverse transform back to time domain
            vec S_time = real(ifft(new_freq));

            // d) Rank Replacement: Match ranks of S_time to target exact values
            uvec sort_idx = sort_index(S_time);
            vec updated_col(n_rows);
            updated_col(sort_idx) = Y_sorted.col(j);
            
            // Update the surrogate column
            S.col(j) = updated_col;
        }
    }

    return wrap(S);
}

// ============================================================
// Step 3: Multivariate Xi CCF (Cross-Correlation) Wrapper
// ============================================================
//' Compute Lagged Xi Cross-Correlation with MIAAFT Surrogates
//'
//' @param x A numeric vector (Variable 1)
//' @param y A numeric vector (Variable 2)
//' @param max_lag The maximum lag to compute (computes from -max_lag to +max_lag)
//' @param n_surr Number of MIAAFT surrogates to generate
//' @return A list containing lags, original xi values, and surrogate xi values.
//' @export
// [[Rcpp::export]]
List compute_xi_ccf_cpp(NumericVector x, NumericVector y, int max_lag, int n_surr) {
    vec xv = as<vec>(x);
    vec yv = as<vec>(y);
    int n = xv.n_elem;
    
    // Calculate number of lags (-max_lag to +max_lag)
    int n_lags = 2 * max_lag + 1;
    vec lags = regspace<vec>(-max_lag, max_lag);
    
    vec xi_original(n_lags, fill::zeros);
    mat xi_surrogates(n_lags, n_surr, fill::zeros);
    
    // --- 1. Compute Original CCF ---
    for(int i = 0; i < n_lags; i++) {
        int k = lags(i);
        if (k >= 0) {
            vec x_lagged = xv.subvec(0, n - k - 1);
            vec y_target = yv.subvec(k, n - 1);
            xi_original(i) = xi_coefficient(x_lagged, y_target);
        } else {
            int abs_k = -k; // Negative lag means Y leads X
            vec x_lagged = xv.subvec(abs_k, n - 1);
            vec y_target = yv.subvec(0, n - abs_k - 1);
            xi_original(i) = xi_coefficient(x_lagged, y_target);
        }
    }
    
    // --- 2. Generate Surrogates and Compute CCF ---
    // Combine x and y into a 2-column matrix
    mat Z(n, 2);
    Z.col(0) = xv;
    Z.col(1) = yv;
    NumericMatrix Z_rcpp = wrap(Z);
    
    for(int s = 0; s < n_surr; s++) {
        // Generate 1 MIAAFT surrogate preserving cross-correlation
        NumericMatrix Z_surr_rcpp = generate_miaaft_surrogate_cpp(Z_rcpp, 100);
        mat Z_surr = as<mat>(Z_surr_rcpp);
        
        vec x_surr = Z_surr.col(0);
        vec y_surr = Z_surr.col(1);
        
        for(int i = 0; i < n_lags; i++) {
            int k = lags(i);
            if (k >= 0) {
                vec x_lagged = x_surr.subvec(0, n - k - 1);
                vec y_target = y_surr.subvec(k, n - 1);
                xi_surrogates(i, s) = xi_coefficient(x_lagged, y_target);
            } else {
                int abs_k = -k;
                vec x_lagged = x_surr.subvec(abs_k, n - 1);
                vec y_target = y_surr.subvec(0, n - abs_k - 1);
                xi_surrogates(i, s) = xi_coefficient(x_lagged, y_target);
            }
        }
    }
    
    return List::create(
        Named("lags") = lags,
        Named("xi_original") = xi_original,
        Named("xi_surrogates") = xi_surrogates
    );
}