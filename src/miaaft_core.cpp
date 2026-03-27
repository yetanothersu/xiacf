#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// ============================================================
// Helper 1: FFT / IFFT and Complex Number Core Test
// ============================================================
// A simple function to test if Armadillo's FFT and complex number handling
// work correctly. It takes a vector, performs FFT, extracts amplitudes,
// and performs IFFT to reconstruct the signal.

//' Test FFT and IFFT in C++
//'
//' @param x A numeric vector.
//' @return A list containing the original, amplitudes, phases, and reconstructed vector.
//' @export
// [[Rcpp::export]]
List test_fft_cpp(NumericVector x) {
    // Convert R vector to Armadillo vector (Matching xi_core.cpp style)
    vec y = as<vec>(x); 
    
    // 1. Perform Fast Fourier Transform (Output is Complex Vector)
    cx_vec Y_freq = fft(y);
    
    // 2. Extract Amplitudes (Magnitudes)
    vec amplitudes = abs(Y_freq);
    
    // 3. Extract Phases
    vec phases = arg(Y_freq);
    
    // 4. Reconstruct the signal using Inverse FFT (IFFT)
    vec y_reconstructed = real(ifft(Y_freq));
    
    // Return the results to R for inspection
    return List::create(
        Named("original") = y,
        Named("amplitudes") = amplitudes,
        Named("phases") = phases,
        Named("reconstructed") = y_reconstructed
    );
}