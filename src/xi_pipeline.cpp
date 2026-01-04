#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// ============================================================
// Helper 1: Rの乱数シードと同期するシャッフル関数
// ============================================================
void shuffle_indices_with_r_seed(uvec& idx) {
    GetRNGstate();
    int n = idx.n_elem;
    for (int i = n - 1; i > 0; i--) {
        int j = floor(R::runif(0, 1) * (i + 1));
        uword tmp = idx[i];
        idx[i] = idx[j];
        idx[j] = tmp;
    }
    PutRNGstate();
}

// ============================================================
// Helper 2: ランダムタイブレーク付きランク計算
// ============================================================
uvec rank_random_ties_r_sync(vec x) {
    uword n = x.n_elem;
    uvec idx = sort_index(x);
    vec x_sorted = x(idx);
    
    for (uword i = 0; i < n; ) {
        uword j = i + 1;
        while (j < n && x_sorted[j] == x_sorted[i]) {
            j++;
        }
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
    
    double xi = 1.0 - (3.0 * term1) / (std::pow((double)n, 2.0) - 1.0);
    return xi;
}

// ============================================================
// IAAFT Surrogate Generator
// ============================================================
// [[Rcpp::export]]
arma::mat generate_iaaft_surrogates(arma::vec x, int n_surr, int max_iter = 100) {
    int n = x.n_elem;
    mat surrogates(n, n_surr);
    
    // std::complexを明示
    cx_vec X_fft = fft(std::complex<double>(0,0) * zeros(n) + x);
    vec amplitudes = abs(X_fft); 
    
    vec x_sorted = sort(x);
    GetRNGstate();

    for(int s = 0; s < n_surr; s++) {
        vec surr = shuffle(x);
        
        for(int iter = 0; iter < max_iter; iter++) {
            cx_vec S_fft = fft(std::complex<double>(0,0) * zeros(n) + surr);
            
            // 【修正箇所】ここです！ absの結果は実数ベクトル(vec)で受けます
            vec S_abs = abs(S_fft); 
            
            // vec + double は vec になり、cx_vec / vec は要素ごとの除算として機能します
            cx_vec phases = S_fft / (S_abs + 1e-10); 
            
            cx_vec S_new_fft = phases % amplitudes; 
            vec surr_spec = real(ifft(S_new_fft));
            
            uvec r = sort_index(surr_spec);
            uvec rank_indices = sort_index(r); 
            surr = x_sorted(rank_indices);
        }
        surrogates.col(s) = surr;
    }
    
    PutRNGstate();
    return surrogates;
}

// ============================================================
// Pipeline Wrapper
// ============================================================
//' @export
// [[Rcpp::export]]
List run_xi_test_cpp(NumericVector y_in, int max_lag, int n_surr) {
    vec y = as<vec>(y_in);
    int n = y.n_elem;
    
    vec xi_original(max_lag);
    mat xi_surrogates(max_lag, n_surr);
    
    for(int k = 1; k <= max_lag; k++) {
        if(n > k) {
            vec y_lagged = y.subvec(0, n - k - 1);
            vec y_target = y.subvec(k, n - 1);
            xi_original(k-1) = xi_coefficient(y_lagged, y_target);
        }
    }
    
    mat surrogates = generate_iaaft_surrogates(y, n_surr);
    
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