#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Shuffle indices using R's random number generator
//' @noRd
void shuffle_indices_with_r_seed(uvec& idx) {
    int n = idx.n_elem;
    for (int i = n - 1; i > 0; i--) {
        int j = floor(R::runif(0, 1) * (i + 1));
        uword tmp = idx[i];
        idx[i] = idx[j];
        idx[j] = tmp;
    }
}

//' Calculate ranks with random tie-breaking
//' @noRd
uvec rank_random_ties_r_sync(vec x) {
    uword n = x.n_elem;
    uvec idx = sort_index(x);
    vec x_sorted = x(idx);
    
    for (uword i = 0; i < n; ) {
        uword j = i + 1;
        while (j < n && x_sorted[j] == x_sorted[i]) {
            j++;
        }
        // Break ties randomly
        if (j > i + 1) {
            uvec tie_idx = idx.subvec(i, j - 1);
            shuffle_indices_with_r_seed(tie_idx);
            idx.subvec(i, j - 1) = tie_idx;
        }
        i = j;
    }
    
    uvec r(n);
    for (uword i = 0; i < n; i++) {
        r[idx[i]] = i + 1;
    }
    return r;
}

//' Compute Chatterjee's Xi Coefficient (Internal)
//' @noRd
double xi_coefficient(arma::vec x, arma::vec y) {
    int n = x.n_elem;
    if (n < 2) return NA_REAL;
    
    uvec r = rank_random_ties_r_sync(y);
    uvec ord = sort_index(x);
    uvec sorted_r = r(ord);
    
    double num = 0.0;
    for (int i = 0; i < n - 1; i++) {
        double diff = (double)sorted_r[i+1] - (double)sorted_r[i];
        num += std::abs(diff);
    }
    
    double den = (double)n * ((double)n - 1.0) / 2.0;
    return 1.0 - (num / den);
}

//' Compute Chatterjee's Xi coefficient (Exported to R)
//' @param x A numeric vector.
//' @param y A numeric vector.
//' @return The Xi coefficient.
//' @export
// [[Rcpp::export(name = "xi_coefficient")]]
double xi_coefficient_export(NumericVector x, NumericVector y) {
    return xi_coefficient(as<vec>(x), as<vec>(y));
}