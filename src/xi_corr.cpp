// xi_corr.cpp
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// ヘルパー：ランク計算 (平均ランクなどの複雑な処理は割愛し、高速化のため単純ランクを使用)
std::vector<double> get_ranks(const std::vector<double>& v) {
    int n = v.size();
    std::vector<std::pair<double, int>> v_idx(n);
    for(int i = 0; i < n; ++i) {
        v_idx[i] = {v[i], i};
    }
    
    // 値でソート
    std::sort(v_idx.begin(), v_idx.end());
    
    // ランクを元の位置に戻す
    std::vector<double> ranks(n);
    for(int i = 0; i < n; ++i) {
        ranks[v_idx[i].second] = i + 1; // 1-based rank
    }
    return ranks;
}

// Chatterjee係数単体の計算
// [[Rcpp::export]]
double calc_xi_cpp(NumericVector x, NumericVector y) {
    int n = x.size();
    if (n < 2) return 0.0;

    // 1. xに基づいてyを並べ替える
    std::vector<std::pair<double, double>> xy(n);
    for(int i = 0; i < n; ++i) {
        xy[i] = {x[i], y[i]};
    }
    std::sort(xy.begin(), xy.end()); // x (first) でソートされる

    // 2. ソートされたyを取り出す
    std::vector<double> y_sorted(n);
    for(int i = 0; i < n; ++i) {
        y_sorted[i] = xy[i].second;
    }

    // 3. ランク計算
    std::vector<double> r = get_ranks(y_sorted);

    // 4. 隣接ランクの差の絶対値の和
    double sum_diff = 0.0;
    for(int i = 0; i < n - 1; ++i) {
        sum_diff += std::abs(r[i+1] - r[i]);
    }

    // 5. 公式適用
    double xi = 1.0 - (3.0 * sum_diff) / (std::pow(n, 2) - 1.0);
    return xi;
}

// ラグごとのXiを計算してベクトルで返す関数
// [[Rcpp::export]]
NumericVector xi_correlogram_cpp(NumericVector ts, int max_lag) {
    NumericVector results(max_lag);
    int n = ts.size();
    
    for(int k = 1; k <= max_lag; ++k) {
        // ラグデータの作成 (Rの x[1:(n-k)] と x[(1+k):n])
        // C++なのでインデックスは 0 から n-1
        int len = n - k;
        NumericVector x_lag(len);
        NumericVector y_lag(len);
        
        for(int i = 0; i < len; ++i) {
            x_lag[i] = ts[i];      // x_t
            y_lag[i] = ts[i + k];  // x_{t+k}
        }
        
        results[k - 1] = calc_xi_cpp(x_lag, y_lag);
    }
    
    return results;
}