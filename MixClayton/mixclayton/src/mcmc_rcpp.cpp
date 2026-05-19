#include <Rcpp.h>
using namespace Rcpp;

inline double clayton_log_density(const std::vector<double>& u, double theta) {
  if (theta < 1e-10) return 0.0;

  const int m = u.size();
  double log_prod = 0.0;
  double sum_power = 0.0;
  double log_const = 0.0;

  for (int j = 0; j < m; ++j) {
    if (u[j] <= 0.0 || u[j] >= 1.0) return R_NegInf;
    const double log_uj = std::log(u[j]);
    log_prod += log_uj;
    sum_power += std::exp(-theta * log_uj);
  }

  for (int j = 1; j < m; ++j) {
    log_const += std::log(1.0 + j * theta);
  }

  const double term = sum_power - m + 1.0;

  if (term <= 0.0 || !R_finite(term)) return R_NegInf;

  return log_const -
    (theta + 1.0) * log_prod -
    (1.0 / theta + m) * std::log(term);
}

// [[Rcpp::export]]
NumericMatrix clayton_rot_logdens_cpp(const NumericMatrix& U, const NumericVector& theta) {
  const int n = U.nrow();
  const int m = U.ncol();
  const int K = theta.size();
  const int expected_K = 1 << m;

  if (K != expected_K) {
    stop("Length of `theta` must equal 2^m, where m = ncol(U).");
  }

  NumericMatrix out(n, K);
  std::vector<double> u_rot(m);

  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < K; ++k) {

      for (int j = 0; j < m; ++j) {
        const double uj = U(i, j);
        const bool reflected = k & (1 << (m - 1 - j));
        u_rot[j] = reflected ? 1.0 - uj : uj;
      }

      out(i, k) = clayton_log_density(u_rot, theta[k]);
    }
  }

  return out;
}

// [[Rcpp::export]]
IntegerMatrix sample_z_cpp(const NumericMatrix& U,
                           const NumericVector& theta,
                           const NumericVector& pi) {
  const int n = U.nrow();
  const int K = theta.size();
  IntegerMatrix Z(n, K);
  NumericMatrix ll = clayton_rot_logdens_cpp(U, theta);

  for (int i = 0; i < n; ++i) {
    std::vector<double> weights(K);
    double total = 0.0;

    for (int k = 0; k < K; ++k) {
      weights[k] = pi[k] * std::exp(ll(i, k));
      total += weights[k];
    }

    if (!(total > 0.0) || !R_finite(total)) {
      total = 0.0;
      for (int k = 0; k < K; ++k) {
        weights[k] = pi[k];
        total += pi[k];
      }
    }

    for (int k = 0; k < K; ++k) weights[k] /= total;

    std::vector<int> draw(K);
    R::rmultinom(1, weights.data(), K, draw.data());
    for (int k = 0; k < K; ++k) Z(i, k) = draw[k];
  }

  return Z;
}

// [[Rcpp::export]]
NumericVector theta_component_loglik_cpp(const NumericMatrix& U,
                                         const NumericVector& theta,
                                         const IntegerMatrix& Z) {
  const int n = U.nrow();
  const int m = U.ncol();
  const int K = theta.size();
  NumericVector sums(K);
  std::vector<double> u_rot(m);

  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < K; ++k) {
      if (Z(i, k) == 1) {
        for (int j = 0; j < m; ++j) {
          const double uj = U(i, j);
          const bool reflected = k & (1 << (m - 1 - j));
          u_rot[j] = reflected ? 1.0 - uj : uj;
        }

        sums[k] += clayton_log_density(u_rot, theta[k]);
        break;
      }
    }
  }

  return sums;
}

// [[Rcpp::export]]
NumericVector z_col_sums_cpp(const IntegerMatrix& Z) {
  const int K = Z.ncol();
  NumericVector sums(K);
  const int n = Z.nrow();

  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < K; ++k) {
      sums[k] += Z(i, k);
    }
  }

  return sums;
}
