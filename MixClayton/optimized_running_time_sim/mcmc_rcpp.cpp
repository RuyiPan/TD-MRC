#include <Rcpp.h>
using namespace Rcpp;

inline double clayton_log_density(double u1, double u2, double theta) {
  if (theta < 1e-10) return 0.0;
  if (u1 <= 0.0 || u2 <= 0.0) return R_NegInf;
  const double log_u1 = std::log(u1);
  const double log_u2 = std::log(u2);
  const double term = std::exp(-theta * log_u1) + std::exp(-theta * log_u2) - 1.0;
  if (term <= 0.0 || !R_finite(term)) return R_NegInf;
  return std::log(theta + 1.0) -
    (theta + 1.0) * (log_u1 + log_u2) -
    (1.0 / theta + 2.0) * std::log(term);
}

// [[Rcpp::export]]
NumericMatrix clayton_rot_logdens_cpp(const NumericMatrix& U, const NumericVector& theta) {
  const int n = U.nrow();
  NumericMatrix out(n, 4);

  for (int i = 0; i < n; ++i) {
    const double u = U(i, 0);
    const double v = U(i, 1);

    out(i, 0) = clayton_log_density(u,       v,       theta[0]);
    out(i, 1) = clayton_log_density(1.0 - u, v,       theta[1]);
    out(i, 2) = clayton_log_density(1.0 - u, 1.0 - v, theta[2]);
    out(i, 3) = clayton_log_density(u,       1.0 - v, theta[3]);
  }

  return out;
}

// [[Rcpp::export]]
IntegerMatrix sample_z_cpp(const NumericMatrix& U,
                           const NumericVector& theta,
                           const NumericVector& pi) {
  const int n = U.nrow();
  IntegerMatrix Z(n, 4);
  NumericMatrix ll = clayton_rot_logdens_cpp(U, theta);

  for (int i = 0; i < n; ++i) {
    double weights[4];
    double total = 0.0;

    for (int k = 0; k < 4; ++k) {
      weights[k] = pi[k] * std::exp(ll(i, k));
      total += weights[k];
    }

    if (!(total > 0.0) || !R_finite(total)) {
      for (int k = 0; k < 4; ++k) weights[k] = pi[k];
      total = pi[0] + pi[1] + pi[2] + pi[3];
    }
    for (int k = 0; k < 4; ++k) weights[k] /= total;

    int draw[4] = {0, 0, 0, 0};
    R::rmultinom(1, weights, 4, draw);
    for (int k = 0; k < 4; ++k) {
      Z(i, k) = draw[k];
    }
  }

  return Z;
}

// [[Rcpp::export]]
NumericVector theta_component_loglik_cpp(const NumericMatrix& U,
                                         const NumericVector& theta,
                                         const IntegerMatrix& Z) {
  const int n = U.nrow();
  NumericVector sums(4);

  for (int i = 0; i < n; ++i) {
    const double u = U(i, 0);
    const double v = U(i, 1);

    if (Z(i, 0) == 1) {
      sums[0] += clayton_log_density(u, v, theta[0]);
    } else if (Z(i, 1) == 1) {
      sums[1] += clayton_log_density(1.0 - u, v, theta[1]);
    } else if (Z(i, 2) == 1) {
      sums[2] += clayton_log_density(1.0 - u, 1.0 - v, theta[2]);
    } else if (Z(i, 3) == 1) {
      sums[3] += clayton_log_density(u, 1.0 - v, theta[3]);
    }
  }

  return sums;
}

// [[Rcpp::export]]
NumericVector z_col_sums_cpp(const IntegerMatrix& Z) {
  NumericVector sums(4);
  const int n = Z.nrow();

  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < 4; ++k) {
      sums[k] += Z(i, k);
    }
  }

  return sums;
}
