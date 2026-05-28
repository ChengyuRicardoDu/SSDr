// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <functional>
#include <limits>

using namespace Rcpp;
using namespace arma;

static inline bool normalize_factor_with_sigma(arma::mat& factor,
                                               arma::mat& Sigma,
                                               double eps,
                                               bool allow_zero = false) {
  arma::vec sigma_diag = Sigma.diag();
  if (sigma_diag.n_elem != factor.n_cols) {
    return false;
  }

  for (uword j = 0; j < factor.n_cols; ++j) {
    double nrm = arma::norm(factor.col(j), 2);
    if (!std::isfinite(nrm) || nrm <= eps) {
      if (!allow_zero) {
        return false;
      }
      factor.col(j).zeros();
      sigma_diag(j) = 0.0;
      continue;
    }
    factor.col(j) /= nrm;
    sigma_diag(j) *= nrm;
  }

  Sigma = arma::diagmat(sigma_diag);
  return true;
}

static double line_search(const arma::mat& current,
                          const arma::mat& grad,
                          double step_size_init,
                          const std::function<double(const arma::mat&)>& obj,
                          int line_search_steps,
                          double tau = 0.6) {
  double f0 = obj(current);
  double best_alpha = 0.0;
  double best_f1 = f0;

  std::vector<double> alphas;
  alphas.reserve(line_search_steps);
  alphas.push_back(step_size_init);

  int n_large = line_search_steps / 2;
  int n_small = line_search_steps - 1 - n_large;

  double a_large = step_size_init;
  for (int i = 0; i < n_large; ++i) {
    a_large /= tau;
    if (std::isfinite(a_large)) {
      alphas.push_back(a_large);
    }
  }

  double a_small = step_size_init;
  for (int i = 0; i < n_small; ++i) {
    a_small *= tau;
    if (std::isfinite(a_small)) {
      alphas.push_back(a_small);
    }
  }

  std::sort(alphas.begin(), alphas.end(), std::greater<double>());

  for (double alpha : alphas) {
    arma::mat cand = current - alpha * grad;
    double f1 = obj(cand);
    if (std::isfinite(f1) && f1 < best_f1) {
      best_f1 = f1;
      best_alpha = alpha;
    }
  }

  return best_alpha;
}

static double poisson_nll(const arma::mat& low_rank,
                          const arma::mat& X,
                          const arma::mat& B) {
  if (low_rank.n_rows != X.n_rows || low_rank.n_cols != X.n_cols ||
      B.n_rows != X.n_rows || B.n_cols != X.n_cols) {
    stop("Dimension mismatch in `poisson_nll`.");
  }

  arma::mat E = B % arma::exp(low_rank);
  double term1 = arma::accu(E);
  double term2 = arma::accu(low_rank % X);
  double term3 = arma::accu(X % arma::log(B + 1e-16));

  return term1 - term2 - term3;
}

// [[Rcpp::export]]
Rcpp::List ssdr_p_optimize_cpp(const arma::mat& X,
                               arma::mat B,
                               const arma::vec& Sigma_svd,
                               const arma::mat& V_svd,
                               const arma::mat& K,
                               const arma::mat& A_init,
                               double lambda,
                               double step_size_init,
                               int max_iter,
                               double tol = 1e-3,
                               int line_search_steps = 50) {
  const uword n = X.n_rows;
  const uword p = X.n_cols;
  const uword rank = A_init.n_cols;
  const double np_inv = 1.0 / double(n * p);
  const double n_inv = 1.0 / double(n);
  const double norm_eps = std::sqrt(arma::datum::eps);

  if (n == 0 || p == 0) {
    stop("`X` must be a non-empty matrix.");
  }
  if (K.n_rows != n || K.n_cols != n) {
    stop("`K` must be an n x n matrix where n = nrow(X).");
  }
  if (B.n_rows != n || B.n_cols != p) {
    stop("`B` must have the same dimensions as `X`.");
  }
  if (V_svd.n_rows != p || V_svd.n_cols != rank || Sigma_svd.n_elem != rank) {
    stop("Initial SVD dimensions do not match `X` and `A_init`.");
  }
  if (!X.is_finite() || !B.is_finite() || !Sigma_svd.is_finite() ||
      !V_svd.is_finite() || !K.is_finite() || !A_init.is_finite()) {
    stop("Inputs to `ssdr_p_optimize_cpp` contain non-finite values.");
  }
  if (!std::isfinite(lambda) || !std::isfinite(step_size_init) || !std::isfinite(tol)) {
    stop("`lambda`, `step_size_init`, and `tol` must be finite.");
  }
  if (max_iter < 1 || line_search_steps < 1) {
    stop("`max_iter` and `line_search_steps` must be >= 1.");
  }
  if (lambda < 0.0 || step_size_init <= 0.0 || tol < 0.0) {
    stop("`lambda` and `tol` must be non-negative; `step_size_init` must be positive.");
  }

  arma::mat A = A_init;
  arma::mat V = V_svd;
  arma::mat Sigma = diagmat(Sigma_svd);

  bool ok = normalize_factor_with_sigma(A, Sigma, norm_eps, true);
  if (!ok) {
    stop("Dimension mismatch while normalizing initial factors.");
  }

  arma::mat U = K * A;

  auto objective = [&](const arma::mat& A_eval,
                       const arma::mat& V_eval,
                       const arma::mat& Sigma_eval) {
    arma::mat low_rank = (K * A_eval) * Sigma_eval * V_eval.t();
    double loss = poisson_nll(low_rank, X, B) * np_inv;
    double reg = (lambda * trace(A_eval.t() * (K * A_eval))) * n_inv;
    return loss + reg;
  };

  double f = objective(A, V, Sigma);
  int iterations = 0;

  for (int cycle = 0; cycle < max_iter; ++cycle) {
    bool improved = false;
    iterations = cycle + 1;

    for (int step = 0; step < 10; ++step) {
      arma::mat E = B % arma::exp(U * Sigma * V.t());
      arma::mat gradA =
        (-K * X * (V * Sigma.t()) + K * (E * (V * Sigma.t()))) * np_inv
        + 2.0 * (lambda * (K * A)) * n_inv;

      auto objA = [&](const arma::mat& Acan) {
        return objective(Acan, V, Sigma);
      };

      double alpha = line_search(A, gradA, step_size_init, objA, line_search_steps);
      if (alpha <= 0.0) {
        break;
      }

      arma::mat Atrial = A - alpha * gradA;
      double f2 = objA(Atrial);
      if (f2 < f - tol * std::abs(f)) {
        A = std::move(Atrial);
        U = K * A;
        f = objective(A, V, Sigma);
        improved = true;
      } else {
        break;
      }
    }

    for (int step = 0; step < 10; ++step) {
      arma::mat E = B % arma::exp(U * Sigma * V.t());
      arma::mat gradV = (-X.t() * U * Sigma + E.t() * U * Sigma) * np_inv;

      auto objV = [&](const arma::mat& Vcan) {
        return objective(A, Vcan, Sigma);
      };

      double alphaV = line_search(V, gradV, step_size_init, objV, line_search_steps);
      if (alphaV <= 0.0) {
        break;
      }

      arma::mat Vtrial = V - alphaV * gradV;
      double f2 = objV(Vtrial);
      if (f2 < f - tol * std::abs(f)) {
        V = std::move(Vtrial);
        f = objective(A, V, Sigma);
        improved = true;
      } else {
        break;
      }
    }

    {
      double f_before_offset = f;
      arma::mat low_rank = U * Sigma * V.t();
      arma::mat Mtmp = B % arma::exp(low_rank);

      for (int it = 0; it < 20; ++it) {
        arma::vec sumX_row = sum(X, 1);
        arma::vec sumM_row = arma::clamp(sum(Mtmp, 1), norm_eps, arma::datum::inf);
        arma::vec r_scale = sumX_row / sumM_row;
        Mtmp.each_col() %= r_scale;
        B.each_col() %= r_scale;

        arma::rowvec sumX_col = sum(X, 0);
        arma::rowvec sumM_col = arma::clamp(sum(Mtmp, 0), norm_eps, arma::datum::inf);
        arma::rowvec c_scale = sumX_col / sumM_col;
        Mtmp.each_row() %= c_scale;
        B.each_row() %= c_scale;
      }

      double f_after_offset = objective(A, V, Sigma);
      if (f_after_offset < f_before_offset - tol * std::abs(f_before_offset)) {
        improved = true;
        f = f_after_offset;
      }
    }

    {
      arma::mat A_old = A;
      arma::mat V_old = V;
      arma::mat Sigma_old = Sigma;
      double f_before_norm = f;

      bool okA = normalize_factor_with_sigma(A, Sigma, norm_eps, true);
      bool okV = normalize_factor_with_sigma(V, Sigma, norm_eps, true);

      if (!okA || !okV) {
        A = A_old;
        V = V_old;
        Sigma = Sigma_old;
      } else {
        U = K * A;
        double f_after_norm = objective(A, V, Sigma);
        if (std::isfinite(f_after_norm) &&
            f_after_norm <= f_before_norm + 1e-10 * (1.0 + std::abs(f_before_norm))) {
          f = f_after_norm;
        } else {
          A = A_old;
          V = V_old;
          Sigma = Sigma_old;
          U = K * A;
          f = objective(A, V, Sigma);
        }
      }
    }

    if (!improved) {
      break;
    }
  }

  U = K * A;
  return Rcpp::List::create(
    Rcpp::Named("U") = U,
    Rcpp::Named("V") = V,
    Rcpp::Named("Sigma") = Sigma,
    Rcpp::Named("objective") = f,
    Rcpp::Named("iterations") = iterations
  );
}
