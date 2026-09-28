// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

using namespace Rcpp;

static inline double frobenius_loss(const arma::mat& X,
                                    const arma::mat& U,
                                    const arma::vec& sigma,
                                    const arma::mat& V) {
  arma::mat R = X - (U * arma::diagmat(sigma)) * V.t();
  return arma::accu(arma::square(R));
}

static arma::mat update_a_block(const arma::vec& sigma,
                                const arma::mat& X,
                                const arma::mat& V,
                                const arma::mat& K,
                                double lambda) {
  arma::mat residual = X;
  const double n = residual.n_rows;
  const double p = residual.n_cols;
  const int rank = sigma.n_elem;
  arma::mat A(K.n_rows, rank, arma::fill::zeros);

  for (int k = 0; k < rank; ++k) {
    double sigma_k = sigma[k];
    if (sigma_k == 0.0) {
      stop("`sigma` must be nonzero.");
    }

    arma::vec v = V.col(k);
    arma::mat M = sigma_k * K + lambda * p * arma::eye<arma::mat>(n, n) / sigma_k;
    arma::colvec rhs = residual * v;

    arma::colvec a_k = arma::solve(M, rhs, arma::solve_opts::fast);

    A.col(k) = a_k;
    residual -= (K * a_k) * (sigma_k * v.t());
  }

  return A;
}

// [[Rcpp::export]]
Rcpp::List ssdr_f_optimize_cpp(const arma::mat& X,
                               const arma::vec& Sigma_svd,
                               const arma::mat& V_svd,
                               const arma::mat& K,
                               double lambda,
                               int max_iter = 100,
                               double tol = 1e-2) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  const int rank = Sigma_svd.n_elem;
  const double np = (double)n * p;
  const double n_inv = 1.0 / (double)n;

  arma::vec sigma = Sigma_svd;
  arma::mat Sigma = arma::diagmat(sigma);
  arma::mat V = V_svd;
  arma::mat A(n, rank, arma::fill::zeros);
  arma::mat U(n, rank, arma::fill::zeros);

  double best_objective = std::numeric_limits<double>::infinity();
  arma::mat best_U = U;
  arma::mat best_V = V;
  arma::vec best_sigma = sigma;
  int iterations = 0;

  for (int iter = 0; iter < max_iter; ++iter) {
    // (a) A-block: closed-form ridge solve for A, then U = K A
    A = update_a_block(sigma, X, V, K, lambda);
    U = K * A;

    // (b) V-block: least-squares update of V given A and Sigma
    arma::mat M = Sigma.t() * A.t() * K * K * A * Sigma;
    arma::mat RHS = X.t() * K * A * Sigma;

    arma::mat Vt = arma::solve(M.t(), RHS.t(), arma::solve_opts::fast);
    V = Vt.t();

    // (c) Normalize V columns, folding the norms into Sigma
    for (int k = 0; k < rank; ++k) {
      double norm_val = arma::norm(V.col(k), 2);
      if (norm_val > arma::datum::eps) {
        V.col(k) /= norm_val;
        Sigma.col(k) *= norm_val;
      }
    }

    sigma = arma::diagvec(Sigma);

    // (d) Penalized Frobenius objective and convergence check
    double objective = frobenius_loss(X, U, sigma, V) / np
      + lambda * n_inv * arma::trace(A.t() * K * A);

    iterations = iter + 1;
    if (objective < best_objective) {
      if (std::isfinite(best_objective)) {
        double rel_change = std::abs(best_objective - objective) / std::abs(best_objective);
        if (rel_change < tol) {
          break;
        }
      }

      best_objective = objective;
      best_U = U;
      best_V = V;
      best_sigma = sigma;
    } else {
      break;
    }
  }

  double reconstruction_loss = frobenius_loss(X, best_U, best_sigma, best_V) / np;

  return Rcpp::List::create(
    Rcpp::Named("U") = best_U,
    Rcpp::Named("V") = best_V,
    Rcpp::Named("Sigma") = arma::diagmat(best_sigma),
    Rcpp::Named("objective") = best_objective,
    Rcpp::Named("reconstruction_loss") = reconstruction_loss,
    Rcpp::Named("iterations") = iterations
  );
}
