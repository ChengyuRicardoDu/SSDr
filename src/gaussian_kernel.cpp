// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::export]]
arma::mat ssdr_gaussian_kernel_cpp(const arma::mat& coordinates, double bandwidth) {
  if (coordinates.n_rows == 0 || coordinates.n_cols == 0) {
    Rcpp::stop("`coordinates` must be a non-empty numeric matrix.");
  }
  if (!coordinates.is_finite()) {
    Rcpp::stop("`coordinates` contains non-finite values.");
  }
  if (!std::isfinite(bandwidth) || bandwidth <= 0.0) {
    Rcpp::stop("`bandwidth` must be a positive finite number.");
  }

  arma::uword n = coordinates.n_rows;
  arma::mat D(n, n, arma::fill::zeros);

  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = i + 1; j < n; ++j) {
      double d2 = arma::accu(arma::square(coordinates.row(i) - coordinates.row(j)));
      D(i, j) = d2;
      D(j, i) = d2;
    }
  }

  return arma::exp(-D / (2.0 * bandwidth * bandwidth));
}
