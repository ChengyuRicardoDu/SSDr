// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat gaussian_kernel(const arma::mat& coordinates, double bandwidth,
                          Rcpp::Nullable<Rcpp::NumericMatrix> landmarks = R_NilValue) {
  arma::uword n = coordinates.n_rows;

  if (landmarks.isNotNull()) {
    arma::mat points = Rcpp::as<arma::mat>(landmarks.get());
    arma::mat D(n, points.n_rows);
    for (arma::uword i = 0; i < n; ++i) {
      for (arma::uword j = 0; j < points.n_rows; ++j) {
        D(i, j) = arma::accu(arma::square(coordinates.row(i) - points.row(j)));
      }
    }
    return arma::exp(-D / (2.0 * bandwidth * bandwidth));
  }

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
