// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List ssdr_truncated_svd_cpp(const arma::mat& X, int rank) {
  if (X.n_rows == 0 || X.n_cols == 0) {
    stop("`X` must be a non-empty numeric matrix.");
  }
  if (!X.is_finite()) {
    stop("`X` contains non-finite values.");
  }
  if (rank < 1) {
    stop("`rank` must be a positive integer.");
  }

  arma::mat U;
  arma::mat V;
  arma::vec s;

  bool ok = arma::svd(U, s, V, X);
  if (!ok) {
    stop("SVD failed to converge.");
  }

  rank = std::min((int)s.n_elem, rank);

  return Rcpp::List::create(
    Rcpp::Named("U") = U.cols(0, rank - 1),
    Rcpp::Named("s") = s.subvec(0, rank - 1),
    Rcpp::Named("V") = V.cols(0, rank - 1)
  );
}
