// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List truncated_svd(const arma::mat& X, int rank) {
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
