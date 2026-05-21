#' Spatially Smoothed Dimension Reduction with Frobenius Loss (SSDr-F)
#'
#' Fits the kernel Frobenius formulation of SSDr. The embedding is written as
#' `U = K A`, where `K` is a Gaussian kernel over spatial coordinates.
#'
#' @param X Numeric matrix with spots or cells in rows and features in columns.
#'   For the manuscript SSDr-F analysis this is typically a PCA-reduced,
#'   log-normalized expression matrix.
#' @param coords Numeric matrix or data frame of spatial coordinates aligned
#'   with the rows of `X`.
#' @param rank Positive integer. Target embedding dimension.
#' @param bandwidth Positive numeric. Gaussian kernel bandwidth after
#'   global min-max normalization of `coords`.
#' @param lambda Non-negative numeric smoothness weight. If `NULL`, a
#'   data-dependent default is used.
#' @param center Logical. If `TRUE`, double-center both `X` and the kernel
#'   before fitting.
#' @param max_iter Positive integer. Maximum number of outer iterations.
#' @param tol Non-negative numeric. Relative convergence tolerance.
#' @param verbose Logical. If `TRUE`, print optimizer progress.
#' @param ... Deprecated argument aliases: `data_pixel`, `r`, `sigma`, and
#'   `double_center`.
#'
#' @return An `ssdr_fit` object. The embedding is stored in `fit$U`.
#' @export
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(80), nrow = 20)
#' coords <- cbind(runif(20), runif(20))
#' fit <- ssdr_f(X, coords, rank = 2, bandwidth = 0.2, max_iter = 2)
#' dim(fit$U)
ssdr_f <- function(X,
                   coords = NULL,
                   rank = NULL,
                   bandwidth = NULL,
                   lambda = NULL,
                   center = FALSE,
                   max_iter = 100,
                   tol = 1e-2,
                   verbose = FALSE,
                   ...) {
  call <- match.call()
  dots <- list(...)

  alias <- ssdr_arg_alias(dots, "coords", "data_pixel", coords, missing(coords))
  coords <- alias$value
  dots <- alias$dots
  alias <- ssdr_arg_alias(dots, "rank", "r", rank, missing(rank))
  rank <- alias$value
  dots <- alias$dots
  alias <- ssdr_arg_alias(dots, "bandwidth", "sigma", bandwidth, missing(bandwidth))
  bandwidth <- alias$value
  dots <- alias$dots
  alias <- ssdr_arg_alias(dots, "center", "double_center", center, missing(center))
  center <- alias$value
  dots <- alias$dots
  ssdr_check_unused_dots(dots)

  X <- ssdr_as_numeric_matrix(X, "X")
  if (is.null(coords)) stop("`coords` is required.", call. = FALSE)
  coords <- ssdr_normalize_coords_global(coords)
  if (nrow(X) != nrow(coords)) {
    stop("`nrow(X)` must equal `nrow(coords)`.", call. = FALSE)
  }
  if (is.null(rank)) stop("`rank` is required.", call. = FALSE)
  if (is.null(bandwidth)) stop("`bandwidth` is required.", call. = FALSE)
  rank <- ssdr_check_rank(rank, X)
  bandwidth <- ssdr_check_scalar(bandwidth, "bandwidth", lower = 0, strict = TRUE)
  if (!is.null(lambda)) {
    lambda <- ssdr_check_scalar(lambda, "lambda", lower = 0)
  }
  max_iter <- ssdr_check_integer_scalar(max_iter, "max_iter")
  tol <- ssdr_check_scalar(tol, "tol", lower = 0)
  center <- isTRUE(center)
  verbose <- isTRUE(verbose)

  n <- nrow(X)
  p <- ncol(X)
  K <- ssdr_gaussian_kernel_cpp(coords, bandwidth)

  if (center) {
    H <- diag(n) - matrix(1 / n, n, n)
    J <- diag(p) - matrix(1 / p, p, p)
    X_fit <- H %*% X %*% J
    K_fit <- H %*% K %*% H
  } else {
    X_fit <- X
    K_fit <- K
  }

  svd_result <- ssdr_truncated_svd_cpp(X_fit, rank)
  svd_s <- svd_result$s
  svd_V <- svd_result$V
  if (any(!is.finite(svd_s)) || any(svd_s <= sqrt(.Machine$double.eps))) {
    stop("`rank` is too large for the effective rank of `X`; reduce `rank`.", call. = FALSE)
  }

  if (is.null(lambda)) {
    ratio <- 0.1 / 0.9
    lambda <- ratio * (min(svd_s)^2) / (p * n) * sum(diag(K_fit))
    if (verbose) message(sprintf("Default lambda: %.5e", lambda))
  }

  result <- ssdr_f_optimize_cpp(
    X = X_fit,
    Sigma_svd = svd_s,
    V_svd = svd_V,
    K = K_fit,
    lambda = lambda,
    max_iter = max_iter,
    tol = tol,
    verbose = verbose
  )

  new_ssdr_fit(
    method = "ssdr_f",
    U = result$U,
    V = result$V,
    sigma = diag(result$Sigma),
    objective = result$objective,
    iterations = result$iterations,
    parameters = list(
      rank = rank,
      bandwidth = bandwidth,
      lambda = lambda,
      center = center,
      max_iter = max_iter,
      tol = tol
    ),
    diagnostics = list(reconstruction_loss = result$reconstruction_loss),
    call = call
  )
}
