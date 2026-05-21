#' Spatially Smoothed Dimension Reduction with Kernel Poisson Loss (SSDr-P)
#'
#' Fits the supplementary kernel-Poisson formulation of SSDr. This method is
#' retained as an exported experimental interface because it is described in
#' the supplementary algorithm, but the main manuscript workflow emphasizes
#' `ssdr_f()` and `ssdr_nn()`.
#'
#' @param X Non-negative numeric count matrix with spots or cells in rows and
#'   genes or features in columns.
#' @param coords Numeric matrix or data frame of spatial coordinates aligned
#'   with the rows of `X`.
#' @param rank Positive integer. Target embedding dimension.
#' @param bandwidth Positive numeric. Gaussian kernel bandwidth after
#'   column-wise min-max normalization of `coords`.
#' @param lambda Non-negative numeric smoothness weight.
#' @param step_size_init Positive numeric initial line-search step size.
#' @param max_iter Positive integer. Maximum number of optimization cycles.
#' @param line_search_steps Positive integer. Number of line-search candidates.
#' @param tol Non-negative numeric convergence tolerance.
#' @param verbose Logical. If `TRUE`, print optimizer progress.
#' @param ... Deprecated argument aliases: `data_pixel`, `r`, `sigma`,
#'   `initial_alpha`, and `ls_candidates`.
#'
#' @return An `ssdr_fit` object. The embedding is stored in `fit$U`.
#' @export
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rpois(120, lambda = 3), nrow = 20)
#' coords <- cbind(runif(20), runif(20))
#' fit <- ssdr_p(X, coords, rank = 2, bandwidth = 0.2, max_iter = 1)
#' dim(fit$U)
ssdr_p <- function(X,
                   coords = NULL,
                   rank = 5,
                   bandwidth = 0.2,
                   lambda = 1e5,
                   step_size_init = 1e-5,
                   max_iter = 20,
                   line_search_steps = 50,
                   tol = 1e-3,
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
  alias <- ssdr_arg_alias(
    dots,
    "step_size_init",
    "initial_alpha",
    step_size_init,
    missing(step_size_init)
  )
  step_size_init <- alias$value
  dots <- alias$dots
  alias <- ssdr_arg_alias(
    dots,
    "line_search_steps",
    "ls_candidates",
    line_search_steps,
    missing(line_search_steps)
  )
  line_search_steps <- alias$value
  dots <- alias$dots
  ssdr_check_unused_dots(dots)

  X <- ssdr_as_numeric_matrix(X, "X", nonnegative = TRUE)
  if (is.null(coords)) stop("`coords` is required.", call. = FALSE)
  coords <- ssdr_normalize_coords(coords)
  if (nrow(X) != nrow(coords)) {
    stop("`nrow(X)` must equal `nrow(coords)`.", call. = FALSE)
  }
  rank <- ssdr_check_rank(rank, X)
  bandwidth <- ssdr_check_scalar(bandwidth, "bandwidth", lower = 0, strict = TRUE)
  lambda <- ssdr_check_scalar(lambda, "lambda", lower = 0)
  step_size_init <- ssdr_check_scalar(step_size_init, "step_size_init", lower = 0, strict = TRUE)
  max_iter <- ssdr_check_integer_scalar(max_iter, "max_iter")
  line_search_steps <- ssdr_check_integer_scalar(line_search_steps, "line_search_steps")
  tol <- ssdr_check_scalar(tol, "tol", lower = 0)
  verbose <- isTRUE(verbose)

  n <- nrow(X)
  p <- ncol(X)
  K_raw <- ssdr_gaussian_kernel_cpp(coords, bandwidth)
  logX <- log(X + 1)

  r_mean <- rowMeans(logX)
  c_mean <- colMeans(logX)
  g_mean <- mean(r_mean)
  alpha <- r_mean - g_mean
  beta <- c_mean - g_mean
  B <- exp(g_mean) * outer(exp(alpha), exp(beta), "*")

  H <- diag(n) - matrix(1 / n, n, n)
  J <- diag(p) - matrix(1 / p, p, p)
  K <- H %*% K_raw %*% H

  svd_result <- ssdr_truncated_svd_cpp(H %*% logX %*% J, rank)
  svd_U <- svd_result$U
  svd_s <- svd_result$s
  svd_V <- svd_result$V

  gamma <- 1e-10
  K_reg <- K + diag(gamma, nrow(K))
  A_init <- solve(K_reg, svd_U)

  projection_error <- norm(K %*% A_init - svd_U, "F") / norm(svd_U, "F")
  if (verbose) {
    message(sprintf("relative projection error = %.5e", projection_error))
  }

  col_norms <- sqrt(colSums(A_init^2)) / n
  col_norms[col_norms <= sqrt(.Machine$double.eps)] <- 1
  A_init <- sweep(A_init, 2, col_norms, "/")
  svd_s <- svd_s * col_norms

  result <- ssdr_p_optimize_cpp(
    X = X,
    B = B,
    Sigma_svd = svd_s,
    V_svd = svd_V,
    K = K,
    A_init = A_init,
    lambda = lambda,
    step_size_init = step_size_init,
    max_iter = max_iter,
    tol = tol,
    line_search_steps = line_search_steps,
    verbose = verbose
  )

  new_ssdr_fit(
    method = "ssdr_p",
    U = result$U,
    V = result$V,
    sigma = diag(result$Sigma),
    objective = result$objective,
    iterations = result$iterations,
    parameters = list(
      rank = rank,
      bandwidth = bandwidth,
      lambda = lambda,
      step_size_init = step_size_init,
      max_iter = max_iter,
      line_search_steps = line_search_steps,
      tol = tol
    ),
    diagnostics = list(projection_error = projection_error),
    call = call
  )
}
