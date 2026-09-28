#' SSDr-F
#'
#' Fit a spatial embedding using a Gaussian kernel and squared Frobenius loss.
#'
#' @param X Numeric matrix with spots in rows and features in columns.
#' @param coords Numeric coordinate matrix or data frame, in the same row order as `X`.
#' @param rank Positive integer embedding dimension.
#' @param bandwidth Positive kernel bandwidth after global min-max scaling of `coords`.
#' @param lambda Non-negative penalty weight; `NULL` uses the data-dependent default.
#' @param center If `TRUE`, double-center `X` and the kernel.
#' @param max_iter Positive integer iteration limit.
#' @param tol Non-negative relative convergence tolerance.
#'
#' @return An `ssdr_fit` object with the embedding in `U`.
#' @export
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(80), nrow = 20)
#' coords <- cbind(runif(20), runif(20))
#' fit <- ssdr_f(X, coords, rank = 2, bandwidth = 0.2, max_iter = 2)
#' dim(fit$U)
ssdr_f <- function(X,
                   coords,
                   rank,
                   bandwidth,
                   lambda = NULL,
                   center = FALSE,
                   max_iter = 100,
                   tol = 1e-2) {
  call <- match.call()

  X <- as.matrix(X)
  coords <- as.matrix(coords)
  if (!is.numeric(X)) {
    stop("`X` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(coords)) {
    stop("`coords` must be numeric.", call. = FALSE)
  }
  if (nrow(X) != nrow(coords)) {
    stop("`X` and `coords` must have the same number of rows.", call. = FALSE)
  }

  if (rank < 1 || rank > min(dim(X)) || rank != floor(rank)) {
    stop("`rank` must be an integer between 1 and min(dim(X)).", call. = FALSE)
  }
  rank <- as.integer(rank)

  if (bandwidth <= 0) {
    stop("`bandwidth` must be positive.", call. = FALSE)
  }
  if (!is.null(lambda) && lambda < 0) {
    stop("`lambda` must be non-negative.", call. = FALSE)
  }
  if (max_iter < 1 || max_iter != floor(max_iter)) {
    stop("`max_iter` must be a positive integer.", call. = FALSE)
  }
  max_iter <- as.integer(max_iter)
  if (tol < 0) {
    stop("`tol` must be non-negative.", call. = FALSE)
  }
  center <- isTRUE(center)

  coord_min <- min(coords)
  coord_max <- max(coords)
  coords <- (coords - coord_min) / (coord_max - coord_min)

  n <- nrow(X)
  p <- ncol(X)
  K <- gaussian_kernel(coords, bandwidth)

  if (center) {
    H <- diag(n) - matrix(1 / n, n, n)
    J <- diag(p) - matrix(1 / p, p, p)
    X_fit <- H %*% X %*% J
    K_fit <- H %*% K %*% H
  } else {
    X_fit <- X
    K_fit <- K
  }

  svd_result <- truncated_svd(X_fit, rank)
  svd_s <- svd_result$s
  svd_V <- svd_result$V

  if (is.null(lambda)) {
    ratio <- 0.1 / 0.9
    lambda <- ratio * (min(svd_s)^2) / (p * n) * sum(diag(K_fit))
  }

  result <- ssdr_f_optimize_cpp(
    X = X_fit,
    Sigma_svd = svd_s,
    V_svd = svd_V,
    K = K_fit,
    lambda = lambda,
    max_iter = max_iter,
    tol = tol
  )

  structure(
    list(
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
    ),
    class = "ssdr_fit"
  )
}
