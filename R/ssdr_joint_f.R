kernel_exact <- function(coords, bandwidth, p_lambda) {
  n <- nrow(coords)
  K <- gaussian_kernel(coords, bandwidth)
  H <- diag(n) - matrix(1 / n, n, n)
  kernel <- H %*% K %*% H
  kernel <- (kernel + t(kernel)) / 2

  eig <- eigen(kernel, symmetric = TRUE)
  d <- pmax(eig$values, 0)
  list(
    vectors = eig$vectors,
    values = d,
    shrinkage = d / (d + p_lambda)
  )
}

kernel_nystrom <- function(coords, bandwidth, p_lambda, n_landmarks) {
  # Landmarks
  ids <- sample.int(nrow(coords), min(n_landmarks, nrow(coords)))
  landmarks <- coords[ids, , drop = FALSE]

  # Kernel matrices
  C <- gaussian_kernel(coords, bandwidth, landmarks)
  W <- C[ids, , drop = FALSE]

  # Nystrom basis
  eig <- eigen(W, symmetric = TRUE)
  keep <- eig$values > 1e-10 * eig$values[1L]
  Q <- eig$vectors[, keep, drop = FALSE]
  d <- eig$values[keep]
  inverse_root <- Q %*% diag(1 / sqrt(d), nrow = length(d))
  B <- C %*% inverse_root
  n <- nrow(B)
  H <- diag(n) - matrix(1 / n, n, n)
  B <- H %*% B

  sv <- truncated_svd(B, min(dim(B)))
  d <- sv$s^2
  list(
    vectors = sv$U,
    values = d,
    shrinkage = d / (d + p_lambda)
  )
}

#' Joint SSDr-F
#'
#' @param Y List of expression matrices.
#' @param coords List of spatial coordinate matrices.
#' @param rank Number of embedding components.
#' @param bandwidth Gaussian kernel bandwidth.
#' @param lambda Spatial penalty weight.
#' @param solver `"nystrom"` for approximation or `"exact"` for the full kernel.
#' @param n_landmarks Number of landmarks per sample for `"nystrom"`.
#' @param landmark_seed Seed for landmark selection.
#'
#' @return Sample embeddings in `U` and shared loadings in `V`.
#' @export
ssdr_joint_f <- function(Y,
                         coords,
                         rank,
                         bandwidth,
                         lambda,
                         solver = c("nystrom", "exact"),
                         n_landmarks = 200L,
                         landmark_seed = 1L) {
  call <- match.call()

  if (!is.list(Y) || !is.list(coords) ||
      length(Y) < 2 || length(coords) != length(Y)) {
    stop("`Y` and `coords` must be lists of the same length, with at least two samples.",
         call. = FALSE)
  }
  Y <- lapply(Y, as.matrix)
  coords <- lapply(coords, as.matrix)
  sample_names <- names(Y)
  feature_names <- colnames(Y[[1L]])
  n_features <- ncol(Y[[1L]])
  for (i in seq_along(Y)) {
    if (!is.numeric(Y[[i]]) || !is.numeric(coords[[i]]) ||
        ncol(Y[[i]]) != n_features || nrow(coords[[i]]) != nrow(Y[[i]])) {
      stop("Use numeric matrices with matching feature counts and coordinate rows.",
           call. = FALSE)
    }
  }

  solver <- match.arg(solver)
  if (rank < 1 || rank > n_features || rank != floor(rank) ||
      bandwidth <= 0 || lambda <= 0 ||
      (solver == "nystrom" && (n_landmarks < 1 || n_landmarks != floor(n_landmarks)))) {
    stop("Invalid rank, bandwidth, lambda or landmark count.", call. = FALSE)
  }

  n_samples <- length(Y)
  sample_means <- vector("list", n_samples)
  X <- Y
  scaled_coords <- coords
  for (i in seq_len(n_samples)) {
    n <- nrow(Y[[i]])
    H <- diag(n) - matrix(1 / n, n, n)
    sample_means[[i]] <- colMeans(Y[[i]])
    X[[i]] <- H %*% Y[[i]]
    dimnames(X[[i]]) <- dimnames(Y[[i]])

    coord_min <- apply(coords[[i]], 2, min)
    coord_max <- apply(coords[[i]], 2, max)
    coord_scale <- max(coord_max - coord_min)
    scaled_coords[[i]] <- H %*% coords[[i]] / coord_scale
    dimnames(scaled_coords[[i]]) <- dimnames(coords[[i]])
  }
  names(sample_means) <- sample_names

  p_lambda <- n_features * lambda
  M <- matrix(0, n_features, n_features)
  smoothed_X <- vector("list", n_samples)
  Z <- vector("list", n_samples)
  d <- vector("list", n_samples)

  if (solver == "nystrom") set.seed(landmark_seed)
  for (i in seq_len(n_samples)) {
    spectrum <- if (solver == "exact") {
      kernel_exact(scaled_coords[[i]], bandwidth, p_lambda)
    } else {
      kernel_nystrom(scaled_coords[[i]], bandwidth, p_lambda, n_landmarks)
    }

    Q <- spectrum$vectors
    d[[i]] <- spectrum$values
    D <- diag(spectrum$shrinkage, nrow = length(d[[i]]))
    Z[[i]] <- t(Q) %*% X[[i]]
    DZ <- D %*% Z[[i]]
    smoothed_X[[i]] <- Q %*% DZ
    M <- M + t(Z[[i]]) %*% DZ / (n_samples * nrow(X[[i]]))
  }
  M <- (M + t(M)) / 2

  eig <- eigen(M, symmetric = TRUE)
  V <- eig$vectors[, seq_len(rank), drop = FALSE]
  for (j in seq_len(ncol(V))) {
    anchor <- which.max(abs(V[, j]))
    if (V[anchor, j] < 0) V[, j] <- -V[, j]
  }
  rownames(V) <- feature_names
  colnames(V) <- paste0("SSDr", seq_len(rank))

  U <- vector("list", n_samples)
  objective <- 0
  for (i in seq_len(n_samples)) {
    U[[i]] <- smoothed_X[[i]] %*% V
    dimnames(U[[i]]) <- list(rownames(Y[[i]]), colnames(V))
    residual <- X[[i]] - U[[i]] %*% t(V)
    scores <- Z[[i]] %*% V
    A <- diag(d[[i]] / (d[[i]] + p_lambda)^2, nrow = length(d[[i]]))
    penalty <- sum(scores * (A %*% scores))
    n <- nrow(X[[i]])
    objective <- objective + sum(residual^2) / (n_samples * n * n_features) +
      lambda * penalty / (n_samples * n)
  }
  names(U) <- sample_names

  structure(
    list(
      method = "ssdr_joint_f",
      U = U,
      V = V,
      sample_means = sample_means,
      objective = objective,
      parameters = list(
        rank = rank,
        bandwidth = bandwidth,
        lambda = lambda,
        solver = solver,
        n_landmarks = if (solver == "nystrom") n_landmarks else NULL,
        landmark_seed = if (solver == "nystrom") landmark_seed else NULL
      ),
      call = call
    ),
    class = "ssdr_fit"
  )
}
