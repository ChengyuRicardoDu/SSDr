#' SSDr-NN
#'
#' Fit a spatial embedding using a coordinate neural network and Poisson loss.
#'
#' @param X Count matrix with spots in rows and features in columns.
#' @param coords Spatial coordinates in the same row order as `X`.
#' @param rank Number of embedding components.
#' @param bandwidth Gaussian bandwidth.
#' @param lambda Weight of the spatial smoothness penalty.
#' @param hidden_dim Units per hidden layer.
#' @param hidden_layers Number of hidden layers; 0 gives an affine map.
#' @param learning_rate Adam learning rate.
#' @param max_iter Maximum number of training epochs.
#' @param patience Epochs without sufficient loss improvement before stopping.
#' @param tol Tolerance for improvement in training loss.
#' @param seed Random seed for torch.
#'
#' @return An `ssdr_fit` object with the embedding in `U`.
#' @export
#'
#' @examples
#' \dontrun{
#' X <- matrix(rpois(200, lambda = 3), nrow = 20)
#' coords <- cbind(runif(20), runif(20))
#' fit <- ssdr_nn(X, coords, rank = 2, max_iter = 5, patience = 2)
#' dim(fit$U)
#' }
ssdr_nn <- function(X,
                    coords,
                    rank = 5,
                    bandwidth = 0.2,
                    lambda = 1e-2,
                    hidden_dim = 64,
                    hidden_layers = 2,
                    learning_rate = 5e-3,
                    max_iter = 500,
                    patience = 50,
                    tol = 1e-5,
                    seed = 1) {
  call <- match.call()

  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("Package `torch` is required for ssdr_nn().", call. = FALSE)
  }

  X <- as.matrix(X)
  coords <- as.matrix(coords)
  if (!is.numeric(X) || any(X < 0)) {
    stop("`X` must be numeric and non-negative.", call. = FALSE)
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
  if (lambda < 0) {
    stop("`lambda` must be non-negative.", call. = FALSE)
  }
  if (hidden_dim < 1 || hidden_dim != floor(hidden_dim)) {
    stop("`hidden_dim` must be a positive integer.", call. = FALSE)
  }
  hidden_dim <- as.integer(hidden_dim)
  if (hidden_layers < 0 || hidden_layers != floor(hidden_layers)) {
    stop("`hidden_layers` must be a non-negative integer.", call. = FALSE)
  }
  hidden_layers <- as.integer(hidden_layers)
  if (learning_rate <= 0) {
    stop("`learning_rate` must be positive.", call. = FALSE)
  }
  if (max_iter < 1 || max_iter != floor(max_iter)) {
    stop("`max_iter` must be a positive integer.", call. = FALSE)
  }
  max_iter <- as.integer(max_iter)
  if (patience < 1 || patience != floor(patience)) {
    stop("`patience` must be a positive integer.", call. = FALSE)
  }
  patience <- as.integer(patience)
  if (tol < 0) {
    stop("`tol` must be non-negative.", call. = FALSE)
  }
  seed <- as.integer(seed)

  for (j in seq_len(ncol(coords))) {
    coord_min <- min(coords[, j])
    coord_max <- max(coords[, j])
    if (coord_max == coord_min) {
      coords[, j] <- 0
    } else {
      coords[, j] <- (coords[, j] - coord_min) / (coord_max - coord_min)
    }
  }

  n <- nrow(X)
  p <- ncol(X)
  d <- ncol(coords)

  # Spatial weights
  sqdist <- as.matrix(stats::dist(coords))^2
  W <- exp(-sqdist / (2 * bandwidth^2))
  diag(W) <- 0
  W <- W / rowSums(W)

  torch::torch_manual_seed(seed)

  counts <- torch::torch_tensor(X, dtype = torch::torch_float())
  positions <- torch::torch_tensor(coords, dtype = torch::torch_float())
  weights <- torch::torch_tensor(W, dtype = torch::torch_float())

  # Coordinate network
  layers <- list()
  input_dim <- d
  for (i in seq_len(hidden_layers)) {
    layers <- c(layers, list(
      torch::nn_linear(input_dim, hidden_dim),
      torch::nn_relu()
    ))
    input_dim <- hidden_dim
  }
  layers <- c(layers, list(torch::nn_linear(input_dim, rank)))
  U_net <- do.call(torch::nn_sequential, layers)

  # Model parameters
  V_raw <- torch::nn_parameter(torch::torch_randn(c(p, rank)) * 0.05)
  log_sigma <- torch::nn_parameter(torch::torch_zeros(c(rank)))
  b0 <- torch::nn_parameter(torch::torch_zeros(c(1)))
  alpha_raw <- torch::nn_parameter(torch::torch_zeros(c(n, 1)))
  beta_raw <- torch::nn_parameter(torch::torch_zeros(c(1, p)))

  # Poisson model and spatial penalty
  compute_loss <- function() {
    U_raw <- U_net(positions)
    U <- U_raw - U_raw$mean(dim = 1, keepdim = TRUE)
    V <- V_raw - V_raw$mean(dim = 1, keepdim = TRUE)
    sigma <- log_sigma$exp()
    alpha <- alpha_raw - alpha_raw$mean()
    beta <- beta_raw - beta_raw$mean()

    low_rank <- (U * sigma$unsqueeze(1))$matmul(V$t())
    eta <- b0 + alpha + beta + low_rank
    mu <- eta$exp()
    nll <- (mu - counts * eta)$mean()

    norm2 <- (U * U)$sum(dim = 2, keepdim = TRUE)
    distance2 <- norm2 + norm2$t() - 2 * U$matmul(U$t())
    smooth_pen <- (weights * distance2)$sum() / n
    loss <- nll + lambda * smooth_pen

    list(loss = loss, U = U, V = V, sigma = sigma)
  }

  # Adam optimization
  params <- c(
    U_net$parameters,
    list(V_raw, log_sigma, b0, alpha_raw, beta_raw)
  )

  optimizer <- torch::optim_adam(params = params, lr = learning_rate)

  best_loss <- Inf
  best_epoch <- 0
  no_improve <- 0
  best_fit <- NULL

  for (epoch in seq_len(max_iter)) {
    optimizer$zero_grad()

    out <- compute_loss()
    loss <- out$loss
    loss$backward()
    optimizer$step()

    current_loss <- loss$item()

    if (is.finite(current_loss) &&
        (!is.finite(best_loss) || current_loss < best_loss - tol * (1 + abs(best_loss)))) {
      best_loss <- current_loss
      best_epoch <- epoch
      no_improve <- 0
      best_fit <- lapply(out, function(x) x$detach())
    } else {
      no_improve <- no_improve + 1
    }

    if (no_improve >= patience) {
      break
    }
  }

  if (is.null(best_fit)) {
    stop("Training did not produce a finite objective.", call. = FALSE)
  }

  # Return the selected fit
  U_final <- as.matrix(best_fit$U)
  V_final <- as.matrix(best_fit$V)
  sigma_final <- as.numeric(best_fit$sigma)

  for (k in seq_len(ncol(U_final))) {
    idx_max <- which.max(abs(U_final[, k]))
    if (U_final[idx_max, k] < 0) {
      U_final[, k] <- -U_final[, k]
      V_final[, k] <- -V_final[, k]
    }
  }

  diagnostics <- list(best_epoch = best_epoch)

  structure(
    list(
      method = "ssdr_nn",
      U = U_final,
      V = V_final,
      sigma = sigma_final,
      objective = best_loss,
      iterations = epoch,
      parameters = list(
        rank = rank,
        bandwidth = bandwidth,
        lambda = lambda,
        hidden_dim = hidden_dim,
        hidden_layers = hidden_layers,
        learning_rate = learning_rate,
        max_iter = max_iter,
        patience = patience,
        tol = tol,
        seed = seed
      ),
      diagnostics = diagnostics,
      call = call
    ),
    class = "ssdr_fit"
  )
}
