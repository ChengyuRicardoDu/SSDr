ssdr_arg_alias <- function(dots, target, alias, current, target_missing = is.null(current)) {
  if (!alias %in% names(dots)) {
    return(list(value = current, dots = dots))
  }
  if (!target_missing) {
    stop("Use only one of `", target, "` and deprecated `", alias, "`.", call. = FALSE)
  }
  warning("`", alias, "` is deprecated; use `", target, "`.", call. = FALSE)
  current <- dots[[alias]]
  dots[[alias]] <- NULL
  list(value = current, dots = dots)
}

ssdr_check_unused_dots <- function(dots) {
  if (length(dots) > 0) {
    stop("Unused argument(s): ", paste(names(dots), collapse = ", "), call. = FALSE)
  }
}

ssdr_as_numeric_matrix <- function(x, name, nonnegative = FALSE) {
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("`", name, "` must be numeric.", call. = FALSE)
  }
  if (nrow(x) == 0 || ncol(x) == 0) {
    stop("`", name, "` must be non-empty.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("`", name, "` contains non-finite values.", call. = FALSE)
  }
  if (nonnegative && any(x < 0)) {
    stop("`", name, "` must be non-negative.", call. = FALSE)
  }
  x
}

ssdr_normalize_coords <- function(coords) {
  coords <- ssdr_as_numeric_matrix(coords, "coords")
  out <- coords
  for (j in seq_len(ncol(coords))) {
    mn <- min(coords[, j])
    mx <- max(coords[, j])
    if (abs(mx - mn) < .Machine$double.eps) {
      out[, j] <- 0
    } else {
      out[, j] <- (coords[, j] - mn) / (mx - mn)
    }
  }
  out
}

ssdr_check_rank <- function(rank, X) {
  if (!is.numeric(rank) || length(rank) != 1 || !is.finite(rank)) {
    stop("`rank` must be a positive integer.", call. = FALSE)
  }
  rank <- as.integer(rank)
  if (rank < 1) {
    stop("`rank` must be a positive integer.", call. = FALSE)
  }
  if (rank > min(dim(X))) {
    stop("`rank` cannot exceed min(dim(X)).", call. = FALSE)
  }
  rank
}

ssdr_check_scalar <- function(x, name, lower = -Inf, strict = FALSE) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x)) {
    stop("`", name, "` must be a finite numeric scalar.", call. = FALSE)
  }
  ok <- if (strict) x > lower else x >= lower
  if (!ok) {
    op <- if (strict) ">" else ">="
    stop("`", name, "` must be ", op, " ", lower, ".", call. = FALSE)
  }
  x
}

ssdr_check_integer_scalar <- function(x, name, lower = 1) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x)) {
    stop("`", name, "` must be an integer scalar.", call. = FALSE)
  }
  x <- as.integer(x)
  if (x < lower) {
    stop("`", name, "` must be >= ", lower, ".", call. = FALSE)
  }
  x
}
