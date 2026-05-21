#' Create an SSDr Fit Object
#'
#' Internal constructor used by the SSDr method wrappers.
#'
#' @param method Character scalar identifying the fitted method.
#' @param U Numeric matrix containing the low-dimensional embedding.
#' @param V Optional numeric matrix containing feature loadings.
#' @param sigma Optional numeric vector of component scales.
#' @param objective Optional final objective value.
#' @param iterations Optional number of optimization iterations.
#' @param parameters List of user-facing tuning parameters.
#' @param diagnostics List of additional diagnostics.
#' @param call Matched call.
#'
#' @return An object of class `ssdr_fit`.
#' @noRd
new_ssdr_fit <- function(method,
                         U,
                         V = NULL,
                         sigma = NULL,
                         objective = NULL,
                         iterations = NULL,
                         parameters = list(),
                         diagnostics = list(),
                         call = NULL) {
  structure(
    list(
      method = method,
      U = U,
      V = V,
      sigma = sigma,
      objective = objective,
      iterations = iterations,
      parameters = parameters,
      diagnostics = diagnostics,
      call = call
    ),
    class = c("ssdr_fit", "list")
  )
}

#' Print an SSDr Fit
#'
#' @param x An `ssdr_fit` object.
#' @param ... Additional arguments, currently unused.
#'
#' @return `x`, invisibly.
#' @export
print.ssdr_fit <- function(x, ...) {
  dims <- if (is.matrix(x$U)) paste0(nrow(x$U), " x ", ncol(x$U)) else "unknown"
  cat("<ssdr_fit>\n")
  cat("  method: ", x$method, "\n", sep = "")
  cat("  U:      ", dims, "\n", sep = "")
  if (!is.null(x$objective)) {
    cat("  objective: ", format(x$objective, digits = 6), "\n", sep = "")
  }
  if (!is.null(x$iterations)) {
    cat("  iterations: ", x$iterations, "\n", sep = "")
  }
  invisible(x)
}
