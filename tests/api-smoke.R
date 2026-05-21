library(ssdr)

set.seed(1)

coords <- cbind(runif(12), runif(12))
X <- matrix(rnorm(60), nrow = 12)

fit_f <- ssdr_f(X, coords, rank = 2, bandwidth = 0.3, lambda = 0.1, max_iter = 2)
stopifnot(inherits(fit_f, "ssdr_fit"))
stopifnot(identical(dim(fit_f$U), c(12L, 2L)))

counts <- matrix(rpois(60, lambda = 3), nrow = 12)
fit_p <- ssdr_p(
  counts,
  coords,
  rank = 2,
  bandwidth = 0.3,
  lambda = 10,
  step_size_init = 1e-6,
  max_iter = 1,
  line_search_steps = 5
)
stopifnot(inherits(fit_p, "ssdr_fit"))
stopifnot(identical(dim(fit_p$U), c(12L, 2L)))
