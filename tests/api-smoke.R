library(ssdr)

set.seed(1)

coords <- cbind(runif(12), runif(12))
X <- matrix(rnorm(60), nrow = 12)

fit_f <- ssdr_f(X, coords, rank = 2, bandwidth = 0.3, lambda = 0.1, max_iter = 2)
stopifnot(inherits(fit_f, "ssdr_fit"))
stopifnot(identical(dim(fit_f$U), c(12L, 2L)))

coord_df <- data.frame(x = c(10, 20), y = c(100, 200))
global_coords <- ssdr:::ssdr_normalize_coords_global(coord_df)
column_coords <- ssdr:::ssdr_normalize_coords(coord_df)
expected_global <- matrix(
  c(0, 0.0526315789473684, 0.473684210526316, 1),
  nrow = 2
)
stopifnot(isTRUE(all.equal(unname(global_coords), expected_global)))
stopifnot(!isTRUE(all.equal(global_coords, column_coords)))

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
