# SSDr

Spatially Smoothed Dimension Reduction (SSDr) provides R implementations of
spatially regularized dimension-reduction methods for spatial transcriptomics.

## Installation

```r
remotes::install_github("ChengyuRicardoDu/SSDr")
```

`ssdr_nn()` requires the optional R `torch` package.

## Methods

| Function | Manuscript name | Input | Objective |
|---|---|---|---|
| `ssdr_f()` | SSDr-F | normalized expression or PCA scores | Frobenius loss with RKHS smoothness |
| `ssdr_p()` | SSDr-P | count matrix | kernel Poisson loss with RKHS smoothness |
| `ssdr_nn()` | SSDr-NN | count matrix | neural Poisson loss with graph smoothness |

`ssdr_p()` is retained as an exported supplementary/experimental method. The
main manuscript workflow emphasizes `ssdr_f()` and `ssdr_nn()`.

## Common Arguments

The public API uses the same names across methods where possible:

- `X`: spots or cells in rows and features in columns.
- `coords`: spatial coordinates aligned with `nrow(X)`.
- `rank`: target embedding dimension.
- `bandwidth`: Gaussian spatial bandwidth after coordinate normalization.
- `lambda`: smoothness penalty weight.
- `max_iter`: maximum optimizer iterations or training epochs.
- `tol`: convergence tolerance.

All methods return an `ssdr_fit` object. The embedding is always available as
`fit$U`.

## Minimal Example

```r
library(ssdr)

set.seed(1)
X <- matrix(rnorm(120), nrow = 30)
coords <- cbind(runif(30), runif(30))

fit <- ssdr_f(X, coords, rank = 2, bandwidth = 0.2, max_iter = 3)
fit
embedding <- fit$U
```

For count data:

```r
counts <- matrix(rpois(150, lambda = 3), nrow = 30)
fit_p <- ssdr_p(counts, coords, rank = 2, bandwidth = 0.2, max_iter = 1)
```
