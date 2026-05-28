# SSDr

[GitHub repository for reviewers](https://github.com/ChengyuRicardoDu/SSDr)

Spatially Smoothed Dimension Reduction (SSDr) is an R package for learning
spatially smooth low-dimensional embeddings from spatial transcriptomics data.
Supplementary analysis code and manuscript examples are maintained separately
from this package repository.
The package provides the two main methods used in the manuscript:

- `ssdr_f()` for continuous or preprocessed expression features.
- `ssdr_nn()` for raw UMI count matrices.

`ssdr_p()` is also exported for the kernel Poisson formulation described in the
Supplementary Methods.

## Installation

```r
install.packages("remotes")
remotes::install_github("ChengyuRicardoDu/SSDr")
```

`ssdr_nn()` requires the optional R `torch` package and a working torch backend:

```r
install.packages("torch")
torch::install_torch()
```

## Input Format

All SSDr functions use the same basic input layout:

- `X`: an `n x p` matrix with spots in rows and genes or features in columns.
- `coords`: an `n x d` matrix or data frame of spatial coordinates matched to
  the rows of `X`.
- `rank`: target embedding dimension.
- `bandwidth`: Gaussian spatial bandwidth after coordinate normalization.
- `lambda`: smoothness penalty weight.

Each function returns an `ssdr_fit` object. The learned embedding is stored in
`fit$U`.

## SSDr-F

`ssdr_f()` fits the kernel Frobenius formulation. It is intended for continuous
or preprocessed inputs, such as log-normalized expression features or PCA
scores.

```r
library(ssdr)

set.seed(1)
n <- 50
p <- 20

X <- matrix(rnorm(n * p), nrow = n)
coords <- cbind(
  x = runif(n),
  y = runif(n)
)

fit_f <- ssdr_f(
  X = X,
  coords = coords,
  rank = 2,
  bandwidth = 0.2,
  lambda = NULL,
  max_iter = 50
)

embedding_f <- fit_f$U
dim(embedding_f)
fit_f$parameters
```

When `lambda = NULL`, `ssdr_f()` uses the package default smoothness weight.
Spatial coordinates are globally min-max normalized inside the function before
the Gaussian kernel is constructed.

## SSDr-NN

`ssdr_nn()` fits the neural Poisson formulation. It should be used with a
non-negative count matrix on the original count scale.

```r
library(ssdr)

set.seed(1)
n <- 50
p <- 20

counts <- matrix(rpois(n * p, lambda = 3), nrow = n)
coords <- cbind(
  x = runif(n),
  y = runif(n)
)

fit_nn <- ssdr_nn(
  X = counts,
  coords = coords,
  rank = 2,
  bandwidth = 0.2,
  lambda = 0.01,
  hidden_dim = 32,
  hidden_layers = 2,
  learning_rate = 0.005,
  max_iter = 100,
  patience = 20,
  seed = 1
)

embedding_nn <- fit_nn$U
dim(embedding_nn)
fit_nn$diagnostics$loss_trace
```

For manuscript-scale analyses, `rank`, `bandwidth`, `lambda`, network size and
learning rate should be selected according to the analysis design. The example
above is intentionally small.

## Citation

If you use SSDr, please cite the accompanying manuscript and the software
repository:

```r
citation("ssdr")
```
