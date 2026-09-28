# SSDr

Spatially Smoothed Dimension Reduction (SSDr) estimates spatially smooth embeddings.

[GitHub repository](https://github.com/ChengyuRicardoDu/SSDr)

## Installation

```r
install.packages("remotes")
remotes::install_github("ChengyuRicardoDu/SSDr")
```

## SSDr-F

```r
X <- matrix(rnorm(200), 20, 10)
coords <- cbind(runif(20), runif(20))

fit_f <- ssdr_f(X, coords, rank = 2, bandwidth = 0.2, max_iter = 50)
```

## SSDr-NN

```r
counts <- matrix(rpois(200, lambda = 3), 20, 10)
coords <- cbind(runif(20), runif(20))

fit_nn <- ssdr_nn(counts, coords, rank = 2, max_iter = 100)
```

## Joint SSDr-F

```r
Y <- list(matrix(rnorm(200), 20, 10),
          matrix(rnorm(240), 24, 10))
coords_joint <- list(cbind(runif(20), runif(20)),
                     cbind(runif(24), runif(24)))

fit_joint <- ssdr_joint_f(Y, coords_joint, rank = 2, bandwidth = 0.2,
                          lambda = 0.02, solver = "exact")
```

Example analyses: [SSDr-analysis](https://github.com/ChengyuRicardoDu/SSDr-analysis).
