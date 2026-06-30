# VAPOR

**V**ector **A**utoregressive **P**recision fact**OR** models.

VAPOR jointly estimates vector autoregression (VAR) dynamics and a sparse
innovation precision matrix in a single Bayesian Gibbs sampler. The innovation
precision matrix uses a low-rank plus diagonal factor decomposition,
Omega = Lambda Lambda^T + Delta, so its off-diagonal structure encodes the
contemporaneous conditional dependence of the innovations after temporal
dynamics are removed.

The motivation: a sequential approach (fit a VAR, then estimate a precision
matrix on the residuals) ignores estimation uncertainty in the VAR
coefficients. VAPOR updates the VAR coefficients and the precision matrix
jointly, propagating uncertainty between them. This helps most when the VAR is
hard to estimate: strong autocorrelation, small samples, or high dimension.

## Installation

```r
# install.packages("remotes")
remotes::install_github("peterdunson/VAPOR")
```

The package depends only on `Rcpp`, `RcppArmadillo`, and `RcppProgress`. No GSL
or OpenMP system libraries are required.

## Usage

```r
library(VAPOR)

# Y: T x p matrix of stationary series (rows = time points)
set.seed(1)
Y <- matrix(rnorm(200 * 5), 200, 5)

fit <- vapor(Y, lags = 1, k = 2, nrun = 5000, burn = 2500, thin = 5)

fit                        # summary
fit$Omega                  # posterior mean precision matrix
partial_correlations(fit)  # partial correlation matrix
```

## Model

For a VAR(1):

- y_t = B y_{t-1} + e_t,  e_t ~ N(0, Omega^{-1})
- Omega = Lambda Lambda^T + Delta  (low-rank plus diagonal)
- Dirichlet-Laplace prior on the rows of B (Kastner and Huber, 2020)
- Dirichlet-Laplace prior on the loadings Lambda (Chandra and others, 2026)
- Dirichlet process prior on the diagonal Delta

The sampler iterates three blocks: sample B given Omega, sample Omega given the
VAR residuals, and update the Dirichlet process concentration.

## Main functions

- `vapor()` fits the model and returns a `"vapor"` object.
- `partial_correlations()` returns the partial correlation matrix from a fit.
- `posterior_mean_omega()` builds the posterior mean precision matrix from
  stored draws.

## Arguments

`vapor()` takes a stationary `T x p` matrix `Y` and the following key controls:

- `lags`: VAR lag order (default 1)
- `k`: number of precision factors (default 5)
- `nrun`, `burn`, `thin`: MCMC controls
- `jitter`: optional Gaussian jitter for numerical stability; chosen
  automatically if `NULL`
- `seed`: optional seed

## References

- Chandra, N. K., Müller, P., and Sarkar, A. (2026). Bayesian Scalable
  Precision Factor Analysis for Gaussian Graphical Models. *Bayesian Analysis*,
  21(1), 75-103.
- Kastner, G. and Huber, F. (2020). Sparse Bayesian vector autoregressions in
  huge dimensions. *Journal of Forecasting*, 39(7), 1142-1165.

## License

MIT
