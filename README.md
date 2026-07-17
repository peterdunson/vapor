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

When the number of series is close to the number of time points (T ~ p), the
full VAR coefficient matrix is over-parameterized. VAPOR therefore also provides
a **reduced-rank** parameterization of the dynamics, B = U V^T with a small rank
r, which regularizes the dynamics and is the recommended choice in that regime.

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

### Full-rank vs reduced-rank dynamics

By default `rank = "full"`, which uses the full-rank sampler with a
Dirichlet-Laplace prior on the VAR coefficients B. When T is close to p, set
`rank` to a small integer to use the reduced-rank parameterization B = U V^T,
which shrinks the dynamics and avoids overfitting:

```r
# reduced-rank dynamics (recommended when T is close to p)
fit_rr <- vapor(Y, lags = 1, k = 2, rank = 3,
                prior_B = "dl", nrun = 5000, burn = 2500, thin = 5)
```

`prior_B` controls the prior on the reduced-rank factor V: `"dl"` (default)
applies Dirichlet-Laplace shrinkage, `"ridge"` a constant-variance Gaussian
prior, and `"mgps"` a Multiplicative Gamma Process prior whose shrinkage
increases with column index -- intended for use with a generous rank ceiling,
letting the effective rank emerge from the posterior. `prior_B` is ignored when
`rank = "full"`.

### Selecting a sparse graph

`select_graph()` thresholds the precision matrix into a conditional independence
graph using posterior false discovery rate (pFDR) control. It returns a 0/1
adjacency matrix of selected edges.

```r
sel <- select_graph(fit, alpha = 0.1, beta = 0.85)

sum(sel$P) / 2   # number of edges
sel$P            # adjacency matrix
sel$pFDR         # realized false discovery rate
```

For each pair of variables, the procedure asks how often across the posterior
draws the absolute partial correlation exceeds a magnitude threshold. Pairs
exceeding the probability cutoff `beta` are declared edges, subject to keeping
the expected proportion of false edges below `alpha`.

## Model

For a VAR(1):

- y_t = B y_{t-1} + e_t,  e_t ~ N(0, Omega^{-1})
- Omega = Lambda Lambda^T + Delta  (low-rank plus diagonal)
- B is either full rank, with a Dirichlet-Laplace prior on its rows
  (Kastner and Huber, 2020), or reduced rank, B = U V^T, with a
  Dirichlet-Laplace, ridge, or Multiplicative Gamma Process prior on V
- Dirichlet-Laplace prior on the loadings Lambda (Chandra, Müller, and Sarkar, 2026)
- Dirichlet process prior on the diagonal Delta

The sampler iterates three blocks: sample B given Omega, sample Omega given the
VAR residuals, and update the Dirichlet process concentration. Under the
reduced-rank parameterization the first block samples U given V, then V given U,
with the residuals e_t = y_t - U V^T y_{t-1} passed to the remaining blocks.

## Main functions

- `vapor()` fits the model and returns a `"vapor"` object.
- `select_graph()` performs pFDR edge selection on a fit.
- `partial_correlations()` returns the partial correlation matrix from a fit.
- `posterior_mean_omega()` builds the posterior mean precision matrix from
  stored draws.

## Arguments

`vapor()` takes a stationary `T x p` matrix `Y` and the following key controls:

- `lags`: VAR lag order (default 1)
- `k`: number of precision factors (default 5)
- `rank`: rank of the VAR coefficient matrix B. `"full"` (default) uses the
  full-rank sampler; an integer `r <= min(p, lags*p)` uses the reduced-rank
  parameterization B = U V^T
- `prior_B`: prior on the reduced-rank factor V (`"dl"`, `"ridge"`, or
  `"mgps"`); ignored when `rank = "full"`
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
