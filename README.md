# VAPOR

**V**ector **A**utoregressive **P**recision fact**OR** models.

VAPOR is a joint Bayesian model for the conditional independence structure of
high-dimensional time series. A VAR(1) removes the temporal dynamics, and a
low-rank-plus-diagonal (LRD) decomposition of the innovation precision matrix,
Omega = Lambda Lambda^T + Delta, estimates the contemporaneous conditional
dependence that remains. The VAR coefficients B and Omega are updated in a
single Gibbs sampler, so uncertainty in each informs the other. A sequential
approach, which fits the VAR first and applies a precision estimator to the
residuals, treats the fitted residuals as data and ignores the estimation
uncertainty in the VAR step.

The LRD form does not impose exact zeros in Omega. It concentrates the
off-diagonal mass into a small number of directions, shrinking most partial
correlations toward zero, and the reported graph is obtained by a separate
posterior selection step (see [Selecting a graph](#selecting-a-graph)).

When the number of series approaches or exceeds the number of time points,
VAPOR replaces B with a reduced-rank factorization B = U V^T of fixed rank r.

This repository holds the R package and the code that reproduces every table
and figure in the paper and its supplements.

## Installation

```r
# install.packages("remotes")
remotes::install_github("peterdunson/vapor")
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
partial_correlations(fit)  # posterior mean partial correlations
```

### Full-rank and reduced-rank dynamics

By default `rank = "full"`, which places a Dirichlet-Laplace prior on each row
of B. Setting `rank` to a small integer uses the reduced-rank parameterization
B = U V^T:

```r
fit_rr <- vapor(Y, lags = 1, k = 2, rank = 3,
                prior_B = "dl", nrun = 5000, burn = 2500, thin = 5)
```

Under reduced rank, U has a ridge prior and `prior_B` sets the prior on V:
`"dl"` (default) for the Dirichlet-Laplace prior, or `"ridge"` for a
constant-variance Gaussian prior. `prior_B` is ignored when `rank = "full"`.

The paper uses r = 3 in both applications and in Simulations E and F, and the
full-rank parameterization in Simulations A through D.

## Model

For a VAR(1):

- y_t = B y_{t-1} + e_t, e_t ~ N(0, Omega^{-1})
- Omega = Lambda Lambda^T + Delta, with Lambda a p x k loadings matrix and
  Delta = diag(delta_1^2, ..., delta_p^2) the idiosyncratic precisions

**Priors.**

- **VAR coefficients B (full rank).** Dirichlet-Laplace prior on each row of
  B, following the equation-by-equation scheme of Kastner and Huber (2020),
  with shape a_B = 1/(p K_x), the inverse of the number of VAR coefficients,
  where K_x = p * lags.
- **Reduced rank.** Ridge prior on U. On V, either a ridge prior or the same
  Dirichlet-Laplace prior applied row by row over the K_x rows of V, with shape
  a_V = 1/(K_x r).
- **Loadings Lambda.** Dirichlet-Laplace prior following Chandra, Müller, and
  Sarkar (2026), with local shape a_L = 0.5.
- **Idiosyncratic precisions.** Dirichlet process prior on delta_1^2, ...,
  delta_p^2. The base measure is Gamma(1, 1) (`as = bs = 1`), which has mean
  one on the scale of standardized data. Chandra, Müller, and Sarkar (2026)
  use Gamma(0.1, 0.1).
- **DP concentration.** Gamma(0.1, 0.1) (`a_dir = b_dir = 0.1`), updated by
  the auxiliary-variable scheme of Escobar and West (1995)
  (`dp_prior = "escobar_west"`, the default). The conjugate Stirling-gamma
  prior of Zito et al. (2026) is available as `dp_prior = "stirling_gamma"`.

**Sampler.** Each iteration cycles through three blocks: sample B given Omega;
sample Lambda and Delta given the VAR residuals e_t = y_t - B y_{t-1}; and
update the DP concentration and cluster allocations. Under reduced rank the
first block samples U given V, then V given U.

**Numerical safeguard.** After each draw, every idiosyncratic precision is set
to max(delta_j^2, 1e-8), as a safeguard for the Cholesky factorization in the
update of Lambda. The fit records how often the bound was reached:
`fit$floor_hits` per series, `fit$floor_first` for the first such iteration,
and `fit$floor_iters` for the number of iterations in which any series was
bounded. The bound is not reached in any of the paper's application fits.

```r
fit$floor_iters   # 0 if the bound was never reached
```

## Selecting a graph

Omega is dense, so a conditional independence graph is selected from the
posterior. The paper uses two rules.

**Credible interval (the application).** An edge is reported when the central
95% posterior credible interval of its partial correlation excludes zero. The
package has no function for this rule; it takes a few lines on a fit's draws:

```r
# per-draw partial correlations
pc_draws <- function(fit) {
   d <- dim(fit$Lambda); p <- d[1]; n <- d[3]
   out <- array(0, dim = c(p, p, n))
   for (s in seq_len(n)) {
      om <- tcrossprod(fit$Lambda[, , s])
      diag(om) <- diag(om) + fit$residual_prec[s, ]
      sd <- sqrt(diag(om))
      out[, , s] <- -om / outer(sd, sd)
   }
   out
}

PC <- pc_draws(fit)
lo <- apply(PC, c(1, 2), quantile, 0.025)
hi <- apply(PC, c(1, 2), quantile, 0.975)
P  <- (lo > 0 | hi < 0) * 1
diag(P) <- 0

sum(P) / 2   # number of edges
```

In the paper the draws of four chains are pooled before the quantiles are
taken, and the reported partial correlations are the posterior means of the
per-draw partial correlations over the same pooled draws.

**pFDR (the simulations).** `select_graph()` implements the posterior
false-discovery-rate procedure of Chandra, Müller, and Sarkar (2026), which
controls the positive false discovery rate of Storey (2003). The simulations
use `alpha = 0.1` and `beta = 0.85`:

```r
sel <- select_graph(fit, alpha = 0.1, beta = 0.85)

sum(sel$P) / 2   # number of edges
sel$P            # adjacency matrix, zero on the diagonal
sel$pFDR         # realized pFDR
```

The function defaults are `alpha = 0.1`, `beta = 0.7`, and
`epsilon_grid = seq(0.05, 0.3, length.out = 15)`.

## Main functions

- `vapor()` fits the model and returns a `"vapor"` object.
- `select_graph()` performs pFDR edge selection on a fit.
- `partial_correlations()` returns the posterior mean partial correlation
  matrix, with `NA` on the diagonal.
- `posterior_mean_omega()` builds the posterior mean precision matrix from the
  stored draws (`Lambda` and `residual_prec`).
- `vapor_warm_start()` builds starting values from a least-squares VAR fit; it
  is used by `vapor()` when `init = "warm"`.

## Arguments

`vapor()` takes a stationary `T x p` matrix `Y` and the following controls:

- `lags`: VAR lag order (default 1)
- `k`: number of precision factors (default 5)
- `rank`: `"full"` (default) or an integer `r <= min(p, lags * p)` for the
  reduced-rank parameterization
- `prior_B`: prior on V under reduced rank, `"dl"` (default) or `"ridge"`
- `nrun`, `burn`, `thin`: defaults 5000, `floor(nrun / 2)`, and 5. Draws are
  thinned first and the burn-in is discarded from the thinned chain, so the
  defaults keep 500 draws
- `as`, `bs`: shape and rate of the DP base measure, both 1
- `a_dir`, `b_dir`: shape and rate of the prior on the DP concentration, both
  0.1
- `dp_prior`: `"escobar_west"` (default) or `"stirling_gamma"`
- `init`: `"random"` (default) or `"warm"`; `init_perturb` (default 0.15)
  scales the perturbation of the warm start
- `jitter`: Gaussian jitter added to `Y`. If `NULL` (default) it is `1e-3`
  when `T/p < 5` and `1e-6` otherwise
- `seed`: optional seed
- `a_B`, `b_B`: Dirichlet-Laplace shape and rate for B under full rank,
  defaults `1 / (p * K_x)` and 0.5
- `a_V`, `b_V`: Dirichlet-Laplace shape and rate for V under reduced rank,
  defaults `1 / (K_x * r)` and 0.5
- `s_U2`: ridge prior variance for U (default 1)
- `s_V2`: ridge prior variance for V when `prior_B = "ridge"` (default 1)
- `a_L`, `b_L`: Dirichlet-Laplace shape and rate for Lambda, both 0.5

## Reproducing the paper

```
applications/   the FRED-QD application at p = 29 and p = 100
simulation/     Simulations A to F
```

**Data.** The FRED-QD database of McCracken and Ng (2020), May 2026 vintage,
included at `applications/data/2026-05-qd.csv`. Source:
<https://www.stlouisfed.org/research/economists/mccracken/fred-databases/>.
Each series is transformed with the transformation code distributed with the
database. Both panels train on 1996Q1 to 2020Q4 (100 quarters) and test on
2021Q1 to 2025Q4 (20 quarters), and each series is standardized with its
training mean and standard deviation. Because source agencies revise
historical values, a different vintage gives different estimates.
`applications/fredqd_groups.csv` maps each series to its FRED-QD group, and
`applications/varlist_descriptions.csv` supplies the descriptions in the
supplement's variable list.

**Settings.**

| | p = 29 | p = 100 |
| --- | --- | --- |
| rank r | 3 | 3 |
| factors k | 10 | 25 |
| chains (seeds) | 4 (101 to 104) | 4 (101 to 104) |
| iterations, burn-in, thinning | 100,000, 20,000, 80 | 100,000, 20,000, 80 |
| series held out per test quarter | 8 | 30 |
| held-out sets | 100 | 100 |

Kastner's model at p = 29 uses seeds 101, 102, 104 and 105, because seed 103
stopped with a Cholesky decomposition error. At p = 100 all eight
100,000-iteration Kastner chains (seeds 101 to 108) stopped with Cholesky
errors, so it runs four chains of 20,000 iterations with a burn-in of 10,000
and thinning of 10.

The simulations run 25 replicates per scenario. The Gibbs samplers for
PrecFact, VAR + PrecFact and VAPOR run 5000 iterations with a burn-in of 2500
and thinning of 5, with k = 5.

**Running it.** `applications/run_all.sh` runs the application and the
simulations in dependency order:

```sh
bash applications/run_all.sh
```

Each stage runs as its own `Rscript` process, and a failure stops the stages
that depend on it. The simulations can also be run on their own:

```sh
Rscript simulation/final_simulations.R          # Simulations A to E
Rscript simulation/run_sim_f_setting.R 50  35   # Simulation F, one setting
Rscript simulation/run_sim_f_setting.R 100 70
Rscript simulation/run_sim_f_setting.R 150 105
Rscript simulation/combine_sim_f.R              # combine the three settings
Rscript simulation/final_sim_tables.R           # supplement tables
```

**Competitor inputs.** Four objects under `applications/results/fits/` hold the
competitor fits that `run_all.sh` reads: `reused_p29.rds` and `reused_p100.rds`
(Kastner, VAR + Glasso, the sample covariance and the random walk, with the
panel they were fitted on) and `kastner_p100_20k.rds` and
`kastner_p100_20k_draws.rds` (Kastner at p = 100). `applications/rebuild_reused.R`
and `applications/rebuild_kastner_p100_20k.R` regenerate them:

```sh
REBUILD_INPUTS=1 bash applications/run_all.sh
```

At p = 29 the rebuild reproduces the stored object exactly.

**Outputs.**

| | |
| --- | --- |
| `applications/results/figures/` | the eight application figures |
| `applications/results/tables/` | `varlist_p100.tex`, `results_summary.md` |
| `applications/results/tables/paper_checks.txt` | checks for numbers stated in the text |
| `applications/results/tables/ppc_comparison_p29.rds` and `_p100.rds` | posterior predictive check for the competing models |
| `simulation/results/figures/` | `sim_a_heatmap_rho95_T100.png`, `sim_e_heatmap_p150_T150.png` |
| `simulation/results/tables/` | `sim_tables.tex` |
| `*/results/fits/`, `*/results/logs/` | per-stage `.rds` and logs |

`results_summary.md` collects every number reported in the paper.

**Run times** on an Apple M2 laptop:

| stage | time |
| --- | --- |
| `final_p29.R` | 16 chains, 0.4 to 1.1 min each |
| `final_p100.R` | 16 chains, 3.4 to 14.0 min each |
| `final_sweeps.R 29` | 9 min |
| `final_sweeps.R 100` | 130 min |
| Simulation F, p = 50, 100, 150 | 27, 74, 291 min |

The three Simulation F settings are independent and can run in parallel.

**Software.** R 4.5.1 on aarch64-apple-darwin20.

| package | version | | package | version |
| --- | --- | --- | --- | --- |
| VAPOR | 0.1.0 | | ggplot2 | 4.0.3 |
| PrecFact | 1.0 | | igraph | 2.2.1 |
| bayesianVARs | 0.1.8 | | pheatmap | 1.0.13 |
| factorstochvol | 1.1.2 | | posterior | 1.6.1 |
| glasso | 1.11 | | scales | 1.4.0 |
| MASS | 7.3-65 | | gridExtra | 2.3 |
| mvtnorm | 1.3-3 | | Rcpp | 1.1.2 |
| tidyverse | 2.0.0 | | RcppArmadillo | 15.4.0.1 |

All but VAPOR and PrecFact are on CRAN. PrecFact links against GSL
(macOS: `brew install gsl`). The version used is the fork at
`peterdunson/PrecFact`; the authors' repository is `noirritchandra/PrecFact`.

```r
remotes::install_github("peterdunson/PrecFact", ref = "72059b0")
```

## Citation

Dunson, P. (2026). Bayesian precision factor models for high-dimensional
macroeconomic time series. Working paper.

## References

- Chandra, N. K., Müller, P., and Sarkar, A. (2026). Bayesian scalable
  precision factor analysis for Gaussian graphical models. *Bayesian Analysis*,
  21(1), 75-103.
- Escobar, M. D. and West, M. (1995). Bayesian density estimation and
  inference using mixtures. *Journal of the American Statistical Association*,
  90(430), 577-588.
- Friedman, J., Hastie, T., and Tibshirani, R. (2008). Sparse inverse
  covariance estimation with the graphical lasso. *Biostatistics*, 9(3),
  432-441.
- Gruber, L., Haan, S., and Kastner, G. (2026). *bayesianVARs: MCMC Estimation
  of Bayesian Vectorautoregressions*. R package.
  <https://doi.org/10.32614/CRAN.package.bayesianVARs>
- Kastner, G. and Huber, F. (2020). Sparse Bayesian vector autoregressions in
  huge dimensions. *Journal of Forecasting*, 39(7), 1142-1165.
- McCracken, M. and Ng, S. (2020). FRED-QD: A quarterly database for
  macroeconomic research. NBER Working Paper 26872.
- Storey, J. D. (2003). The positive false discovery rate: a Bayesian
  interpretation and the q-value. *Annals of Statistics*, 31(6), 2013-2035.
- Zito, A., Rigon, T., and Dunson, D. B. (2026). Bayesian nonparametric
  modeling of latent partitions via Stirling-gamma priors. *Bayesian Analysis*.

## License

MIT
