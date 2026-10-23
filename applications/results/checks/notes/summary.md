# Does the 1e-8 residual-precision floor bind in the simulation fits?

Replicate 1 of the first (smallest) setting of each of `sim_a.R` ... `sim_g.R`,
fitted twice with `vapor()`: once at the package defaults `as = bs = 0.1` and
once at `as = bs = 1`. Everything else — the DGP, the seed, `k`, `rank`,
`prior_B`, `nrun`, `burn`, `thin`, `init` — is the script's own.

## Method

`sim_*.R` are neither modified nor run. For each script the harness
(`applications/sim_base_measure.R`) parses the file and evaluates only the
top-level expressions **before its first top-level `for` loop** — the
libraries, the constants and the DGP and helper definitions, 22 to 35
expressions per script — then calls that script's own DGP for replicate 1 and
fits `vapor()` with the script's own argument list, changing only `as` and
`bs`. One chain per fit, which is what the simulation scripts themselves run,
so the reported R-hat is the **within-chain split statistic**, not a
between-chain one.

Metric is the one the script reports for VAPOR: `||cov2cor(Omega_hat) -
cov2cor(Omega_true)||_F` for sim_a, sim_b, sim_d, sim_e, sim_f and sim_g, and
held-out MSPE for sim_c, which has no single true Omega.

## Result

**The floor never binds, at either base measure, in any of the fourteen fits.**
`at_floor`, `< 1e-6`, `iterations with any series floored` and
`series ever floored` are 0 in every row. The smallest idiosyncratic precision
seen anywhere is 9.724e-04 (sim_b at the default), roughly five orders of
magnitude above the 1e-8 floor.

| script | setting | base measure | min residual_prec | at floor | < 1e-6 | iters any floor | series floored | max R-hat | median R-hat | min bulk ESS | metric |
|---|---|---|---|---|---|---|---|---|---|---|---|
| sim_a | rho=0.10, p=20, T=200 | 0.1 (default) | 7.595e-01 | 0 | 0 | 0 | 0 | 1.019 | 1.001 | 114 | 0.6512 |
| sim_a | rho=0.10, p=20, T=200 | 1 | 7.908e-01 | 0 | 0 | 0 | 0 | 1.020 | 1.001 | 120 | 0.6530 |
| sim_b | rho=0.10, p=20, T=200 | 0.1 (default) | 9.724e-04 | 0 | 0 | 0 | 0 | 1.832 | 1.010 | 1 | 0.8400 |
| sim_b | rho=0.10, p=20, T=200 | 1 | 1.066e-02 | 0 | 0 | 0 | 0 | 1.093 | 1.002 | 16 | 0.7799 |
| sim_c | sparsity=0.01, p=20, T=200 | 0.1 (default) | 1.118e-03 | 0 | 0 | 0 | 0 | 1.332 | 1.005 | 2 | 0.8391 (MSPE) |
| sim_c | sparsity=0.01, p=20, T=200 | 1 | 4.999e-03 | 0 | 0 | 0 | 0 | 1.303 | 1.001 | 3 | 0.8409 (MSPE) |
| sim_d | val=0.5, p=20, T=200 | 0.1 (default) | 8.292e-01 | 0 | 0 | 0 | 0 | 1.023 | 1.001 | 82 | 0.4902 |
| sim_d | val=0.5, p=20, T=200 | 1 | 8.281e-01 | 0 | 0 | 0 | 0 | 1.020 | 1.001 | 131 | 0.4831 |
| sim_e | p=10, k=2, T=200 | 0.1 (default) | 7.111e-01 | 0 | 0 | 0 | 0 | 1.013 | 1.000 | 234 | 0.1946 |
| sim_e | p=10, k=2, T=200 | 1 | 6.261e-01 | 0 | 0 | 0 | 0 | 1.011 | 1.000 | 183 | 0.1929 |
| sim_f | p=50, T=50 | 0.1 (default) | 9.035e-01 | 0 | 0 | 0 | 0 | 1.049 | 1.002 | 44 | 4.2896 |
| sim_f | p=50, T=50 | 1 | 3.799e-01 | 0 | 0 | 0 | 0 | 1.094 | 1.002 | 10 | 4.2695 |
| sim_g | p=50, T=35 | 0.1 (default) | 7.479e-01 | 0 | 0 | 0 | 0 | 1.061 | 1.001 | 21 | 4.4147 |
| sim_g | p=50, T=35 | 1 | 5.920e-01 | 0 | 0 | 0 | 0 | 1.038 | 1.001 | 51 | 4.3963 |

## Where the two base measures differ

Accuracy is unchanged to within 0.07 in every script; the largest move is
sim_b, 0.8400 to 0.7799.

Mixing differs in the two misspecified DGPs. sim_b (sparse Omega, no factor
structure) goes from max R-hat 1.832 and min bulk ESS 1 at the default to
1.093 and 16; sim_c (factor stochastic volatility) from 1.332 / 2 to 1.303 / 3.
sim_f moves the other way, 1.049 / 44 to 1.094 / 10. The remaining four are
within 0.01 of each other on max R-hat.

## Caveat

One replicate and one chain per script. These are the smallest settings in each
script; the application fits that do hit the floor are p = 29 and p = 100 with
T = 100, which is a smaller T/p than any of the settings here except sim_f and
sim_g.

Scope: `applications/sim_base_measure.R`, written for this check. Nothing in
`simulation/`, the package defaults or any existing result was modified.
