#' Warm-start initial values for a VAPOR fit
#'
#' Builds starting values for Lambda, Delta, and (reduced-rank only) V from a
#' fast ridge-VAR point estimate and a shrunk residual covariance. Values are
#' perturbed so that chains started from different seeds remain dispersed,
#' keeping between-chain convergence diagnostics meaningful.
#'
#' @param Y Numeric T x p matrix (already jittered, as passed to the sampler).
#' @param lags,k,rank As in \code{vapor}.
#' @param ridge Ridge penalty for the VAR point estimate.
#' @param shrink Shrinkage of the residual covariance toward a scaled identity.
#' @param perturb Relative SD of the multiplicative perturbation.
#' @return A list with elements Lambda, delta, and V (V is NULL for full rank).
#' @export
vapor_warm_start <- function(Y, lags = 1, k = 5, rank = "full",
                             ridge = 1e-2, shrink = 0.2, perturb = 0.15) {
   Y <- as.matrix(Y); T_obs <- nrow(Y); p <- ncol(Y); Teff <- T_obs - lags
   
   Ymat <- Y[(lags + 1):T_obs, , drop = FALSE]
   X <- do.call(cbind, lapply(seq_len(lags),
                              function(l) Y[(lags + 1 - l):(T_obs - l), , drop = FALSE]))
   Kx <- ncol(X)
   
   # ridge VAR point estimate
   Bhat <- t(solve(crossprod(X) + ridge * Teff * diag(Kx), crossprod(X, Ymat)))
   E <- Ymat - X %*% t(Bhat)
   
   # shrunk residual covariance, inverted
   S  <- crossprod(E) / Teff
   S  <- (1 - shrink) * S + shrink * mean(diag(S)) * diag(p)
   Om <- solve(S)
   Om <- (Om + t(Om)) / 2
   
   # LRD split: Delta = half the diagonal, Lambda from top-k eigenpairs
   delta <- 0.5 * diag(Om)
   R  <- Om - diag(delta)
   ev <- eigen(R, symmetric = TRUE)
   kk <- min(k, sum(ev$values > 0))
   Lambda <- matrix(0, p, k)
   if (kk > 0)
      Lambda[, seq_len(kk)] <- ev$vectors[, seq_len(kk), drop = FALSE] %*%
      diag(sqrt(ev$values[seq_len(kk)]), kk)
   
   # reduced-rank V from the SVD of Bhat
   Vinit <- NULL
   if (!identical(rank, "full")) {
      r  <- as.integer(rank)
      sv <- svd(Bhat, nu = r, nv = r)
      Vinit <- sv$v %*% diag(sqrt(sv$d[seq_len(r)]), r)
   }
   
   # perturb so chains from different seeds start dispersed
   pert <- function(M) M * (1 + perturb * matrix(stats::rnorm(length(M)),
                                                 nrow(M), ncol(M)))
   Lambda <- pert(Lambda)
   if (!is.null(Vinit)) Vinit <- pert(Vinit)
   delta  <- pmax(delta * (1 + perturb * stats::rnorm(p)), 1e-8)
   
   list(Lambda = Lambda, delta = delta, V = Vinit)
}

#' Fit a VAPOR model
#'
#' Joint Bayesian estimation of VAR dynamics and a sparse innovation precision
#' matrix.
#'
#' @param Y A numeric T x p matrix of stationary series (rows = time).
#' @param lags Integer VAR lag order.
#' @param k Integer number of precision factors.
#' @param rank Rank of the VAR coefficient matrix B. `"full"` (default) uses
#'   the original full-rank sampler with a Dirichlet-Laplace prior on B. An
#'   integer `r <= min(p, lags*p)` uses a reduced-rank parameterization
#'   B = U V^T, with `prior_B` controlling the prior on V.
#' @param nrun,burn,thin MCMC controls.
#' @param a_L,b_L Dirichlet-Laplace shape/rate for the precision loadings.
#' @param a_B,b_B Dirichlet-Laplace shape and fixed global rate for the
#'   full-rank VAR coefficient matrix B. Used only when rank = "full". The
#'   default `a_B = 1 / (ncol(Y) * ncol(Y) * lags)` is the reciprocal of the
#'   number of VAR coefficients, p * K_x with K_x = p * lags, as in the paper.
#' @param a_V,b_V Dirichlet-Laplace shape and fixed global rate for the
#'   reduced-rank factor V. Used only when rank is numeric and
#'   prior_B = "dl". The default `a_V = 1 / (ncol(Y) * lags * rank)` is the
#'   reciprocal of the number of entries of V, K_x * r. Its default mentions
#'   `rank`, so it is evaluated only when `rank` is an integer.
#' @param prior_B Prior on the reduced-rank factor V. Ignored when
#'   rank = "full" (always Dirichlet-Laplace on the full B). `"dl"`
#'   (default) uses Dirichlet-Laplace shrinkage on V; `"ridge"` uses a
#'   constant-variance Gaussian prior.
#' @param s_U2,s_V2 Ridge prior variances for U and V. Only used when rank
#'   is numeric; s_V2 only applies when prior_B = "ridge".
#' @param as,bs Dirichlet process base-measure shape/rate.
#' @param a_dir,b_dir Dirichlet process concentration prior shape/rate.
#' @param dp_prior Update for the DP concentration parameter.
#'   `"escobar_west"` (default) uses the auxiliary-variable scheme;
#'   `"stirling_gamma"` uses the conjugate Stirling-gamma prior.
#' @param sg_a,sg_b,sg_m Stirling-gamma hyperparameters. Only used when
#'   dp_prior = "stirling_gamma"; sg_m = 0 defaults to p.
#' @param init Initialization scheme. `"warm"` (default) starts Lambda,
#'   Delta, and V from a fast ridge-VAR point estimate, perturbed per seed;
#'   `"random"` uses the sampler's internal random initialization.
#' @param init_perturb Relative SD of the warm-start perturbation. Larger
#'   values disperse chains further, keeping Gelman-Rubin diagnostics honest.
#' @param jitter SD of Gaussian jitter for stability; NULL picks automatically.
#' @param seed Optional integer seed.
#' @return An object of class "vapor".
#' @export
vapor <- function(Y, lags = 1, k = 5, rank = "full",
                  nrun = 5000, burn = floor(nrun / 2), thin = 5,
                  a_L = 0.5, b_L = 0.5,
                  a_B = 1 / (ncol(Y) * ncol(Y) * lags), b_B = 0.5,
                  a_V = 1 / (ncol(Y) * lags * rank),    b_V = 0.5,
                  prior_B = c("dl", "ridge"),
                  s_U2 = 1.0, s_V2 = 1.0,
                  as = 1.0, bs = 1.0, a_dir = 0.1, b_dir = 0.1,
                  dp_prior = c("escobar_west", "stirling_gamma"),
                  sg_a = 3.0, sg_b = 1.0, sg_m = 0,
                  init = c("random", "warm"), init_perturb = 0.15,
                  jitter = NULL, seed = NULL) {
   
   Y <- as.matrix(Y)
   if (!is.numeric(Y))      stop("`Y` must be numeric.")
   if (anyNA(Y))            stop("`Y` contains missing values.")
   if (burn >= nrun)        stop("`burn` must be less than `nrun`.")
   
   prior_B  <- match.arg(prior_B)
   dp_prior <- match.arg(dp_prior)
   init     <- match.arg(init)

   #  every DL hyperparameter must be one positive finite number.  a_V is
   #  checked in the reduced-rank branch instead: its default refers to `rank`
   #  and forcing it here would error whenever rank = "full".
   chk_pos <- function(x, nm){
      if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
         stop(sprintf("`%s` must be a single positive finite number.", nm))
      invisible(TRUE) }
   for (nm in c("a_L", "b_L", "b_B", "b_V", "s_U2", "s_V2",
                "as", "bs", "a_dir", "b_dir"))
      chk_pos(get(nm), nm)
   
   p <- ncol(Y); T_obs <- nrow(Y)
   if (!is.null(seed)) set.seed(seed)
   
   if (is.null(jitter)) jitter <- if (T_obs / p < 5) 1e-3 else 1e-6
   if (jitter > 0) Y <- Y + matrix(stats::rnorm(T_obs * p, 0, jitter), T_obs, p)
   
   ws <- if (init == "warm")
      vapor_warm_start(Y, lags = lags, k = k, rank = rank, perturb = init_perturb)
   else list(Lambda = NULL, delta = NULL, V = NULL)
   
   if (identical(rank, "full")) {
      
      chk_pos(a_B, "a_B")
      res <- var_precision_sampler(
         Y = Y, lags = as.integer(lags), k = as.integer(k),
         nrun = as.integer(nrun), burn = as.integer(burn), thin = as.integer(thin),
         a_L = a_L, b_L = b_L, a_B = a_B, b_B = b_B,
         as = as, bs = bs, a_dir = a_dir, b_dir = b_dir,
         dp_prior = dp_prior, sg_a = sg_a, sg_b = sg_b, sg_m = sg_m,
         Lambda_init = ws$Lambda, delta_init = ws$delta
      )
      rank_used <- "full"
      
   } else {
      
      if (!is.numeric(rank) || rank < 1)
         stop("`rank` must be \"full\" or a positive integer.")
      max_rank <- min(p, lags * p)
      if (rank > max_rank)
         stop(sprintf("`rank` must be <= min(p, lags*p) = %d.", max_rank))
      chk_pos(a_V, "a_V")
      
      res <- var_precision_sampler_rr(
         Y = Y, lags = as.integer(lags), k = as.integer(k), rank = as.integer(rank),
         nrun = as.integer(nrun), burn = as.integer(burn), thin = as.integer(thin),
         a_L = a_L, b_L = b_L, prior_B = prior_B,
         s_U2 = s_U2, s_V2 = s_V2, a_V = a_V, b_V = b_V,
         as = as, bs = bs, a_dir = a_dir, b_dir = b_dir,
         dp_prior = dp_prior, sg_a = sg_a, sg_b = sg_b, sg_m = sg_m,
         Lambda_init = ws$Lambda, delta_init = ws$delta, V_init = ws$V
      )
      rank_used <- as.integer(rank)
   }
   if (is.null(res)) stop("Sampler aborted.")
   
   n_mc <- dim(res$Lambda)[3]
   keep <- seq.int(floor(burn / thin) + 1L, n_mc)
   
   Lambda   <- res$Lambda[, , keep, drop = FALSE]
   ps       <- res$residuals[keep, , drop = FALSE]
   B        <- res$B[, , keep, drop = FALSE]
   loglik   <- res$loglik[, keep, drop = FALSE]
   dir_prec <- res$dir_prec[keep]
   nclust   <- res$nclust[keep]
   
   V_kept   <- NULL
   if (!identical(rank, "full")) {
      V_kept   <- res$V[, , keep, drop = FALSE]
   }
   
   Omega <- posterior_mean_omega(Lambda, ps)
   
   #  per-series count of the iterations at which the residual-precision floor
   #  bound, and the first such iteration (0 = never)
   fh <- as.integer(res$floor_hits); ff <- as.integer(res$floor_first)
   names(fh) <- names(ff) <- colnames(Y)
   #  number of ITERATIONS in which at least one series was floored; this is
   #  not recoverable from fh, whose per-series counts cannot be intersected
   fi <- as.integer(res$floor_iters)

   out <- list(Lambda = Lambda, residual_prec = ps, B = B, Omega = Omega,
               V = V_kept, loglik = loglik,
               dir_prec = dir_prec, nclust = nclust,
               floor_hits = fh, floor_first = ff, floor_iters = fi,
               dims = list(T = T_obs, p = p, k = k, lags = lags,
                           n_kept = length(keep)),
               settings = list(nrun = nrun, burn = burn, thin = thin,
                               rank = rank_used,
                               prior_B = if (identical(rank, "full")) "dl" else prior_B,
                               init = init, init_perturb = init_perturb))
   class(out) <- "vapor"
   out
}

#' Posterior mean precision matrix
#' @param Lambda p x k x n array of loading draws.
#' @param residual_prec n x p matrix of idiosyncratic precision draws.
#' @return A p x p matrix.
#' @export
posterior_mean_omega <- function(Lambda, residual_prec) {
   d <- dim(Lambda); p <- d[1]; n <- d[3]
   acc <- matrix(0, p, p)
   for (i in seq_len(n)) {
      om <- tcrossprod(Lambda[, , i])
      diag(om) <- diag(om) + residual_prec[i, ]
      acc <- acc + om
   }
   acc / n
}