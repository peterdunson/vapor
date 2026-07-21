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
#' @param a_B,b_B Dirichlet-Laplace shape/rate for the VAR coefficients
#'   (full B when rank = "full", or V when rank is numeric and
#'   prior_B = "dl").
#' @param prior_B Prior on the reduced-rank factor V. Ignored when
#'   rank = "full" (always Dirichlet-Laplace on the full B). `"dl"`
#'   (default) uses Dirichlet-Laplace shrinkage on V; `"ridge"` uses a
#'   constant-variance Gaussian prior.
#' @param s_U2,s_V2 Ridge prior variances for U and V. Only used when rank
#'   is numeric; s_V2 only applies when prior_B = "ridge".
#' @param as,bs Dirichlet process base-measure shape/rate.
#' @param a_dir,b_dir Dirichlet process concentration prior shape/rate.
#' @param jitter SD of Gaussian jitter for stability; NULL picks automatically.
#' @param seed Optional integer seed.
#' @return An object of class "vapor".
#' @export
vapor <- function(Y, lags = 1, k = 5, rank = "full",
                  nrun = 5000, burn = floor(nrun / 2), thin = 5,
                  a_L = 0.5, b_L = 0.5,
                  a_B = 0.5, b_B = 0.5,
                  prior_B = c("dl", "ridge"),
                  s_U2 = 1.0, s_V2 = 1.0,
                  as = 0.1, bs = 0.1, a_dir = 0.1, b_dir = 0.1,
                  dp_prior = c("escobar_west", "stirling_gamma"),
                  sg_a = 3.0, sg_b = 1.0, sg_m = 0,
                  jitter = NULL, seed = NULL) {
   
   Y <- as.matrix(Y)
   if (!is.numeric(Y))      stop("`Y` must be numeric.")
   if (anyNA(Y))            stop("`Y` contains missing values.")
   if (burn >= nrun)        stop("`burn` must be less than `nrun`.")
   
   prior_B <- match.arg(prior_B)
   dp_prior <- match.arg(dp_prior)
   
   p <- ncol(Y); T_obs <- nrow(Y)
   if (!is.null(seed)) set.seed(seed)
   
   if (is.null(jitter)) jitter <- if (T_obs / p < 5) 1e-3 else 1e-6
   if (jitter > 0) Y <- Y + matrix(stats::rnorm(T_obs * p, 0, jitter), T_obs, p)
   
   if (identical(rank, "full")) {
      
      res <- var_precision_sampler(
         Y = Y, lags = as.integer(lags), k = as.integer(k),
         nrun = as.integer(nrun), burn = as.integer(burn), thin = as.integer(thin),
         a_L = a_L, b_L = b_L, a_B = a_B, b_B = b_B,
         as = as, bs = bs, a_dir = a_dir, b_dir = b_dir,
         dp_prior = dp_prior, sg_a = sg_a, sg_b = sg_b, sg_m = sg_m
      )
      rank_used <- "full"
      
   } else {
      
      if (!is.numeric(rank) || rank < 1)
         stop("`rank` must be \"full\" or a positive integer.")
      max_rank <- min(p, lags * p)
      if (rank > max_rank)
         stop(sprintf("`rank` must be <= min(p, lags*p) = %d.", max_rank))
      
      res <- var_precision_sampler_rr(
         Y = Y, lags = as.integer(lags), k = as.integer(k), rank = as.integer(rank),
         nrun = as.integer(nrun), burn = as.integer(burn), thin = as.integer(thin),
         a_L = a_L, b_L = b_L, prior_B = prior_B,
         s_U2 = s_U2, s_V2 = s_V2, a_V = a_B, b_V = b_B,
         as = as, bs = bs, a_dir = a_dir, b_dir = b_dir,
         dp_prior = dp_prior, sg_a = sg_a, sg_b = sg_b, sg_m = sg_m
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
   
   out <- list(Lambda = Lambda, residual_prec = ps, B = B, Omega = Omega,
               V = V_kept, loglik = loglik,
               dir_prec = dir_prec, nclust = nclust,
               dims = list(T = T_obs, p = p, k = k, lags = lags,
                           n_kept = length(keep)),
               settings = list(nrun = nrun, burn = burn, thin = thin,
                               rank = rank_used,
                               prior_B = if (identical(rank, "full")) "dl" else prior_B))
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