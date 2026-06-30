#' Fit a VAPOR model
#'
#' Joint Bayesian estimation of VAR dynamics and a sparse innovation precision
#' matrix.
#'
#' @param Y A numeric T x p matrix of stationary series (rows = time).
#' @param lags Integer VAR lag order.
#' @param k Integer number of precision factors.
#' @param nrun,burn,thin MCMC controls.
#' @param a_L,b_L Dirichlet-Laplace shape/rate for the loadings.
#' @param b_B Dirichlet-Laplace rate for the VAR coefficients.
#' @param as,bs Dirichlet process base-measure shape/rate.
#' @param a_dir,b_dir Dirichlet process concentration prior shape/rate.
#' @param jitter SD of Gaussian jitter for stability; NULL picks automatically.
#' @param seed Optional integer seed.
#' @return An object of class "vapor".
#' @export
vapor <- function(Y, lags = 1, k = 5, nrun = 5000, burn = floor(nrun / 2),
                  thin = 5, a_L = 0.5, b_L = 0.5, b_B = 0.5,
                  as = 0.1, bs = 0.1, a_dir = 0.1, b_dir = 0.1,
                  jitter = NULL, seed = NULL) {
   
   Y <- as.matrix(Y)
   if (!is.numeric(Y))      stop("`Y` must be numeric.")
   if (anyNA(Y))            stop("`Y` contains missing values.")
   if (burn >= nrun)        stop("`burn` must be less than `nrun`.")
   
   p <- ncol(Y); T_obs <- nrow(Y)
   if (!is.null(seed)) set.seed(seed)
   
   if (is.null(jitter)) jitter <- if (T_obs / p < 5) 1e-3 else 1e-6
   if (jitter > 0) Y <- Y + matrix(stats::rnorm(T_obs * p, 0, jitter), T_obs, p)
   
   res <- var_precision_sampler(
      Y = Y, lags = as.integer(lags), k = as.integer(k),
      nrun = as.integer(nrun), burn = as.integer(burn), thin = as.integer(thin),
      a_L = a_L, b_L = b_L, a_B = b_B, b_B = b_B,
      as = as, bs = bs, a_dir = a_dir, b_dir = b_dir
   )
   if (is.null(res)) stop("Sampler aborted.")
   
   n_mc <- dim(res$Lambda)[3]
   keep <- seq.int(floor(burn / thin) + 1L, n_mc)
   
   Lambda <- res$Lambda[, , keep, drop = FALSE]
   ps     <- res$residuals[keep, , drop = FALSE]
   B      <- res$B[, , keep, drop = FALSE]
   
   Omega <- posterior_mean_omega(Lambda, ps)
   
   out <- list(Lambda = Lambda, residual_prec = ps, B = B, Omega = Omega,
               dims = list(T = T_obs, p = p, k = k, lags = lags,
                           n_kept = length(keep)),
               settings = list(nrun = nrun, burn = burn, thin = thin))
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