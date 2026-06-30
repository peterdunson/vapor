#' @export
print.vapor <- function(x, ...) {
   cat("VAPOR fit\n")
   cat(sprintf("  T = %d, p = %d | lags = %d, factors = %d\n",
               x$dims$T, x$dims$p, x$dims$lags, x$dims$k))
   cat(sprintf("  MCMC: %d iters, %d burn, thin %d (%d kept)\n",
               x$settings$nrun, x$settings$burn, x$settings$thin, x$dims$n_kept))
   invisible(x)
}

#' @export
summary.vapor <- function(object, ...) {
   R <- -stats::cov2cor(object$Omega); diag(R) <- NA
   pc <- R[upper.tri(R)]
   cat("VAPOR fit summary\n")
   cat(sprintf("  T = %d, p = %d | lags = %d, factors = %d\n",
               object$dims$T, object$dims$p, object$dims$lags, object$dims$k))
   cat(sprintf("  partial corr: mean |pc| = %.3f, max |pc| = %.3f\n",
               mean(abs(pc)), max(abs(pc))))
   invisible(object)
}

#' Partial correlation matrix from a vapor fit
#' @param object A "vapor" object.
#' @return A p x p partial correlation matrix.
#' @export
partial_correlations <- function(object) {
   R <- -stats::cov2cor(object$Omega); diag(R) <- NA
   R
}