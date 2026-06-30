#' Select a sparse precision graph via pFDR thresholding
#'
#' Applies the posterior false discovery rate (pFDR) edge-selection procedure to
#' a fitted VAPOR model. For each off-diagonal entry of the precision matrix this
#' computes the posterior probability that the absolute partial correlation
#' exceeds a threshold, then selects edges so the expected proportion of false
#' edges stays below `alpha`.
#'
#' @param object A fitted `"vapor"` object.
#' @param alpha Target pFDR level. Defaults to 0.1.
#' @param beta Posterior probability cutoff for declaring an edge. Defaults to 0.7.
#' @param epsilon_grid Grid of magnitude thresholds searched. Defaults to
#'   `seq(0.05, 0.3, length.out = 15)`.
#'
#' @return A list with `Omega` (posterior mean precision matrix), `P` (a 0/1
#'   adjacency matrix of selected edges, zero on the diagonal), and `pFDR` (the
#'   realized pFDR).
#' @export
select_graph <- function(object,
                         alpha = 0.1,
                         beta  = 0.7,
                         epsilon_grid = seq(0.05, 0.3, length.out = 15)) {
   if (!inherits(object, "vapor")) stop("`object` must be a 'vapor' fit.")
   
   Lambda <- object$Lambda          # p x k x n
   ps     <- object$residual_prec   # n x p
   n      <- dim(Lambda)[3]
   p      <- dim(Lambda)[1]
   
   upper <- function(x) x[upper.tri(x)]
   
   # |partial correlation| matrices for each draw
   abs_pc <- array(0, dim = c(p, p, n))
   for (i in seq_len(n)) {
      om <- tcrossprod(Lambda[, , i])
      diag(om) <- diag(om) + ps[i, ]
      abs_pc[, , i] <- abs(stats::cov2cor(om))
   }
   
   # search the threshold grid
   results <- lapply(epsilon_grid, function(ct) {
      postprob <- apply(abs_pc, c(1, 2), function(x) mean(x > ct))
      decs     <- postprob > beta
      pfdr     <- sum(upper(decs * (1 - postprob))) / max(sum(upper(decs)), 1)
      list(graph = decs, pFDR = pfdr)
   })
   
   FDRs <- vapply(results, function(r) r$pFDR, numeric(1))
   ok   <- which(FDRs < alpha)
   if (length(ok)) {
      results <- results[ok]; FDRs <- FDRs[ok]
      ind <- which.max(FDRs)
   } else {
      warning(sprintf("No pFDR below alpha = %g. Try increasing beta.", alpha))
      ind <- which.min(FDRs)
   }
   
   P <- results[[ind]]$graph * 1
   diag(P) <- 0
   
   list(
      Omega = object$Omega,
      P     = P,
      pFDR  = FDRs[ind]
   )
}