# ============================================================================
#  Figure: Simulation A partial-correlation heatmaps
#  One replicate, p=20, T=100, rho(B)=0.95
# ============================================================================

library(VAPOR); library(PrecFact); library(bayesianVARs)
library(glasso); library(MASS); library(ggplot2); library(scales)

source("~/Desktop/precision_factor/PrecFact/vignettes/simulate_data_fxns.R")
OUT <- "~/Desktop/precision_factor/PrecFact/external_testing/var_prec"

REP     <- 1        # which replicate to display
NRUN    <- 5000
THIN    <- 5
REORDER <- TRUE     # cluster variables so the true blocks are contiguous

# ── DGP, identical to the Simulation A script ────────────────────────────────
simulate_var1 <- function(T, p, target_radius=0.89, k=5,
                          pr=0.2, val=1.5, seed=NULL) {
   if (!is.null(seed)) set.seed(seed)
   lambda <- simulate_lambda_sparse(k, p, pr=pr, val=val)
   Omega  <- tcrossprod(lambda) + diag(p)
   Sigma  <- chol2inv(chol(Omega))
   B      <- matrix(0, p, p)
   B[sample(p^2, floor(0.1*p^2))] <- rnorm(floor(0.1*p^2), 0, 0.3)
   B      <- B * (target_radius / max(abs(eigen(B, only.values=TRUE)$values)))
   Y      <- matrix(0, T, p)
   eps    <- MASS::mvrnorm(T, rep(0,p), Sigma)
   Y[1,]  <- eps[1,]
   for (t in 2:T) Y[t,] <- B %*% Y[t-1,] + eps[t,]
   list(Y=Y, B=B, Omega_true=Omega)
}

fit_glasso_oracle <- function(Y_input, Omega_true,
                              rho_grid=10^seq(-2, 0, length=20)) {
   S <- cov(scale(Y_input)); best <- list(frob=Inf)
   for (rho in rho_grid) {
      gl <- tryCatch(glasso(S, rho=rho, penalize.diagonal=FALSE),
                     error=function(e) NULL)
      if (is.null(gl)) next
      frob <- norm(cov2cor(gl$wi) - cov2cor(Omega_true), "F")
      if (frob < best$frob) best <- list(Omega=gl$wi, frob=frob)
   }
   best
}

fit_precfact <- function(Y_input, k=5) {
   res <- prec_est(as=0.1, bs=0.1, a=0.5, b=0.5, k=k, nrun=NRUN, thin=THIN,
                   Y=scale(Y_input), alldiff=0, a_dir=0.1, b_dir=0.1)
   fact_prec_processing(res, burnout=(NRUN/2)/THIN, alpha=0.1,
                        beta=0.85, epsilon_grid=seq(0.05,0.1,length=20))
}

get_var_residuals <- function(Y, seed) {
   set.seed(seed)
   Y_jit <- Y + matrix(rnorm(prod(dim(Y)), 0, 1e-6), nrow(Y), ncol(Y))
   mod   <- bvar(Y_jit, lags=1L, draws=2000, burnin=1000, quiet=TRUE)
   B_hat <- apply(mod$PHI, c(1,2), mean)
   mod$Y - mod$X %*% B_hat
}

fit_vapor <- function(Y, k=5, seed=99) {
   has_init <- "init" %in% names(formals(vapor))
   args <- list(Y, lags=1, k=k, rank="full", prior_B="dl",
                nrun=NRUN, burn=NRUN/2, thin=THIN, seed=seed)
   if (has_init) args$init <- "random"
   fit <- do.call(vapor, args)
   fact_prec_processing(list(Lambda=fit$Lambda, residuals=fit$residual_prec),
                        burnout=0, alpha=0.1, beta=0.85,
                        epsilon_grid=seq(0.05,0.1,length.out=20))
}

fit_kastner <- function(Y, k=5, seed=99) {
   set.seed(seed)
   Y_jit <- Y + matrix(rnorm(prod(dim(Y)), 0, 1e-6), nrow(Y), ncol(Y))
   p <- ncol(Y)
   sp <- specify_prior_sigma(M=p, type="factor", factor_factors=k,
                             factor_priorfacloadtype="colwiseng",
                             factor_heteroskedastic=TRUE, quiet=TRUE)
   mod <- bvar(Y_jit, lags=1L, draws=2000, burnin=1000,
               prior_sigma=sp, quiet=TRUE)
   solve(apply(vcov(mod)[1,,,], c(1,2), mean))
}

# ── Fit ──────────────────────────────────────────────────────────────────────
sim   <- simulate_var1(T=100, p=20, target_radius=0.95, k=5, seed=REP)
resid <- get_var_residuals(sim$Y, REP)

mats <- list(
   "True"           = sim$Omega_true,
   "VAPOR"          = fit_vapor(sim$Y, seed=REP)$Omega,
   "VAR + PrecFact" = fit_precfact(resid)$Omega,
   "VAR + Glasso"   = fit_glasso_oracle(resid, sim$Omega_true)$Omega,
   "Kastner"        = fit_kastner(sim$Y, seed=REP)
)

# ── Partial correlation, note the negation ───────────────────────────────────
pcor <- function(Omega) { R <- -cov2cor(Omega); diag(R) <- NA; R }
mats <- lapply(mats, pcor)

p_dim <- ncol(mats[[1]])
ord <- if (REORDER) {
   D <- 1 - abs(replace(mats[["True"]], is.na(mats[["True"]]), 0))
   hclust(as.dist(D), method="average")$order
} else seq_len(p_dim)

df <- do.call(rbind, lapply(names(mats), function(nm) {
   M <- mats[[nm]][ord, ord]
   data.frame(row    = rep(seq_len(p_dim), times=p_dim),
              col    = rep(seq_len(p_dim), each=p_dim),
              value  = as.vector(M),
              method = factor(nm, levels=names(mats)))
}))

lim <- 1

# ── Plot ─────────────────────────────────────────────────────────────────────
g <- ggplot(df, aes(col, row, fill=value)) +
   geom_tile(colour="grey88", linewidth=0.08) +
   facet_wrap(~ method, nrow=1) +
   scale_fill_gradient2(low="#2166AC", mid="#FFFFFF", high="#B2182B",
                        midpoint=0, limits=c(-lim, lim),
                        oob=scales::squish, na.value="grey80",
                        name=NULL,
                        breaks=c(-lim, 0, lim),
                        labels=sprintf("%.1f", c(-lim, 0, lim)),
                        guide=guide_colourbar(title=NULL, ticks=FALSE)) +
   scale_x_continuous(expand=c(0,0)) +
   scale_y_reverse(expand=c(0,0)) +
   coord_fixed() +
   theme_void(base_size=10) +
   theme(strip.text        = element_text(size=10, margin=margin(b=4)),
         panel.background  = element_rect(fill="white", colour=NA),
         plot.background   = element_rect(fill="white", colour=NA),
         panel.border      = element_rect(fill=NA, colour="grey45", linewidth=0.4),
         panel.spacing     = unit(6, "pt"),
         legend.position   = "bottom",
         legend.key.width  = unit(40, "pt"),
         legend.key.height = unit(5, "pt"),
         legend.text       = element_text(size=8),
         legend.margin     = margin(t=2),
         plot.margin       = margin(3, 3, 3, 3))


ggsave(file.path(path.expand(OUT), "sim_a_heatmap_rho95_T100.png"),
       g, width=6.5, height=1.85, units="in", dpi=600, bg="white")
