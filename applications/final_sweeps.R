# ============================================================================
#  Sensitivity sweeps: k, rank, a_L and seeds, at either dimension.
#     Rscript applications/final_sweeps.R <29|100>
#
#  One 100,000-iteration chain per setting.  Settings at the REPORTED
#  configuration (k = CFG$K, rank = CFG$RANK, a_L = 0.5) and seed in 101-104 are
#  NOT refitted: those chains already exist in the main pipeline's output and
#  are read from it, so the sweep costs only the settings that differ.  Seeds
#  105 and 106 are fitted here.
#
#  Fits that do have to run go concurrently, up to physical cores minus one.
#  Every fit's Omega is saved, so no count or MSPE has to be refitted later.
#
#  Writes sweeps_p<dim>_final.rds.  Every setting is in CFG below.
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 200)
DIM <- as.integer(commandArgs(trailingOnly = TRUE)[1])
if (!DIM %in% c(29L, 100L)) stop("usage: Rscript final_sweeps.R <29|100>")

CFG <- list(
   LABEL = sprintf("sensitivity sweeps, p = %d", DIM),
   P = DIM,
   CSV = "/Users/peterdunson/Desktop/vapor/applications/data/2026-05-qd.csv",
   OUT = FITS_DIR,
   OUT_RDS  = sprintf("sweeps_p%d_final.rds", DIM),
   BASE_RDS = sprintf("p%d_final.rds", DIM),   # panel, held-out sets, and the reused chains

   #  reported configuration and the production MCMC length
   K    = if (DIM == 29L) 10L else 25L,
   RANK = 3L, A_L = 0.5,
   NRUN = 100000L, BURN = 20000L, THIN = 80L,
   JIT  = 5e-3,                                 # the config D jitter, both dimensions
   VAR_JIT = 1e-3, WORKERS = NULL,

   #  the grids
   K_GRID    = if (DIM == 29L) c(5L, 8L, 10L, 15L) else c(15L, 25L, 35L),
   RANK_GRID = c(2L, 3L, 4L),
   A_L_GRID  = c(0.3, 0.5, 1.0),
   SEED_GRID = 101:106,
   REUSE_SEEDS = 101:104,                       # read from BASE_RDS, never refitted

   N_HID = 100L, HID_SEED_BASE = 90000L, HID_FRAC = 0.3,
   #  corrected panel: both dimensions now run 1996Q1-2020Q4, matching
   #  final_p29.R and final_p100.R
   WIN   = c(1996 + 0/12, 2025 + 11/12),
   TRAIN = c(1996 + 0/12, 2020 + 11/12),
   TEST  = c(2021 + 2/12, 2025 + 11/12))

BASE <- readRDS(file.path(CFG$OUT, CFG$BASE_RDS))
CFG$VARS <- BASE$vars
stopifnot(identical(BASE$config$K, CFG$K), identical(BASE$config$RANK, CFG$RANK),
          identical(BASE$config$NRUN, CFG$NRUN), identical(BASE$config$THIN, CFG$THIN),
          isTRUE(all.equal(BASE$config$JIT, CFG$JIT)))

D <- build_panel(CFG)
Y_train <- D$train; Y_test <- D$test; vars <- D$vars; p <- D$p
stopifnot(isTRUE(all.equal(unname(Y_train), unname(BASE$Y_train), tolerance = 1e-12)))
hid_sets <- make_hid_sets(CFG$HID_SEED_BASE, CFG$N_HID, p, nrow(Y_test), CFG$HID_FRAC)

#  the grid, with the reported configuration marked
GRID <- rbind(
   data.frame(block="k",    k=CFG$K_GRID,  rank=CFG$RANK, a_L=CFG$A_L, seed=CFG$REUSE_SEEDS[1]),
   data.frame(block="rank", k=CFG$K, rank=CFG$RANK_GRID,  a_L=CFG$A_L, seed=CFG$REUSE_SEEDS[1]),
   data.frame(block="a_L",  k=CFG$K, rank=CFG$RANK, a_L=CFG$A_L_GRID,  seed=CFG$REUSE_SEEDS[1]),
   data.frame(block="seed", k=CFG$K, rank=CFG$RANK, a_L=CFG$A_L, seed=CFG$SEED_GRID))
GRID$setting <- with(GRID, ifelse(block=="k", sprintf("k=%d", k),
   ifelse(block=="rank", sprintf("rank=%d", rank),
   ifelse(block=="a_L", sprintf("a_L=%.1f", a_L), sprintf("seed=%d", seed)))))
GRID$reused <- with(GRID, k == CFG$K & rank == CFG$RANK & a_L == CFG$A_L &
                          seed %in% CFG$REUSE_SEEDS)
#  the reused chains, by seed, from the main pipeline
reuse_by_seed <- setNames(BASE$draws$vapor,
   vapply(BASE$draws$vapor, function(f) as.character(attr(f, "seed")), character(1)))
GRID$reused <- GRID$reused & as.character(GRID$seed) %in% names(reuse_by_seed)
cat(sprintf("%d settings: %d reused from %s, %d to fit\n", nrow(GRID),
            sum(GRID$reused), CFG$BASE_RDS, sum(!GRID$reused)))

one_fit <- function(k, rank, a_L, seed){
   f <- NULL; jj <- CFG$JIT; a <- 0L
   while (is.null(f) && a < 4L){ a <- a + 1L
      f <- tryCatch(vapor(Y_train, lags=1, k=k, rank=rank, prior_B="dl",
                          a_L=a_L, b_L=a_L, nrun=CFG$NRUN, burn=CFG$BURN,
                          thin=CFG$THIN, seed=seed, jitter=jj, init="random"),
                    error = function(e) NULL)
      jj <- jj*2 }
   if (!is.null(f)) attr(f, "attempts") <- a
   f }

#  fit the settings that are not reused, concurrently
todo <- which(!GRID$reused)
W <- n_workers(CFG$WORKERS)
cat(sprintf("fitting %d setting(s) on %d worker(s)\n", length(todo), W))
FITS <- parallel::mclapply(todo, function(i){ g <- GRID[i,]
   t0 <- Sys.time(); f <- one_fit(g$k, g$rank, g$a_L, g$seed)
   list(i = i, fit = f, minutes = as.numeric(difftime(Sys.time(), t0, units="mins")))
}, mc.cores = W, mc.preschedule = FALSE)

summarise_fit <- function(f, g, minutes, source){
   if (is.null(f)) return(list(setting=g$setting, block=g$block, status="failed",
                               k=g$k, rank=g$rank, a_L=g$a_L, seed=g$seed))
   Om <- omega_hat(f); M <- pair_draws(f, p); sel <- ci_indices(M)
   #  item 18: B is stored for every fit, and rho(B) summarised from it
   rho <- pr_explosive(f, p)
   list(setting=g$setting, block=g$block, status="completed", source=source,
        k=g$k, rank=g$rank, a_L=g$a_L, seed=g$seed,
        minutes=minutes, attempts=if (is.null(attr(f,"attempts"))) NA_integer_ else attr(f,"attempts"),
        Omega=Om, pcor=pcor_from_omega(Om),
        B=f$B, median_rho=unname(rho["med"]), pr_rho_gt_1=unname(rho["pr"]),
        edges=pair_labels(vars)[sel], n_edges=length(sel),
        mspe=mean(loss_avg(Om, hid_sets, Y_test, p)),
        floor_hits=f$floor_hits,
        floor_iters=if (is.null(f$floor_iters)) NA_integer_ else f$floor_iters) }

RES <- vector("list", nrow(GRID))
for (r in FITS){ if (inherits(r, "try-error") || is.null(r)) next
   g <- GRID[r$i,]; RES[[r$i]] <- summarise_fit(r$fit, g, r$minutes, "fitted here") }
for (i in which(GRID$reused)){ g <- GRID[i,]
   f <- reuse_by_seed[[as.character(g$seed)]]
   RES[[i]] <- summarise_fit(f, g, attr(f, "minutes"), CFG$BASE_RDS) }

ok <- Filter(function(z) !is.null(z) && identical(z$status, "completed"), RES)
tab <- do.call(rbind, lapply(ok, function(z) data.frame(
   block=z$block, setting=z$setting, seed=z$seed, source=z$source,
   MSPE=round(z$mspe,4), edges=z$n_edges,
   median_rho=round(z$median_rho,3), pr_rho_gt_1=round(z$pr_rho_gt_1,3),
   minutes=round(z$minutes,1),
   series_floored=sum(z$floor_hits > 0))))
cat(sprintf("\n=== p = %d sensitivity sweep ===\n", DIM)); print(tab, row.names=FALSE)

sens_blocks <- Filter(function(z) z$block %in% c("k","rank","a_L"), ok)
alle <- table(unlist(lapply(sens_blocks, function(z) z$edges)))
stab <- data.frame(edge=names(alle), n_settings=as.integer(alle))
stab <- stab[order(-stab$n_settings, stab$edge), ]
cat(sprintf("\n=== edge stability over the %d k/rank/a_L settings ===\n", length(sens_blocks)))
print(head(stab, 40), row.names=FALSE)

edge_sets <- setNames(lapply(ok, function(z) z$edges),
                      vapply(ok, function(z) z$setting, character(1)))
saveRDS(list(config=CFG, grid=GRID, results=RES, table=tab, stability=stab,
             edge_sets=edge_sets,
             run_at=format(Sys.time())), file.path(CFG$OUT, CFG$OUT_RDS))
cat(sprintf("\nsaved %s\n", CFG$OUT_RDS))
