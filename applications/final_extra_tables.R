# ============================================================================
#  Items 2g (rho(B)), 2m (condition number), 4 (p = 100 glasso and the
#  cross-dimension table).  Reads results/fits only; fits nothing.
#     Rscript applications/final_extra_tables.R
#  Writes extra_tables.rds.
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 205)
GROUPS <- file.path(dirname(RES_DIR), "fredqd_groups.csv")
OUT <- list()

for (d in c(29L, 100L)){
   R <- readRDS(file.path(FITS_DIR, sprintf("p%d_final.rds", d)))
   ch <- R$draws$vapor; p <- R$p; vars <- R$vars; key <- as.character(d)

   # ---- 2. spectral radius of B from the four final chains ----------------
   rho <- unlist(lapply(ch, function(f)
      vapply(seq_len(dim(f$B)[3]), function(s) spec_rad(fixB(f$B[,,s], p)), numeric(1))))
   per <- do.call(rbind, lapply(ch, function(f){
      r <- vapply(seq_len(dim(f$B)[3]), function(s) spec_rad(fixB(f$B[,,s], p)), numeric(1))
      data.frame(chain = sprintf("seed %d", attr(f, "seed")), draws = length(r),
                 median_rho = round(median(r), 3), pr_gt_1 = round(mean(r > 1), 4)) }))
   pooled <- data.frame(chain = "pooled", draws = length(rho),
                        median_rho = round(median(rho), 3), pr_gt_1 = round(mean(rho > 1), 4))
   OUT[[sprintf("rho_p%d", d)]] <- rbind(per, pooled)
   cat(sprintf("\n=== 2. p = %d: spectral radius of B, four final chains ===\n", d))
   print(rbind(per, pooled), row.names = FALSE)

   # ---- 2b. the same for the sweep settings, where it can be computed -----
   S <- readRDS(file.path(FITS_DIR, sprintf("sweeps_p%d_final.rds", d)))
   okS <- Filter(function(z) !is.null(z) && identical(z$status,"completed"), S$results)
   reuse <- setNames(lapply(ch, function(f){
      r <- vapply(seq_len(dim(f$B)[3]), function(s) spec_rad(fixB(f$B[,,s], p)), numeric(1))
      c(med = median(r), pr = mean(r > 1)) }),
      vapply(ch, function(f) as.character(attr(f, "seed")), character(1)))
   sw <- do.call(rbind, lapply(okS, function(z){
      got <- identical(z$source, sprintf("p%d_final.rds", d)) &&
             as.character(z$seed) %in% names(reuse)
      data.frame(block = z$block, setting = z$setting, seed = z$seed, source = z$source,
                 MSPE = round(z$mspe, 4), edges = z$n_edges,
                 median_rho = if (got) round(reuse[[as.character(z$seed)]]["med"], 3) else NA_real_,
                 pr_rho_gt_1 = if (got) round(reuse[[as.character(z$seed)]]["pr"], 4) else NA_real_) }))
   OUT[[sprintf("sweep_rho_p%d", d)]] <- sw
   cat(sprintf("\n=== 2b. p = %d: sensitivity sweep with Pr(rho(B) > 1) ===\n", d))
   cat("    B was not saved for settings fitted inside the sweep, so the two\n")
   cat("    right-hand columns are available only for the rows whose chains come\n")
   cat("    from the main pipeline; the rest need a refit.\n")
   print(sw, row.names = FALSE)

   # ---- 3. held-out sets in which VAPOR is lower --------------------------
   hid <- R$hid_sets; Y_test <- R$Y_test
   setloss <- function(oms) rowMeans(vapply(oms, function(Om)
      vapply(hid, function(hl) mean(loss_one(Om, hl, Y_test, p)), numeric(1)),
      numeric(length(hid))))
   OM <- R$omegas
   if (d == 100L){
      kf <- file.path(FITS_DIR, "kastner_p100_20k.rds")
      if (file.exists(kf)){ KK <- readRDS(kf)
         if (length(KK$chains)) OM[["Kastner (VAR+FSV)"]] <- lapply(KK$chains, function(q) q$Omega) } }
   vap <- setloss(OM[["VAPOR (rank 3)"]])
   rw  <- vapply(hid, function(hl) mean(vapply(seq_len(nrow(Y_test)), function(t){
            hh <- hl[[t]]; yp <- if (t == 1) R$Y_train[nrow(R$Y_train), ] else Y_test[t-1, ]
            mean((Y_test[t,hh] - yp[hh])^2) }, numeric(1))), numeric(1))
   cmp <- list(); for (nm in names(OM)) if (nm != "VAPOR (rank 3)" && length(OM[[nm]]))
      cmp[[nm]] <- setloss(OM[[nm]])
   cmp[["Random walk"]] <- rw
   t3 <- do.call(rbind, lapply(names(cmp), function(nm) data.frame(
      Method = nm, sets = length(hid),
      VAPOR_lower = sum(vap < cmp[[nm]]),
      share = round(mean(vap < cmp[[nm]]), 3),
      median_gap = round(median(cmp[[nm]] - vap), 4))))
   OUT[[sprintf("sets_p%d", d)]] <- t3
   cat(sprintf("\n=== 3. p = %d: held-out sets (of %d) in which VAPOR's MSPE is lower ===\n",
               d, length(hid)))
   print(t3, row.names = FALSE)

   # ---- 5. edge stability denominator -------------------------------------
   blocks <- Filter(function(z) z$block %in% c("k","rank","a_L"), okS)
   t5 <- data.frame(n_settings_out_of = length(blocks),
                    settings = paste(vapply(blocks, function(z) z$setting, character(1)),
                                     collapse = ", "))
   OUT[[sprintf("stab_p%d", d)]] <- t5
   cat(sprintf("\n=== 5. p = %d: edge stability is out of %d settings ===\n", d, length(blocks)))
   cat(sprintf("    %s\n", t5$settings))

   # ---- 6. condition number of the training-sample covariance -------------
   Sg <- cov(R$Y_train); ev <- eigen(Sg, symmetric = TRUE, only.values = TRUE)$values
   t6 <- data.frame(p = p, T_train = nrow(R$Y_train), rank = qr(Sg)$rank,
                    min_eigenvalue = signif(min(ev), 4),
                    max_eigenvalue = signif(max(ev), 4),
                    condition = if (min(ev) <= .Machine$double.eps*max(ev)*p)
                       "singular (not of full rank)" else signif(max(ev)/min(ev), 6))
   OUT[[sprintf("cond_p%d", d)]] <- t6
   cat(sprintf("\n=== 6. p = %d: condition number of the training-sample covariance ===\n", d))
   print(t6, row.names = FALSE)
}

# ============================================================================
#  7. the p = 100 graphical lasso (rank-3 first stage, BIC fit)
# ============================================================================
R <- readRDS(file.path(FITS_DIR, "p100_final.rds"))
vars <- R$vars; p <- R$p
gl <- R$glasso; Om <- gl$Omega
g <- read.csv(GROUPS, stringsAsFactors = FALSE)
gc <- setNames(g$group_code, g$series); gn <- setNames(g$group_name, g$series)
ut <- which(upper.tri(matrix(0,p,p)), arr.ind = TRUE)
nz <- Om[cbind(ut[,1], ut[,2])] != 0
same <- gc[vars[ut[,1]]] == gc[vars[ut[,2]]]
t7a <- data.frame(nonzero_off_diagonal = sum(nz),
                  within_group = sum(nz & same),
                  share_within = round(mean(same[nz]), 3),
                  rho = round(gl$rho, 4))
cat("\n=== 7a. p = 100 glasso: nonzero off-diagonals and group membership ===\n")
print(t7a, row.names = FALSE); OUT$glasso_counts <- t7a

ed <- R$long$ci_edges
glpc <- pcor_from_omega(Om)
neg <- ed[ed$pcor < 0, ]
t7b <- data.frame(edge = paste(neg$v1, neg$v2, sep = "--"),
                  VAPOR_pcor = neg$pcor,
                  glasso_precision = signif(Om[cbind(match(neg$v1,vars), match(neg$v2,vars))], 4),
                  glasso_pcor = signif(glpc[cbind(match(neg$v1,vars), match(neg$v2,vars))], 4),
                  exactly_zero = Om[cbind(match(neg$v1,vars), match(neg$v2,vars))] == 0)
cat(sprintf("\n=== 7b. p = 100 glasso at VAPOR's %d negative edges ===\n", nrow(neg)))
print(t7b, row.names = FALSE); OUT$glasso_negative <- t7b

#  pooled VAPOR draws for every pair, for 7c and item 8
DR <- lapply(R$draws$vapor, pair_draws, p = p)
pool <- do.call(rbind, DR)
lab <- pair_labels(vars)
pm <- colMeans(pool); plo <- apply(pool, 2, quantile, .025); phi <- apply(pool, 2, quantile, .975)
sel <- lab[ci_indices(pool)]
glv <- abs(glpc[cbind(ut[,1], ut[,2])])
o <- order(-glv)
o <- o[!(lab[o] %in% sel)][1:25]
t7c <- data.frame(pair = lab[o],
                  glasso_pcor = round(glpc[cbind(ut[o,1], ut[o,2])], 3),
                  VAPOR_mean = round(pm[o], 3),
                  VAPOR_lo = round(plo[o], 3), VAPOR_hi = round(phi[o], 3),
                  same_group = unname(gc[vars[ut[o,1]]] == gc[vars[ut[o,2]]]))
cat("\n=== 7c. p = 100: glasso's 25 largest pairs that VAPOR does not select ===\n")
print(t7c, row.names = FALSE); OUT$glasso_top25 <- t7c

# ============================================================================
#  8. the seven p = 29 edges, located in the pooled p = 100 estimate
# ============================================================================
P29 <- readRDS(file.path(FITS_DIR, "p29_results.rds"))
e29 <- P29$long$ci_edges
rk <- rank(-abs(pm), ties.method = "min")
t8 <- do.call(rbind, lapply(seq_len(nrow(e29)), function(i){
   a <- e29$v1[i]; b <- e29$v2[i]
   if (!(a %in% vars) || !(b %in% vars))
      return(data.frame(edge = paste(a,b,sep="--"), in_p100 = FALSE, rank_of_4950 = NA,
                        p100_mean = NA, p100_lo = NA, p100_hi = NA, selected_at_p100 = NA))
   ia <- match(a, vars); ib <- match(b, vars)
   j <- which(ut[,1] == min(ia,ib) & ut[,2] == max(ia,ib))
   data.frame(edge = paste(a, b, sep = "--"), in_p100 = TRUE, rank_of_4950 = rk[j],
              p100_mean = round(pm[j], 3), p100_lo = round(plo[j], 3),
              p100_hi = round(phi[j], 3), selected_at_p100 = lab[j] %in% sel) }))
t8$p29_mean <- e29$pcor
cat("\n=== 8. the seven p = 29 edges in the pooled p = 100 estimate ===\n")
print(t8, row.names = FALSE); OUT$p29_in_p100 <- t8

saveRDS(OUT, file.path(FITS_DIR, "extra_tables.rds"))
cat("\nsaved extra_tables.rds\n")
