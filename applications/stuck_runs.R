# ============================================================================
#  Item A: runs of 5 or more consecutive retained draws in which a reported
#  edge's partial correlation is EXACTLY identical, for every chain at both
#  dimensions.  Reads saved chains only; fits nothing, plots nothing.
#
#  "Exactly identical" is bit equality of the double, not equality to a
#  tolerance.  For each run the residual precision of both endpoint series is
#  checked against the 1e-8 floor over the draws of that run.
#
#  Writes stuck_runs.rds and appends to results_summary.md.
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 205)
APP <- FITS_DIR
MD  <- file.path(TAB_DIR, "results_summary.md")
MINRUN <- 5L; FLOOR <- 1e-8

find_runs <- function(x, minrun){
   n <- length(x)
   same <- x[-1] == x[-n]                     # bit equality between neighbours
   r <- rle(same)
   out <- NULL; pos <- 1L
   for (i in seq_along(r$lengths)){
      if (isTRUE(r$values[i])){
         start <- pos                          # first draw of the run
         len <- r$lengths[i] + 1L              # draws sharing one value
         if (len >= minrun) out <- rbind(out, data.frame(start = start, length = len)) }
      pos <- pos + r$lengths[i] }
   out }

ALL <- list()
for (d in c(29L, 100L)){
   f <- file.path(APP, sprintf("p%d_final.rds", d)); if (!file.exists(f)) next
   R <- readRDS(f); ch <- R$draws$vapor; vars <- R$vars; p <- R$p
   ed <- R$long$ci_edges
   cat(sprintf("\np = %d: %d chains, %d retained draws each, %d reported edges\n",
               d, length(ch), dim(ch[[1]]$Lambda)[3], nrow(ed)))
   rows <- NULL
   for (ci in seq_along(ch)){
      f_ <- ch[[ci]]; nd <- dim(f_$Lambda)[3]
      #  the partial correlation of every reported edge, per draw
      ij <- cbind(match(ed$v1, vars), match(ed$v2, vars))
      M <- matrix(NA_real_, nd, nrow(ed))
      for (s in seq_len(nd)){
         pc <- pcor_from_omega(omega_one(f_, s))
         M[s,] <- pc[cbind(ij[,1], ij[,2])] }
      for (e in seq_len(nrow(ed))){
         rr <- find_runs(M[,e], MINRUN); if (is.null(rr)) next
         for (k in seq_len(nrow(rr))){
            ix <- rr$start[k]:(rr$start[k] + rr$length[k] - 1L)
            pr <- f_$residual_prec[ix, ij[e,], drop = FALSE]
            rows <- rbind(rows, data.frame(
               chain = sprintf("seed %d", attr(f_, "seed")),
               edge = paste(ed$v1[e], ed$v2[e], sep="--"),
               start_draw = rr$start[k], run_length = rr$length[k],
               value = signif(M[rr$start[k], e], 6),
               v1_at_floor = sum(pr[,1] <= FLOOR), v2_at_floor = sum(pr[,2] <= FLOOR),
               either_at_floor = any(pr <= FLOOR))) } } }
   if (is.null(rows)){
      cat("no run of 5 or more identical consecutive draws in any chain or edge\n")
      ALL[[as.character(d)]] <- data.frame(result = "none found")
   } else {
      rows <- rows[order(-rows$run_length, rows$chain, rows$edge), ]
      print(rows, row.names = FALSE)
      cat(sprintf("\n%d run(s); %d involve at least one endpoint at the floor\n",
                  nrow(rows), sum(rows$either_at_floor)))
      ALL[[as.character(d)]] <- rows } }

cat("\n\n", file = MD, append = TRUE)
for (d in names(ALL))
   md_section(MD, sprintf("p = %s runs of >= %d identical consecutive draws, reported edges",
                          d, MINRUN), ALL[[d]])
saveRDS(ALL, file.path(APP, "stuck_runs.rds"))
cat(sprintf("\nsaved stuck_runs.rds and appended to %s\n", MD))
