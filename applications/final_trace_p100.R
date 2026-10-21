# ============================================================================
#  diagnostics_p100_trace.png, with GDPC1--GPDIC1 appended as one more panel
#  after the existing edge panels.  Same fig_trace() call as the pipeline, so
#  every style element is unchanged; only the edge list is longer.
#     Rscript applications/final_trace_p100.R
#  Reads p100_results.rds; refits nothing.
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
R <- readRDS(file.path(FITS_DIR, "p100_results.rds"))
CFG <- R$config; vars <- R$vars; p <- length(vars)
ed <- R$long$ci_edges
EXTRA <- data.frame(v1 = "GDPC1", v2 = "GPDIC1",
                    pcor = ed$pcor[ed$v1 == "GDPC1" & ed$v2 == "GPDIC1"][1])
if (is.na(EXTRA$pcor)){
   #  not among the selected edges: take its pooled mean from the draws
   DR <- lapply(R$draws$vapor, pair_draws, p = p)
   lab <- pair_labels(vars); j <- match("GDPC1--GPDIC1", lab)
   EXTRA$pcor <- mean(do.call(rbind, DR)[, j]) }
keep <- head(ed, CFG$TRACE_EDGES)
if (!any(keep$v1 == "GDPC1" & keep$v2 == "GPDIC1")) keep <- rbind(keep, EXTRA)
cat(sprintf("panels: log-likelihood, spectral radius, then %d edges:\n", nrow(keep)))
print(data.frame(panel = seq_len(nrow(keep)),
                 edge = paste(keep$v1, keep$v2, sep = "--"),
                 pcor = round(keep$pcor, 4)), row.names = FALSE)
OUT <- file.path(FIG_DIR, CFG$FIG_TRACE)
f <- fig_trace(R$draws$vapor, keep, vars, p, OUT, n_edge = nrow(keep))
cat(sprintf("wrote %s\n", OUT))
