# ============================================================================
#  FRED-QD application, p = 29.  Rscript applications/final_p29.R
#
#  Writes p29_final.rds (raw draws per chain plus every table), the two
#  diagnostic figures, and p29_tables.rds.
#
#  Every setting is in CFG below.  Nothing else in the file sets a parameter.
#  Group purity and the rate block are p = 100 reports and are skipped here
#  (TOP_N_PURITY and RATE_SERIES are NULL).
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 200)

CFG <- list(
   LABEL = "FRED-QD p = 29, rank 3",
   P     = 29L,
   CSV   = "~/Desktop/2026-05-qd.csv",
   OUT   = FITS_DIR,
   OUT_RDS = "p29_final.rds",
   FIG_RDS = "p29_results.rds",
   GROUPS  = "/Users/peterdunson/Desktop/vapor/applications/fredqd_groups.csv",

   K = 10L, RANK = 3L,
   SEEDS = 101:104, NRUN = 100000L, BURN = 20000L, THIN = 80L,
   MAX_REPL = 3L, WORKERS = NULL, USE_CACHE = FALSE,
   JIT = 5e-3,                    # the jitter the config D runs used
   VAR_JIT = 1e-3,
   VAR_DRAWS = 2000L, VAR_BURNIN = 1000L,
   RHO_GRID = 10^seq(-2, 0, length = 30),

   N_HID = 100L, HID_SEED_BASE = 90000L, HID_FRAC = 0.3,

   THRESHOLDS = c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50),
   TOP_N_PURITY = NULL,
   RATE_SERIES  = NULL,
   PPC_N_REP = 400L, PPC_SEED = 7029L, MSPE_REF = NA_real_, MSPE_TOL = 1e-6,   # no reference on the corrected panel
   TRACE_EDGES = 6L,
   FIG_TRACE = "diagnostics_p29_trace.png",
   FIG_PPC   = "diagnostics_p29_ppc.png",

   #  corrected panel: the NA filter now runs on the 29 selected series only,
   #  all of which are complete from 1996Q1, so the window matches p = 100
   WIN   = c(1996 + 0/12, 2025 + 11/12),
   TRAIN = c(1996 + 0/12, 2020 + 11/12),
   TEST  = c(2021 + 2/12, 2025 + 11/12),
   VARS  = c("GDPC1","INDPRO","IPFINAL","IPCONGD","IPBUSEQ","PAYEMS","UEMPMEAN",
      "CLAIMSx","CUMFNS","PCECC96","CMRMTSPLx","BUSINVx","CPIAUCSL","PPICMM","GDPCTPI",
      "OILPRICEx","FEDFUNDS","GS1","GS10","BAA","TB3MS","M2REAL","BUSLOANSx","NONREVSLx",
      "S.P.500","VIXCLSx","GPDIC1","UNRATE","WPSFD49207"))

#  MAIN rows use the full-rank bayesianVARs posterior-mean first stage; the
#  two "(rank 3)" rows are the rank-3-truncated sequential versions, reported
#  in the supplement.
ORDER <- c("PrecFact (no VAR)", "VAR + Glasso", "VAR + PrecFact",
           "Kastner (VAR+FSV)", "Sample covariance", "Random walk",
           "VAPOR (rank 3)", "VAR + Glasso (rank 3)", "VAR + PrecFact (rank 3)")

source("/Users/peterdunson/Desktop/vapor/applications/final_pipeline.R")
