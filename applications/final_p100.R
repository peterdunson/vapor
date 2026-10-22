# ============================================================================
#  FRED-QD application, p = 100.  Rscript applications/final_p100.R
#
#  Writes p100_final.rds (raw draws per chain plus every table), the two
#  diagnostic figures, and p100_tables.rds.  Refits nothing it can read, but
#  this script IS the fit: run it once and the figure scripts follow.
#
#  Every setting is in CFG below.  Nothing else in the file sets a parameter.
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 200)

CFG <- list(
   LABEL = "FRED-QD p = 100, rank 3",
   P     = 100L,
   CSV   = "/Users/peterdunson/Desktop/vapor/applications/data/2026-05-qd.csv",
   OUT   = FITS_DIR,
   OUT_RDS = "p100_final.rds",
   #  the figure scripts read this name, so the pipeline writes it too
   FIG_RDS = "p100_results.rds",
   GROUPS  = "/Users/peterdunson/Desktop/vapor/applications/fredqd_groups.csv",

   #  model
   K = 25L, RANK = 3L,
   #  MCMC, identical for VAPOR, PrecFact, VAR + PrecFact and Kastner
   SEEDS = 101:104, NRUN = 100000L, BURN = 20000L, THIN = 80L,
   MAX_REPL = c(VAPOR = 3L, PrecFact = 3L, `VAR+PrecFact` = 3L,
                `VAR+PrecFact (rank 3)` = 3L),
   #  Kastner needs its own, larger budget: its chains abort inside
   #  bayesianVARs' factor-SV sampler at this length far more often than
   #  VAPOR's do, so a budget sized for VAPOR is spent before Kastner has a
   #  single chain.  Replacements run in parallel with the rest, so the cost of
   #  a failed one is the worker, not the wall clock.
   WORKERS = NULL, USE_CACHE = FALSE,
   JIT = 5e-3,                    # the jitter the config D runs used
   VAR_JIT = 1e-3,                # jitter for the bayesianVARs-based methods
   #  first-stage VAR, deliberately NOT lengthened: changing it would move
   #  B_hat and with it the glasso and VAR + PrecFact residuals
   VAR_DRAWS = 2000L, VAR_BURNIN = 1000L,
   RHO_GRID = 10^seq(-2, 0, length = 30),

   #  held-out evaluation
   N_HID = 100L, HID_SEED_BASE = 90000L, HID_FRAC = 0.3,

   #  reporting
   THRESHOLDS = c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50),
   TOP_N_PURITY = 100L,
   RATE_SERIES = c("FEDFUNDS","TB3MS","GS1","GS5","GS10"),
   PPC_N_REP = 400L, PPC_SEED = 7001L, MSPE_REF = 0.5651570802, MSPE_TOL = 1e-6,   # config D mean, full precision
   TRACE_EDGES = 6L,
   FIG_TRACE = "diagnostics_p100_trace.png",
   FIG_PPC   = "diagnostics_p100_ppc.png",

   #  windows, in decimal years
   WIN   = c(1996 + 0/12, 2025 + 11/12),
   TRAIN = c(1996 + 0/12, 2020 + 11/12),
   TEST  = c(2021 + 2/12, 2025 + 11/12),
   VARS  = c("GDPC1","PCECC96","PCDGx","PCESVx","PCNDx","GPDIC1","FPIx","PNFIx","GCEC1","EXPGSC1",
      "IMPGSC1","DPIC96","INDPRO","IPFINAL","IPCONGD","IPMAT","IPDMAT","IPNMAT","IPBUSEQ","IPMANSICS",
      "IPFUELS","TCU","CUMFNS","PAYEMS","USPRIV","MANEMP","USGOOD","DMANEMP","USCONS","USFIRE","USINFO",
      "USPBS","USLAH","USMINE","USTPU","USGOVT","USTRADE","CE16OV","CIVPART","UNRATE","UNRATESTx",
      "UNRATELTx","UEMPMEAN","UEMPLT5","UEMP27OV","CLAIMSx","HWIURATIOx","AWHMAN","AWHNONAG","AHETPIx",
      "CES2000000008x","HOUST","HOUST5F","PERMIT","HOUSTS","HOUSTW","USSTHPI","SPCS10RSA","CMRMTSPLx",
      "RSAFSx","AMDMNOx","ANDENOx","BUSINVx","ISRATIOx","PCECTPI","PCEPILFE","CPIAUCSL","CPILFESL",
      "PPIACO","WPSFD49207","WPSID61","OILPRICEx","CPITRNSL","CPIMEDSL","FEDFUNDS","TB3MS","GS1","GS5",
      "GS10","AAA","BAA","MORTGAGE30US","GS10TB3Mx","BAA10YM","M1REAL","M2REAL","BOGMBASEREALx","TOTRESNS",
      "BUSLOANSx","CONSUMERx","REALLNx","NONREVSLx","S.P.500","NASDAQCOM","VIXCLSx","NIKKEI225","EXJPUSx",
      "EXUSUKx","EXCAUSx","UMCSENTx"))

#  MAIN rows use the full-rank bayesianVARs posterior-mean first stage; the
#  two "(rank 3)" rows are the rank-3-truncated sequential versions, reported
#  in the supplement.
ORDER <- c("PrecFact (no VAR)", "VAR + Glasso", "VAR + PrecFact",
           "Kastner (VAR+FSV)", "Sample covariance", "Random walk",
           "VAPOR (rank 3)", "VAR + Glasso (rank 3)", "VAR + PrecFact (rank 3)")

source("/Users/peterdunson/Desktop/vapor/applications/final_pipeline.R")
