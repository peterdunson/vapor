# ============================================================================
#  Append the simulation sections to results_summary.md, built only from the
#  current Gamma(1, 1) base-measure outputs.  The application sections are not
#  touched: this appends, and removes any simulation block it wrote before, so
#  it is idempotent.
#     Rscript applications/final_sim_sections.R
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 210)
MD   <- file.path(TAB_DIR, "results_summary.md")
SIMD <- path.expand("~/Desktop/precision_factor/PrecFact/external_testing/var_prec")
MARK <- "# Simulations"
ml <- readLines(MD, warn = FALSE)
cut <- grep(sprintf("^%s$", MARK), ml)[1]
if (!is.na(cut)) ml <- ml[seq_len(cut - 1L)]
writeLines(ml, MD)

sec <- function(title, x, note = NULL, src = NULL){
   if (is.null(x) || (is.data.frame(x) && !nrow(x))) return(invisible(NULL))
   cat(sprintf("\n## %s\n\n", title), file = MD, append = TRUE)
   if (!is.null(note)) cat(sprintf("%s\n\n", note), file = MD, append = TRUE)
   if (!is.null(src))  cat(sprintf("*Source: %s*\n\n", src), file = MD, append = TRUE)
   cat("```\n", file = MD, append = TRUE)
   cat(paste(utils::capture.output(print(x, row.names = FALSE)), collapse = "\n"),
       file = MD, append = TRUE)
   cat("\n```\n", file = MD, append = TRUE) }

A  <- readRDS(file.path(FITS_DIR, "sim_A_results.rds"))
B  <- readRDS(file.path(FITS_DIR, "sim_B_results.rds"))
Cc <- readRDS(file.path(FITS_DIR, "sim_C_results.rds"))
Dd <- readRDS(file.path(FITS_DIR, "sim_D_results.rds"))
Ff <- readRDS(file.path(FITS_DIR, "sim_f_final_results.rds"))
Gg <- readRDS("/Users/peterdunson/Desktop/vapor/simulation/sim_g_parts/sim_g_results.rds")
DIFF <- readRDS(file.path(FITS_DIR, "sim_tables_diff.rds"))$diffs

cat(sprintf("\n\n%s\n\nRerun under the DP base measure Gamma(as = 1, bs = 1), for VAPOR and for PrecFact; a_dir = b_dir = 0.1; residual-precision floor 1e-8. Same seeds, replicates and settings as before.\nPaper naming: Simulation E is sim_f.R and Simulation F is sim_g.R. **Simulation F has three settings, not four: p = 200, T = 140 ran 17 h to replicate 19 of 25 and was cancelled.**\nEvery number below is from this run; nothing is carried over from the Gamma(0.1, 0.1) tables.\n",
            MARK), file = MD, append = TRUE)

#  ---- 1. every cell of every table -----------------------------------------
sec("Simulation A, Experiment A (varying rho(B), T = 200): every cell",
    A$experiment_A$summary, "Mean (sd) over 25 replicates.", "sim_A_results.rds")
sec("Simulation A, Experiment B (varying T, rho(B) = 0.95): every cell",
    A$experiment_B$summary, NULL, "sim_A_results.rds")
sec("Simulation B, Experiment A (varying rho(B), T = 200): every cell",
    B$experiment_A$summary, NULL, "sim_B_results.rds")
sec("Simulation B, Experiment B (varying T, rho(B) = 0.95): every cell",
    B$experiment_B$summary, NULL, "sim_B_results.rds")
sec("Simulation C (factor stochastic volatility DGP): every cell",
    Cc$summary, "Held-out prediction MSPE, mean (sd) over 25 replicates.",
    "sim_C_results.rds")
sec("Simulation D (varying signal strength): every cell",
    Dd$summary, NULL, "sim_D_results.rds")
rr_cells <- function(R) data.frame(
   setting = sprintf("(%d,%d)", R$p, R$T), Method = R$Method,
   Omega_Frob = round(R$Omega_Frob, 4), Omega_Frob_sd = round(R$Omega_Frob_sd, 4),
   B_Frob = round(R$B_Frob, 4), B_Frob_sd = round(R$B_Frob_sd, 4),
   Sens = round(R$Sens, 4), Sens_sd = round(R$Sens_sd, 4),
   Spec = round(R$Spec, 4), Spec_sd = round(R$Spec_sd, 4))
sec("Simulation E (sim_f.R, approximately low-rank B, T >= p): every cell",
    rr_cells(Ff$results),
    "Omega_Frob is ||R_bar - R_0||_F; B_Frob is ||B_hat - B_0||_F/||B_0||_F (VAPOR only); Sens and Spec are at an oracle threshold retaining the true number of edges.",
    "sim_f_final_results.rds")
sec("Simulation F (sim_g.R, T < p at T/p = 0.7): every cell",
    rr_cells(Gg$results), "Three settings; p = 200, T = 140 was cancelled.",
    "sim_g_results.rds")

#  ---- 2. paired comparisons against VAPOR ----------------------------------
pn <- "Diff is the mean paired per-replicate difference against VAPOR, positive favouring VAPOR; SE is its standard error; p is from a paired t-test on 24 degrees of freedom; Wins counts replicates in which VAPOR is closer to the truth."
sec("Simulation A, Experiment A: paired differences in Frobenius distance",
    A$experiment_A$paired, pn, "sim_A_results.rds")
sec("Simulation A, Experiment B: paired differences in Frobenius distance",
    A$experiment_B$paired, pn, "sim_A_results.rds")
sec("Simulation B, Experiment A: paired differences in Frobenius distance",
    B$experiment_A$paired, pn, "sim_B_results.rds")
sec("Simulation B, Experiment B: paired differences in Frobenius distance",
    B$experiment_B$paired, pn, "sim_B_results.rds")
sec("Simulation C: paired differences in held-out MSPE", Cc$paired,
    sub("closer to the truth", "predicting more accurately",
        sub("Frobenius distance", "MSPE", pn)), "sim_C_results.rds")
sec("Simulation D: paired differences in Frobenius distance", Dd$paired, pn,
    "sim_D_results.rds")
sec("Simulation E: paired differences against VAPOR (DL)", Ff$paired,
    "Reference is VAPOR dl r=3. The 'VAPOR ridge r=3' row IS the DL-vs-ridge difference: positive favours DL.",
    "sim_f_final_results.rds")
sec("Simulation F: paired differences against VAPOR (DL)", Gg$paired,
    "Reference is VAPOR dl r=3; the 'VAPOR ridge r=3' row is the DL-vs-ridge difference.",
    "sim_g_results.rds")

#  ---- 3. E and F: coefficient error and oracle-threshold edge recovery ------
coef_tab <- function(R){
   v <- R[grepl("^VAPOR", R$Method), ]
   data.frame(setting = sprintf("(%d,%d)", v$p, v$T), Method = v$Method,
              rel_coef_error = round(v$B_Frob, 4), sd = round(v$B_Frob_sd, 4)) }
edge_tab <- function(R) data.frame(
   setting = sprintf("(%d,%d)", R$p, R$T), Method = R$Method,
   Sens = round(R$Sens, 4), Sens_sd = round(R$Sens_sd, 4),
   Spec = round(R$Spec, 4), Spec_sd = round(R$Spec_sd, 4))
sec("Simulation E: relative coefficient error of VAPOR", coef_tab(Ff$results),
    "||B_hat - B_0||_F / ||B_0||_F, mean (sd) over 25 replicates.", "sim_f_final_results.rds")
sec("Simulation E: edge sensitivity and specificity at the oracle threshold",
    edge_tab(Ff$results),
    "Glasso's values are reported but are degenerate: its off-diagonal mass is near uniform, so the oracle threshold assigns it almost all or almost no edges.",
    "sim_f_final_results.rds")
sec("Simulation F: relative coefficient error of VAPOR", coef_tab(Gg$results),
    NULL, "sim_g_results.rds")
sec("Simulation F: edge sensitivity and specificity at the oracle threshold",
    edge_tab(Gg$results), NULL, "sim_g_results.rds")

#  ---- 4. what changed against the archived Gamma(0.1, 0.1) tables ----------
if (!is.null(DIFF)){
   ct <- as.data.frame(table(DIFF$table)); names(ct) <- c("table","cells_differing")
   sec("Cells differing from the archived Gamma(0.1, 0.1) tables: count per table", ct,
       sprintf("%d cells differ in all. The two PrecFact rows account for %d of them: their base measure moved from 0.1 to 1 along with VAPOR's. Kastner and the Glasso rows move only through the paired differences and win counts, which are computed against VAPOR.",
               nrow(DIFF), sum(grepl("PrecFact", DIFF$method))), "sim_tables_diff.rds")
   bym <- as.data.frame(table(trimws(DIFF$method)))
   names(bym) <- c("method","cells_differing")
   bym <- bym[order(-bym$cells_differing), ]
   sec("Cells differing: count per method", bym, NULL, "sim_tables_diff.rds")
   sec("Every differing cell, archived vs new", DIFF,
       "Formatting is ignored: <1e-04 and $<10^{-4}$ compare equal.", "sim_tables_diff.rds") }

#  ---- provenance ------------------------------------------------------------
sec("Simulation run provenance",
    data.frame(
      simulation = c("A (sim_a)","B (sim_b)","C (sim_c)","D (sim_d)",
                     "E (sim_f)","F (sim_g) p=50","F (sim_g) p=100","F (sim_g) p=150",
                     "F (sim_g) p=200"),
      status = c(rep("complete", 8), "CANCELLED at replicate 19 of 25"),
      finished = c("2026-10-20 17:59","2026-10-20 20:48","2026-10-20 21:01",
                   "2026-10-21 14:31","2026-10-21 20:47","2026-10-22 09:31",
                   "2026-10-21 16:48","2026-10-21 20:24","-"),
      minutes = c(NA, NA, NA, NA, NA, 26.7, 74.3, 290.5, NA)),
    "sim_g's three settings were run as independent processes; sim_g.R re-seeds from the replicate index alone, so per-setting runs reproduce the sequential result exactly (verified: rerunning p = 50 gave results identical at tolerance 0). No floor-hit counts were recorded for this simulation run: the driver that logged them was stopped before writing simulations_floor_report.rds, and the per-setting runner does not wrap vapor().",
    "run logs")
n <- length(readLines(MD, warn = FALSE))
cat(sprintf("results_summary.md: %d lines, %d sections (%d simulation)\n", n,
            sum(grepl("^## ", readLines(MD, warn = FALSE))),
            sum(grepl("^## Simulation|^## Cells differing", readLines(MD, warn = FALSE)))))
