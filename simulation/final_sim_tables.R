# ============================================================================
#  Regenerate every simulation table from this run's outputs, in the exact
#  format of the archived sim_tables.tex, and report every numerically
#  differing cell.
#     Rscript simulation/final_sim_tables.R
#
#  The archived sim_tables.tex is the TEMPLATE: every line that is not a body
#  row is copied verbatim (caption, label, tabular preamble, headers, rules,
#  multirow) and only the numeric cells are replaced, each keeping the decimal
#  places, the $...$ wrapper and the "mean (sd)" shape of the cell it
#  replaces.  Bold is recomputed from the new numbers.
#
#  tab:sim-f now has THREE blocks: p = 200, T = 140 was cancelled.
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 200)
TPL  <- "/Users/peterdunson/Desktop/vapor/simulation/sim_tables_template.tex"
OUT  <- file.path(SIM_TAB, "sim_tables.tex")
tl <- readLines(TPL, warn = FALSE)

is_body <- function(ln) grepl("&", ln) && grepl("\\\\\\\\\\s*$", ln) &&
   !grepl("Scenario *&|\\(p,T\\)\\$ *&|multicolumn", ln)
nums <- function(cell) regmatches(cell, gregexpr("-?[0-9]+\\.?[0-9]*", cell))[[1]]
dps  <- function(cell){ m <- nums(cell)
   vapply(m, function(z){ i <- regexpr("\\.", z); if (i < 0) 0L else nchar(z) - i },
          integer(1), USE.NAMES = FALSE) }
dollar <- function(cell) grepl("\\$-?[0-9.]+\\$", cell)
dashes <- function(cell) grepl("^\\s*-{2,}\\s*$", gsub("\\\\\\\\", "", cell))
fmt <- function(old, v){
   if (is.null(v) || dashes(old)) return(old)
   pre  <- sub("^(\\s*).*", "\\1", old)
   tail <- regmatches(old, regexpr("\\s*(\\\\\\\\)?\\s*$", old)); if (!length(tail)) tail <- ""
   d <- dps(gsub("10\\^\\{?-?4\\}?", "", old)); d <- c(d, 3L, 3L)
   body <- switch(v$type,
      ms  = { a <- sprintf(paste0("%.", d[1], "f"), v$m)
              if (dollar(old)) a <- sprintf("$%s$", a)
              if (is.na(v$s)) a else sprintf("%s (%s)", a, sprintf(paste0("%.", d[2], "f"), v$s)) },
      num = { a <- sprintf(paste0("%.", d[1], "f"), v$v)
              if (dollar(old)) sprintf("$%s$", a) else a },
      str = v$v)
   if (isTRUE(v$bold)) body <- sprintf("\\textbf{%s}", body)
   paste0(pre, body, tail) }
pfmt <- function(s) if (grepl("<", s)) "$<10^{-4}$" else s
skey <- function(x){ x <- gsub("\\\\multirow\\{[0-9]+\\}\\{\\*\\}", "", x)
   v <- regmatches(x, gregexpr("[0-9]+\\.?[0-9]*", x))[[1]]
   if (!length(v)) return(NA_character_)
   paste(sprintf("%g", as.numeric(v)), collapse = ",") }
mkey <- function(x){ x <- gsub("\\\\text(tt|bf)\\{([^}]*)\\}", "\\2", x)
   tolower(gsub("[^A-Za-z0-9+]", "", x)) }
MMAP <- c(vapordl = "VAPOR dl r=3", vaporridge = "VAPOR ridge r=3")
ms_from <- function(s){ z <- as.numeric(regmatches(s, gregexpr("-?[0-9]+\\.?[0-9]*", s))[[1]])
   list(type = "ms", m = z[1], s = if (length(z) > 1) z[2] else NA_real_) }
skey_v <- function(x) vapply(x, skey, character(1), USE.NAMES = FALSE)
mkey_v <- function(x) vapply(x, mkey, character(1), USE.NAMES = FALSE)

A <- readRDS(file.path(SIM_FITS, "sim_A_results.rds"))
B <- readRDS(file.path(SIM_FITS, "sim_B_results.rds"))
Cc<- readRDS(file.path(SIM_FITS, "sim_C_results.rds"))
Dd<- readRDS(file.path(SIM_FITS, "sim_D_results.rds"))
Ee<- readRDS(file.path(SIM_FITS, "sim_e_final_results.rds"))
Ff<- readRDS("/Users/peterdunson/Desktop/vapor/simulation/results/sim_f_parts/sim_f_results.rds")

new_summary <- function(sm, pr, vcols){
   K <- list()
   for (i in seq_len(nrow(sm))){
      k <- paste(skey(sm$Scenario[i]), mkey(sm$Method[i]))
      cl <- lapply(vcols, function(cc) ms_from(sm[[cc]][i]))
      j <- which(skey_v(pr$Scenario) == skey(sm$Scenario[i]) &
                 mkey_v(pr$Method)   == mkey(sm$Method[i]))
      cl[[length(cl)+1L]] <- if (length(j) == 1)
         list(type = "ms", m = pr$Diff[j], s = pr$SE[j]) else NULL
      K[[k]] <- cl }
   K }
new_paired <- function(pr){
   K <- list()
   for (i in seq_len(nrow(pr))){
      K[[paste(skey(pr$Scenario[i]), mkey(pr$Method[i]))]] <-
         list(list(type="num", v=pr$Diff[i]), list(type="num", v=pr$SE[i]),
              list(type="num", v=pr$t[i]), list(type="str", v=pfmt(pr$p[i])),
              list(type="str", v=pr$Wins[i])) }
   K }
`%||%` <- function(a,b) if (is.na(a) || is.null(a)) b else a
new_rr <- function(R){
   K <- list()
   for (i in seq_len(nrow(R))){
      k <- paste(sprintf("%g,%g", R$p[i], R$T[i]),
                 mkey(names(MMAP)[match(R$Method[i], MMAP)] %||% R$Method[i]))
      K[[k]] <- list(list(type="ms", m=R$Omega_Frob[i], s=R$Omega_Frob_sd[i]),
                     if (is.finite(R$B_Frob[i])) list(type="ms", m=R$B_Frob[i], s=R$B_Frob_sd[i]) else NULL,
                     list(type="ms", m=R$Sens[i], s=R$Sens_sd[i]),
                     list(type="ms", m=R$Spec[i], s=R$Spec_sd[i])) }
   K }

DIFFS <- list()
tex_table <- function(label){
   i0 <- grep(sprintf("\\\\label\\{%s\\}", label), tl)[1]
   if (is.na(i0)) return(NULL)
   b <- max(grep("\\\\begin\\{table\\}", tl)[grep("\\\\begin\\{table\\}", tl) < i0])
   e <- min(grep("\\\\end\\{table\\}",   tl)[grep("\\\\end\\{table\\}",   tl) > i0])
   tl[b:e] }
rewrite <- function(label, KEY, bold_col, drop_scen = character(0)){
   tb <- tex_table(label); if (is.null(tb) || is.null(KEY)) return(tb)
   HD <- { h <- tb[grepl("Scenario *&|\\(p,T\\)\\$ *&", tb)][1]
           if (is.na(h)) NULL else trimws(gsub("[\\\\$]|\\\\\\\\\\s*$", "",
                                               strsplit(h, "&", fixed = TRUE)[[1]])) }
   sc <- NA_character_; rows <- vector("list", length(tb)); blk <- integer(length(tb)); bi <- 1L
   for (i in seq_along(tb)){ ln <- tb[i]
      if (grepl("\\\\midrule", ln)) bi <- bi + 1L
      if (grepl("multirow", ln) && !is_body(ln)) sc <- skey(ln)
      if (!is_body(ln)) next
      cs <- strsplit(ln, "&", fixed = TRUE)[[1]]
      s1 <- skey(cs[1]); if (!is.na(s1)) sc <- s1
      rows[[i]] <- list(cs = cs, sc = sc, mk = mkey(cs[2])); blk[i] <- bi }
   #  blocks whose scenario is being dropped
   kill <- vapply(seq_along(tb), function(i){
      r <- rows[[i]]; !is.null(r) && r$sc %in% drop_scen }, logical(1))
   killblk <- unique(blk[kill]); killblk <- killblk[killblk > 0]
   bold <- rep(FALSE, length(tb))
   if (!is.na(bold_col)) for (bb in setdiff(unique(blk[blk > 0]), killblk)){
      ii <- which(blk == bb)
      vv <- vapply(ii, function(i){ v <- KEY[[paste(rows[[i]]$sc, rows[[i]]$mk)]][[bold_col - 2L]]
         if (is.null(v)) NA_real_ else v$m }, numeric(1))
      dd <- dps(rows[[ii[1]]]$cs[bold_col])[1]
      vr <- round(vv, dd)
      if (any(!is.na(vr))) bold[ii[which(vr == min(vr, na.rm = TRUE))]] <- TRUE }
   out <- tb; drop_line <- rep(FALSE, length(tb))
   for (i in seq_along(tb)){
      if (blk[i] %in% killblk){ drop_line[i] <- TRUE; next }
      r <- rows[[i]]; if (is.null(r)) next
      cs <- r$cs; new <- KEY[[paste(r$sc, r$mk)]]
      if (is.null(new)) next
      for (cc in 3:length(cs)){
         v <- if (cc - 2L <= length(new)) new[[cc - 2L]] else NULL
         if (!is.null(v) && identical(cc, bold_col)) v$bold <- bold[i]
         old <- cs[cc]; cs[cc] <- fmt(old, v)
         if (!is.null(v) && !dashes(old)){
            cl <- function(z) gsub("\\\\textbf\\{|\\}|\\$|\\s", "", gsub("\\\\\\\\", "", z))
            a <- cl(old); bq <- cl(cs[cc])
            if (a != bq) DIFFS[[length(DIFFS)+1L]] <<- data.frame(
               table = label, scenario = r$sc,
               method = trimws(gsub("\\\\text(tt|bf)\\{([^}]*)\\}", "\\2", cs[2])),
               column = if (!is.null(HD) && cc <= length(HD)) HD[cc] else as.character(cc),
               archived = trimws(a), new = trimws(bq)) } }
      out[i] <- paste(cs, collapse = "&") }
   #  a killed block also has a \multirow header line and a \midrule before
   #  it; neither is a body row, so they must be dropped by name
   for (i in seq_along(tb)){
      if (grepl("multirow", tb[i]) && !is.na(skey(tb[i])) && skey(tb[i]) %in% drop_scen){
         drop_line[i] <- TRUE
         if (i > 1 && grepl("\\\\midrule", tb[i-1])) drop_line[i-1] <- TRUE } }
   out[!drop_line] }

TBL <- list(
   list("tab:sim-a-autocorr", new_summary(A$experiment_A$summary, A$experiment_A$paired, c("Sens","Spec","Frob")), 5L, character(0)),
   list("tab:sim-a-samplesize", new_summary(A$experiment_B$summary, A$experiment_B$paired, c("Sens","Spec","Frob")), 5L, character(0)),
   list("tab:sim-a-paired-autocorr", new_paired(A$experiment_A$paired), NA_integer_, character(0)),
   list("tab:sim-a-paired-samplesize", new_paired(A$experiment_B$paired), NA_integer_, character(0)),
   list("tab:sim-b-autocorr", new_summary(B$experiment_A$summary, B$experiment_A$paired, c("Sens","Spec","Frob")), 5L, character(0)),
   list("tab:sim-b-samplesize", new_summary(B$experiment_B$summary, B$experiment_B$paired, c("Sens","Spec","Frob")), 5L, character(0)),
   list("tab:sim-c", new_summary(Cc$summary, Cc$paired, c("MSPE")), 3L, character(0)),
   list("tab:sim-d", new_summary(Dd$summary, Dd$paired, c("Sens","Spec","Frob")), 5L, character(0)),
   list("tab:sim-e", new_rr(Ee$results), 3L, character(0)),
   list("tab:sim-f", new_rr(Ff$results), 3L, "200,140"))

lines <- c("%  Simulation tables, regenerated under the Gamma(1, 1) DP base measure.",
           sprintf("%%  %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
           "%  tab:sim-f: p = 200, T = 140 was cancelled; three blocks are reported.", "")
for (t in TBL) lines <- c(lines, sprintf("%% ---- %s ----", t[[1]]),
                          rewrite(t[[1]], t[[2]], t[[3]], t[[4]]), "")
writeLines(lines, OUT)
cat(sprintf("wrote %s (%d lines)\n", OUT, length(lines)))
DD <- if (length(DIFFS)) do.call(rbind, DIFFS) else NULL
cat(sprintf("\n=== cells differing from the archived sim_tables.tex: %d ===\n",
            if (is.null(DD)) 0L else nrow(DD)))
if (!is.null(DD)){
   ct <- as.data.frame(table(DD$table)); names(ct) <- c("table","cells_differing")
   print(ct, row.names = FALSE)
   cat("\n")
   print(DD, row.names = FALSE) }
saveRDS(list(diffs = DD), file.path(SIM_FITS, "sim_tables_diff.rds"))
