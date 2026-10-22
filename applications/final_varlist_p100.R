# ============================================================================
#  Item 9: LaTeX longtable of the series in the p = 100 graph, in the format
#  of Table S8 (tab:varlist).  Descriptions come from varlist_descriptions.csv,
#  table where the series already appears; any series new to the graph is
#  flagged so a description can be written for it.
#     Rscript applications/varlist_table.R
#  Writes results/tables/varlist_p100.tex and prints the table.
# ============================================================================
suppressMessages(source("/Users/peterdunson/Desktop/vapor/applications/final_shared.R"))
options(width = 200)
#  Descriptions, transformation codes and group codes for the series in the
#  graph.  These used to be scraped out of Table S8 in main-77.tex, which made
#  the paper source a build dependency; they now live in a small CSV beside
#  fredqd_groups.csv so this script depends on nothing outside the repo.
DESCCSV <- file.path(dirname(RES_DIR), "varlist_descriptions.csv")
GROUPS  <- file.path(dirname(RES_DIR), "fredqd_groups.csv")
CSV     <- "/Users/peterdunson/Desktop/vapor/applications/data/2026-05-qd.csv"
OUTTEX  <- file.path(TAB_DIR, "varlist_p100.tex")

R  <- readRDS(file.path(FITS_DIR, "p100_final.rds"))
ed <- R$long$ci_edges
ser <- sort(unique(c(ed$v1, ed$v2)))
cat(sprintf("%d edges, %d distinct series\n", nrow(ed), length(ser)))

DD <- read.csv(DESCCSV, stringsAsFactors = FALSE)
desc <- setNames(as.list(DD$description), DD$series)
cat(sprintf("descriptions available: %d (%s)\n", nrow(DD), basename(DESCCSV)))

#  transformation codes straight from the FRED-QD file
raw <- read.csv(path.expand(CSV), nrows = 2)
tc  <- as.numeric(raw[2, -1]); names(tc) <- colnames(raw)[-1]
g   <- read.csv(GROUPS, stringsAsFactors = FALSE)
gc  <- setNames(g$group_code, g$series); gn <- setNames(g$group_name, g$series)

D <- data.frame(Mnemonic = ser,
                Description = vapply(ser, function(s)
                   if (!is.null(desc[[s]])) desc[[s]] else NA_character_,
                   character(1)),
                T = as.integer(tc[ser]), G = as.integer(gc[ser]),
                Group = unname(gn[ser]), stringsAsFactors = FALSE, row.names = NULL)
D$source <- basename(DESCCSV)
miss <- D$Mnemonic[is.na(D$Description)]
#  the CSV also carries T and G; cross-check them against the live sources so
#  a stale CSV cannot silently ship wrong codes
chk <- merge(D[, c("Mnemonic","T","G")], DD, by.x = "Mnemonic", by.y = "series")
badT <- chk$Mnemonic[chk$T != chk$tcode]; badG <- chk$Mnemonic[chk$G != chk$group_code]
if (length(badT) || length(badG))
   stop(sprintf("varlist_descriptions.csv disagrees with the live sources: T %s | G %s",
                paste(badT, collapse=","), paste(badG, collapse=",")))
cat(sprintf("T and G cross-checked against the FRED-QD file and fredqd_groups.csv: %d series OK\n",
            nrow(chk)))

#  a missing description must stop the run: a placeholder in the paper's
#  variable table is worse than no table at all
if (length(miss))
   stop(sprintf("no description for %d series: %s -- add them to NEWDESC",
                length(miss), paste(miss, collapse = ", ")))
print(D, row.names = FALSE)

#  group legend, in the order the codes appear
lg <- unique(D[order(D$G), c("G","Group")])
legend <- paste(sprintf("%d = %s", lg$G, lg$Group), collapse = ", ")

#  two series per row, as in Table S8
n <- nrow(D); half <- ceiling(n/2)
L <- D[seq_len(half), ]; Rt <- D[(half+1):n, ]
while (nrow(Rt) < half){
   pad <- Rt[1, ]; pad[] <- NA; pad$Mnemonic <- ""; pad$Description <- ""
   pad$Group <- ""; pad$source <- ""; Rt <- rbind(Rt, pad) }
esc <- function(x) gsub("&", "\\\\&", x)
rows <- vapply(seq_len(half), function(i)
   sprintf("%s & %s & %s & %s & %s & %s & %s & %s \\\\",
           esc(L$Mnemonic[i]), esc(L$Description[i]),
           ifelse(is.na(L$T[i]), "", L$T[i]), ifelse(is.na(L$G[i]), "", L$G[i]),
           esc(Rt$Mnemonic[i]), esc(Rt$Description[i]),
           ifelse(is.na(Rt$T[i]), "", Rt$T[i]), ifelse(is.na(Rt$G[i]), "", Rt$G[i])),
   character(1))

tex <- c("{\\footnotesize",
   "\\setlength{\\tabcolsep}{4pt}",
   "\\begin{longtable}{l>{\\raggedright\\arraybackslash}p{3.4cm}cc@{\\hskip 1.2em}l>{\\raggedright\\arraybackslash}p{3.4cm}cc}",
   sprintf("\\caption{The %d FRED-QD series appearing in the reported $p = 100$ dependence graph (Figure~\\ref{fig:network_p100}). T is the McCracken--Ng transformation code: 1 = no transformation, 2 = first difference, 5 = first difference of logarithms, 6 = second difference of logarithms. G is the FRED-QD group: %s.}",
           n, legend),
   "\\label{tab:varlist}\\\\", "\\toprule",
   "Mnemonic & Description & T & G & Mnemonic & Description & T & G \\\\",
   "\\midrule", "\\endfirsthead",
   "\\multicolumn{8}{c}{\\tablename\\ \\thetable\\ (continued)}\\\\", "\\toprule",
   "Mnemonic & Description & T & G & Mnemonic & Description & T & G \\\\",
   "\\midrule", "\\endhead", "\\bottomrule", "\\endfoot",
   rows, "\\end{longtable}", "}")
writeLines(tex, OUTTEX)
cat(sprintf("\nwrote %s (%d series, %d body rows)\n", OUTTEX, n, half))
saveRDS(list(table = D, legend = legend, missing = miss, tex = tex),
        file.path(FITS_DIR, "varlist_p100.rds"))
