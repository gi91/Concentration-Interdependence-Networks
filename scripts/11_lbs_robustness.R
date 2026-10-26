# =============================================================================
# 11_lbs_robustness.R -- Online Appendix: banking network on the BIS
# Locational Banking Statistics (residence basis, immediate counterparty,
# cross-border claims of all banks located in each reporting country).
# Portugal and Turkey do not publish bilateral LBS claims, so the comparison
# uses the 17 common reporters and recomputes the CBS network on the same set.
# Output: tabOA_lbs.tex/csv, figOA_lbs.pdf/png
# =============================================================================
source("R/utils.R"); source("R/core_indices.R"); source("R/data_loaders.R")
source("R/figures.R")

THETA <- 0.02; H <- 3
qlab <- function(d) paste0(format(d, "%Y"), "-Q", (as.integer(format(d, "%m")) - 1) %/% 3 + 1)
TABLE_Q <- as.Date(c("2007-01-01", "2009-10-01", "2012-10-01", "2019-10-01", "2025-10-01"))

lbs_cl <- readRDS(file.path(DIR_PROC, "lbs_claims.rds"))
lbs_tt <- readRDS(file.path(DIR_PROC, "lbs_totals.rds"))
# a few LBS cells are negative (reporting adjustments): treated as no exposure
cat("negative LBS bilateral cells set to zero:", sum(lbs_cl$claims < 0), "of", nrow(lbs_cl), "\n")
lbs_cl[claims < 0, claims := 0]
cbs_cl <- readRDS(file.path(DIR_PROC, "bis_claims.rds"))
cbs_tt <- readRDS(file.path(DIR_PROC, "bis_totals.rds"))
REP17 <- sort(intersect(unique(lbs_cl$creditor), unique(cbs_tt$creditor)))
cat("common reporters:", length(REP17), paste(REP17, collapse = " "), "\n")

measures <- function(cl, tt, q) {
  sn <- bis_snapshot(cl, tt, q, REP17)
  if (anyNA(sn$total)) return(NULL)
  A <- bis_adjacency(sn, THETA); w <- sn$w; N <- length(w)
  b <- baseline_index(w, A)
  sd <- structural_decomposition(w, A, H); tb <- sd$table
  sw <- structural_decomposition(w, bis_adjacency(sn, weighted = TRUE), H)
  g <- function(t, c, h) t$Gamma[t$comp == c & t$h == h]
  data.frame(quarter = q, hhi = hhi(w), density = sum(A) / (N * (N - 1)),
             psiD = b$psiD, Gamma1 = b$Gamma1,
             psiF1 = tb$psi[tb$comp == "F" & tb$h == 1],
             share_rec = sd$shares[["reciprocal"]], G_rec = g(tb, "F", 1),
             G_c3 = g(tb, "F", 2), Gw_rec = g(sw$table, "F", 1))
}
QS <- as.Date(sort(intersect(unique(lbs_tt$quarter), unique(cbs_tt$quarter))), origin = "1970-01-01")
QS <- QS[QS >= as.Date("2007-01-01")]
res <- rbind(
  cbind(source = "CBS, ultimate risk", do.call(rbind, lapply(as.list(QS), function(q) measures(cbs_cl, cbs_tt, q)))),
  cbind(source = "LBS, residence", do.call(rbind, lapply(as.list(QS), function(q) measures(lbs_cl, lbs_tt, q)))))
saveRDS(res, file.path(DIR_PROC, "results_lbs.rds"))

# ---- table -----------------------------------------------------------------------
f3 <- function(x) fmt_num(x, 3); f2 <- function(x) fmt_num(x, 2)
vars <- c(hhi = "Herfindahl index of $w$", density = "Density", psiD = "$\\psi_D(w)$",
          Gamma1 = "$\\Gamma_1(w)$", share_rec = "Share of reciprocal exposures in $\\psi_{TOT}$",
          G_rec = "$\\Gamma$ of reciprocal exposures", G_c3 = "$\\Gamma$ of circuits of length 3",
          Gw_rec = "Weighted: $\\Gamma$ of reciprocal exposures")
tab <- do.call(rbind, lapply(unique(res$source), function(s) {
  r <- res[res$source == s, ]; r <- r[match(TABLE_Q, r$quarter), ]
  data.frame(panel = s, Measure = unname(vars),
             t(sapply(names(vars), function(v) (if (v %in% c("hhi", "density", "psiD")) f3 else f2)(r[[v]]))),
             check.names = FALSE)
}))
names(tab)[-(1:2)] <- qlab(TABLE_Q)
write_tex_table(tab[, -1], "tabOA_lbs.tex", panel = tab$panel,
  caption = "International banking network: consolidated versus locational statistics",
  label = "tab:oa-lbs", align = "lccccc", header = c("Measure", qlab(TABLE_Q)),
  notes = paste0("Common sample of ", length(REP17), " reporting countries (Portugal and",
    " Turkey do not publish bilateral locational claims). CBS: foreign claims of",
    " domestic banks on an ultimate-risk basis. LBS: cross-border claims of all banks",
    " located in the country, immediate counterparty, all sectors. Threshold",
    " $\\theta=2\\%$; structural indices with $H=3$ and $\\alpha_h=1/3$."))

# ---- figure -------------------------------------------------------------------------
res$x <- dec_year(res$quarter)
L <- rbind(data.frame(x = res$x, v = res$G_rec, src = res$source, m = "Binary, θ = 2%"),
           data.frame(x = res$x, v = res$Gw_rec, src = res$source, m = "Weighted"))
COL_SRC <- c("CBS, ultimate risk" = COL_MECH[["Reciprocal"]], "LBS, residence" = INK)
p <- ggplot(L, aes(x, v, colour = src)) + geom_crisis() +
  geom_hline(yintercept = 1, colour = MUTED, linewidth = .3, linetype = "22") +
  geom_line(linewidth = .6) + facet_wrap(~ m, nrow = 1) +
  scale_colour_manual(values = COL_SRC) +
  scale_x_continuous(breaks = data_breaks(L$x)) +
  labs(x = NULL, y = expression(Gamma*" of reciprocal exposures"),
       title = "Alignment of reciprocal exposures: consolidated vs locational statistics",
       subtitle = sprintf("Same %d reporting countries; 1 = no alignment with size", length(REP17))) +
  theme_paper()
save_fig(p, "figOA_lbs", width = 7, height = 3.6)
mm <- merge(res[res$source == "CBS, ultimate risk", ], res[res$source == "LBS, residence", ],
            by = "quarter", suffixes = c(".cbs", ".lbs"))
cr <- sapply(c("G_rec", "Gw_rec", "psiD", "share_rec"), function(v)
  cor(mm[[paste0(v, ".cbs")]], mm[[paste0(v, ".lbs")]]))
cat("correlation CBS vs LBS over", nrow(mm), "common quarters:\n")
print(round(cr, 2))
log_msg("LBS robustness written")
