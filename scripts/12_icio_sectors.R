# =============================================================================
# 12_icio_sectors.R -- the production network at the level of world industries
# (OECD ICIO, 50 industries, 1995-2022), as a second level of aggregation of
# the same data used for countries in 02_wiod.R.
#
# Two networks of industries:
#   "All flows"    Z^s[a, b] = inputs supplied by industry a to industry b,
#                  summed over all supplier and user countries, domestic flows
#                  included (the production network of Acemoglu et al., 2012)
#   "Cross-border" only flows between different countries (the trade network
#                  of intermediate products)
# Weights: shares of world value added.  Link a -> b when b sources from a at
# least a share theta of its gross output (A[b, a] = 1), as for countries.
# All flows: theta = 2% (density close to the country network); cross-border:
# theta chosen each year to give the same density (inputs crossing borders
# are an order of magnitude smaller).  Weighted networks need no threshold.
# Output: tabOA_sectors.tex/csv, fig6_levels.pdf/png, fig6b_pairs_sectors.pdf/png
# =============================================================================
source("R/utils.R"); source("R/core_indices.R"); source("R/data_loaders.R")
source("R/fixes.R"); source("R/figures.R")

THETA_ALL <- 0.02
YEARS <- sort(as.integer(sub("^icio_([0-9]{4})\\.rds$", "\\1",
  list.files(DIR_PROC, "^icio_[0-9]{4}\\.rds$"))))
TABLE_YEARS <- intersect(c(1995, 2001, 2007, 2014, 2019, 2022), YEARS)
B_DB <- as.integer(Sys.getenv("NCI_B_SECT", 200))   # DBCM draws for the closure benchmark

# ISIC Rev.4 labels of the 50 ICIO industries (short)
IND_LAB <- c(A01 = "Agriculture", A02 = "Forestry", A03 = "Fishing", B05 = "Coal",
  B06 = "Oil & gas", B07 = "Metal ores", B08 = "Other mining", B09 = "Mining services",
  C10T12 = "Food", C13T15 = "Textiles", C16 = "Wood", C17_18 = "Paper & printing",
  C19 = "Refined petroleum", C20 = "Chemicals", C21 = "Pharmaceuticals", C22 = "Rubber & plastics",
  C23 = "Non-metallic minerals", C24A = "Iron & steel", C24B = "Non-ferrous metals",
  C25 = "Metal products", C26 = "Electronics", C27 = "Electrical equipment", C28 = "Machinery",
  C29 = "Motor vehicles", C301 = "Ships", C302T309 = "Other transport eq.", C31T33 = "Other manuf.",
  D = "Electricity", E = "Water & waste", F = "Construction", G = "Trade", H49 = "Land transport",
  H50 = "Water transport", H51 = "Air transport", H52 = "Warehousing", H53 = "Postal",
  I = "Hotels & restaurants", J58T60 = "Publishing", J61 = "Telecom", J62_63 = "IT services",
  K = "Finance", L = "Real estate", M = "Professional services", N = "Admin. services",
  O = "Public admin.", P = "Education", Q = "Health", R = "Arts", S = "Other services",
  T = "Households")

sector_network <- function(d, cross_border) {
  Z <- d$Z
  if (cross_border) Z[outer(d$country, d$country, "==")] <- 0
  Zs <- t(rowsum(t(rowsum(Z, d$industry, reorder = TRUE)), d$industry, reorder = TRUE))
  X  <- as.numeric(rowsum(d$go, d$industry, reorder = TRUE))
  VA <- as.numeric(rowsum(d$va, d$industry, reorder = TRUE))
  S  <- t(Zs) / X                            # S[b, a] = Z_ab / X_b
  diag(S) <- 0; S[!is.finite(S)] <- 0
  keep <- X > 0 & VA > 0
  list(S = S[keep, keep], w = pmax(VA[keep], 0) / sum(pmax(VA[keep], 0)),
       ind = rownames(Zs)[keep])
}

measure <- function(net, theta) {
  A <- (net$S >= theta) * 1; w <- net$w; N <- length(w)
  b  <- baseline_index(w, A)
  sd <- structural_decomposition(w, A, H = 3); tb <- sd$table
  sw <- structural_decomposition(w, pmin(net$S / max(net$S), 1), H = 3)
  g  <- function(t, c, h) t$Gamma[t$comp == c & t$h == h]
  wo <- within_order_shares(sd)
  # closed share at h = 2 relative to its degree-conditional expectation (Fig. S6)
  fit <- dbcm_fit(A)
  cs_null <- mean(replicate(B_DB, within_order_shares(structural_decomposition(
    w, dbcm_sample(fit), H = 2, alpha = c(.5, .5)))$closed_share[2]), na.rm = TRUE)
  list(row = data.frame(theta = theta, N = N, hhi = hhi(w), density = sum(A) / (N * (N - 1)),
         psiD = b$psiD, Gamma1 = b$Gamma1, psiF1 = tb$psi[tb$comp == "F" & tb$h == 1],
         share_rec = sd$shares[["reciprocal"]], G_rec = g(tb, "F", 1), G_c3 = g(tb, "F", 2),
         G_c4 = g(tb, "F", 3), Gw_rec = g(sw$table, "F", 1), Gw_c3 = g(sw$table, "F", 2),
         closed_h2 = wo$closed_share[2], closed_h2_rel = wo$closed_share[2] / cs_null),
       F1 = sd$mats$F1, F1w = sw$mats$F1)
}

set.seed(12)
rows <- list(); mats <- list()
for (yr in YEARS) {
  d <- load_wiod_year(yr, "icio")
  na <- sector_network(d, FALSE); nx <- sector_network(d, TRUE)
  ma <- measure(na, THETA_ALL)
  th_x <- threshold_for_density(nx$S, ma$row$density)   # same density as all flows
  mx <- measure(nx, th_x)
  rows[[length(rows) + 1]] <- rbind(cbind(year = yr, net = "All flows", ma$row),
                                    cbind(year = yr, net = "Cross-border", mx$row))
  mats[[as.character(yr)]] <- list(all = c(na, list(F1 = ma$F1, F1w = ma$F1w)),
                                   xb = c(nx, list(F1 = mx$F1, F1w = mx$F1w)))
  log_msg("sectors ", yr)
}
ts <- do.call(rbind, rows)
saveRDS(list(ts = ts, mats = mats), file.path(DIR_PROC, "results_icio_sectors.rds"))

# ---- table -------------------------------------------------------------------------
f3 <- function(x) fmt_num(x, 3); f2 <- function(x) fmt_num(x, 2)
vars <- c(hhi = "Herfindahl index of $w$", density = "Density", psiD = "$\\psi_D(w)$",
  Gamma1 = "$\\Gamma_1(w)$", share_rec = "Share of reciprocal links in $\\psi_{TOT}$",
  G_rec = "$\\Gamma$ of reciprocal links", G_c3 = "$\\Gamma$ of circuits of length 3",
  Gw_rec = "Weighted: $\\Gamma$ of reciprocal links", Gw_c3 = "Weighted: $\\Gamma$ of circuits of length 3",
  closed_h2_rel = "Closed share ($h=2$) / degree-conditional expectation")
tab <- do.call(rbind, lapply(c("All flows", "Cross-border"), function(nt) {
  r <- ts[ts$net == nt, ]; r <- r[match(TABLE_YEARS, r$year), ]
  data.frame(panel = nt, Measure = unname(vars),
    t(sapply(names(vars), function(v) (if (v %in% c("hhi", "density", "psiD")) f3 else f2)(r[[v]]))),
    check.names = FALSE)
}))
names(tab)[-(1:2)] <- TABLE_YEARS
write_tex_table(tab[, -1], "tabOA_sectors.tex", panel = tab$panel,
  caption = "Concentration of the world production network at the level of industries",
  label = "tab:oa-sectors", align = paste0("l", strrep("c", length(TABLE_YEARS))),
  header = c("Measure", TABLE_YEARS),
  notes = paste0("OECD ICIO, 2025 edition; 50 world industries weighted by world value added.",
    " All flows: inputs summed over supplier and user countries, threshold $\\theta=2\\%$ of",
    " the user's gross output. Cross-border: flows between different countries only,",
    " threshold chosen each year to give the same density. Weighted rows: input shares",
    " without threshold. Last row: closed share of $\\psi^{(2)}_T$ divided by its mean over ",
    B_DB, " draws of the directed configuration model."))

# ---- figure 6: what survives aggregation ------------------------------------------
rc <- readRDS(file.path(DIR_PROC, "results_icio.rds"))$ts
L <- rbind(
  data.frame(year = rc$year, v = rc$G_rec,  lvl = "Countries (80)", spec = "Binary"),
  data.frame(year = rc$year, v = rc$Gw_rec, lvl = "Countries (80)", spec = "Weighted"),
  data.frame(year = ts$year, v = ts$G_rec,  lvl = paste0("Industries (50), ", tolower(ts$net)), spec = "Binary"),
  data.frame(year = ts$year, v = ts$Gw_rec, lvl = paste0("Industries (50), ", tolower(ts$net)), spec = "Weighted"))
L$lvl <- factor(L$lvl, c("Countries (80)", "Industries (50), all flows", "Industries (50), cross-border"))
COL_LVL <- setNames(c(COL_MECH[["Reciprocal"]], CAT[1], CAT[3]), levels(L$lvl))
p6 <- ggplot(L, aes(year, v, colour = lvl)) +
  geom_hline(yintercept = 1, colour = MUTED, linewidth = .3, linetype = "22") +
  geom_line(linewidth = .65) + facet_wrap(~ spec, nrow = 1) +
  scale_colour_manual(values = COL_LVL) + scale_y_log10(breaks = c(.5, 1, 2, 4, 8)) +
  labs(x = NULL, y = expression(Gamma*" of reciprocal links (log scale)"),
       title = "Alignment of reciprocal links with size: countries and industries",
       subtitle = "Same ICIO data aggregated by country or by world industry; 1 = no alignment") +
  theme_paper()
save_fig(p6, "fig6_levels", width = 7.2, height = 3.8)

# industry pairs driving the alignment of cross-border reciprocal links (binary)
pcs <- do.call(rbind, lapply(TABLE_YEARS[c(1, 3, length(TABLE_YEARS))], function(y) {
  m <- mats[[as.character(y)]]$xb
  lab <- unname(IND_LAB[m$ind]); lab[is.na(lab)] <- m$ind[is.na(lab)]
  p <- pair_contributions(m$w, m$F1, labels = lab, top = 6); G <- attr(p, "Gamma")
  cbind(p, date = sprintf("%d (Γ = %.2f)", y, G), share = p$contribution / G)
}))
write.csv(pcs, file.path(DIR_TAB, "fig6b_pairs_sectors.csv"), row.names = FALSE)
save_fig(fig_pairs(pcs, "Industry pairs driving the alignment of cross-border reciprocal links",
  "World industries, flows between different countries; top-6 pairs, labels: share of Γ", nrow = 3),
  "fig6b_pairs_sectors", width = 6, height = 6)
log_msg("sector-level analysis written")
