# =============================================================================
# 03_bis.R -- Section 7.2: international banking network, Table 3, Figure 3
# =============================================================================
source("R/utils.R"); source("R/core_indices.R"); source("R/data_loaders.R")
source("R/null_models.R"); source("R/fixes.R"); source("R/figures.R")

claims <- readRDS(file.path(DIR_PROC, "bis_claims.rds"))
totals <- readRDS(file.path(DIR_PROC, "bis_totals.rds"))
REPORTERS <- sort(unique(totals$creditor))     # the 19 continuous reporters
QUARTERS  <- sort(unique(totals$quarter[totals$quarter >= as.Date("2007-01-01")]))  # sample of the paper
THETA <- 0.02; THETAS <- c(0.01, 0.02, 0.05)
H <- 3; ALPHA <- rep(1 / 3, 3)
TABLE_Q <- as.Date(c("2007-01-01", "2009-10-01", "2012-10-01", "2019-10-01",
                     "2025-10-01"))
DEST <- c(US = "United States", GB = "United Kingdom", DE = "Germany",
          FR = "France", JP = "Japan")
qlab <- function(d) paste0(format(d, "%Y"), "-Q", (as.integer(format(d, "%m")) - 1) %/% 3 + 1)

gam <- function(tb, comp, h) tb$Gamma[tb$comp == comp & tb$h == h]
psv <- function(tb, comp, h) tb$psi[tb$comp == comp & tb$h == h]

rows <- list(); mats <- list()
for (q in as.list(QUARTERS)) {
  sn <- bis_snapshot(claims, totals, q, REPORTERS)
  w  <- sn$w; A <- bis_adjacency(sn, THETA)
  b  <- baseline_index(w, A)
  tp <- topological_profile(w, A)
  sd <- structural_decomposition(w, A, H, ALPHA); tb <- sd$table
  sw <- structural_decomposition(w, bis_adjacency(sn, weighted = TRUE), H, ALPHA)
  rob <- sapply(THETAS, function(th) {
    t2 <- structural_decomposition(w, bis_adjacency(sn, th), 1, 1)$table
    gam(t2, "F", 1) })
  # propagation: column r of P = distribution of r's claims over countries
  pm <- propagation_measures(w, t(sn$C), 1)
  land <- setNames(pm$wk, REPORTERS)
  rows[[qlab(q)]] <- data.frame(
    quarter = q, hhi = hhi(w), density = sum(A) / (19 * 18),
    psiD = b$psiD, Gamma1 = b$Gamma1, dw = tp$dw, dbar = tp$dbar,
    psiR1 = psv(tb, "R", 1), psiF1 = psv(tb, "F", 1),
    share_rec = sd$shares[["reciprocal"]],
    G_rec = gam(tb, "F", 1), G_c3 = gam(tb, "F", 2), G_c4 = gam(tb, "F", 3),
    Gw_rec = gam(sw$table, "F", 1), Gw_c3 = gam(sw$table, "F", 2),
    Gw_c4 = gam(sw$table, "F", 3),
    G_rec_t1 = rob[1], G_rec_t2 = rob[2], G_rec_t5 = rob[3],
    psiO = pm$psiO, GammaO = pm$GammaO, share_conv = pm$share_conv,
    t(land[names(DEST)]))
  mats[[qlab(q)]] <- list(w = w, A = A, C = sn$C, total = sn$total,
                          R1 = sd$mats$R1, F1 = sd$mats$F1,
                          F1w = sw$mats$F1, R1w = sw$mats$R1,
                          S = bis_adjacency(sn, weighted = TRUE),
                          within = within_order_shares(sd))
}
ts <- do.call(rbind, rows)
saveRDS(list(ts = ts, mats = mats, reporters = REPORTERS),
        file.path(DIR_PROC, "results_bis.rds"))

# ---- Table 3 ---------------------------------------------------------------------
ty <- ts[match(TABLE_Q, ts$quarter), ]
f3 <- function(x) fmt_num(x, 3); f2 <- function(x) fmt_num(x, 2)
pr <- function(a, b, f = f2) paste(f(a), f(b), sep = "; ")
mk <- function(lab, v) setNames(data.frame(lab, t(v)), c("Measure", qlab(TABLE_Q)))
tab3 <- rbind(
  mk("Herfindahl index of $w$", f3(ty$hhi)),
  mk("Density", f3(ty$density)),
  mk("$\\psi_D(w)$", f3(ty$psiD)),
  mk("$\\Gamma_1(w)$", f2(ty$Gamma1)),
  mk("Economic mean distance $\\bar d_w$; average path length $\\bar d$", pr(ty$dw, ty$dbar)),
  mk("Reciprocal component $\\psi^{(1)}_F(w)$", f3(ty$psiF1)),
  mk("Share of reciprocal exposures in $\\psi_{TOT}$", f2(ty$share_rec)),
  mk("$\\Gamma$ of reciprocal exposures", f2(ty$G_rec)),
  mk("$\\Gamma$ of circuits of length 3; of length 4", pr(ty$G_c3, ty$G_c4)),
  mk("Weighted: $\\Gamma$ of reciprocal exposures", f2(ty$Gw_rec)),
  mk("Weighted: $\\Gamma$ of circuits of length 3; of length 4", pr(ty$Gw_c3, ty$Gw_c4)),
  mk("Portfolio overlap $\\psi^{(1)}_O(w)$; $\\Gamma_{O^{(1)}}(w)$",
     paste(f3(ty$psiO), f2(ty$GammaO), sep = "; ")),
  mk("Share of $w^{(1)}$ landing in the United States; United Kingdom", pr(ty$US, ty$GB)))
write_tex_table(tab3, "tab3_bis.tex",
  caption = "Concentration of the international banking network, selected quarters",
  label = "tab:bis", align = "lccccc", header = c("Measure", qlab(TABLE_Q)),
  notes = paste0("BIS Consolidated Banking Statistics, foreign claims on an",
    " ultimate-risk basis, ", length(REPORTERS), " reporting countries weighted by",
    " their foreign claims. Threshold $\\theta=", 100 * THETA, "\\%$; weighted",
    " rows use the exposure shares as link intensities, without threshold.",
    " Structural indices with $H=3$ and $\\alpha_h=1/3$."))

# threshold robustness (text of Section 7.2)
rob_tab <- data.frame(Threshold = c("1\\%", "2\\%", "5\\%"),
  t(sapply(c("G_rec_t1", "G_rec_t2", "G_rec_t5"), function(v) f2(ty[[v]]))))
names(rob_tab)[-1] <- qlab(TABLE_Q)
write_tex_table(rob_tab, "tabOA_bis_threshold.tex",
  caption = "Alignment factor of reciprocal exposures for alternative thresholds",
  label = "tab:oa-bis", align = "lccccc", header = c("$\\theta$", qlab(TABLE_Q)))

# ---- decompositions (Proposition 5) ---------------------------------------------
dec_between <- function(q0, q1, R = "R1", F = "F1") {
  a <- mats[[q0]]; b <- mats[[q1]]
  rbind(cbind(comp = "OW", period = paste(q0, q1, sep = "/"),
              dynamic_decomposition(a$w, a[[R]], b$w, b[[R]])),
        cbind(comp = "RE", period = paste(q0, q1, sep = "/"),
              dynamic_decomposition(a$w, a[[F]], b$w, b[[F]])))
}
dec  <- rbind(dec_between("2007-Q1", "2013-Q4"), dec_between("2013-Q4", "2025-Q4"))
decG <- dec_between("2007-Q2", "2009-Q2")                          # GFC only
decW <- dec_between("2007-Q1", "2013-Q4", "R1w", "F1w")           # weighted
write.csv(rbind(cbind(spec = "binary", dec), cbind(spec = "binary_GFC", decG),
                cbind(spec = "weighted", decW)),
          file.path(DIR_TAB, "fig3c_decomposition_bis.csv"), row.names = FALSE)

# ---- κ-free within-order shares ---------------------------------------------
wos <- do.call(rbind, lapply(qlab(TABLE_Q), function(q) cbind(quarter = q, mats[[q]]$within)))
write_tex_table(data.frame(Order = sprintf("$h=%d$", wos$h), Quarter = wos$quarter,
    Closed = fmt_num(wos$closed_share, 3), Gopen = fmt_num(wos$Gamma_open, 2),
    Gclosed = fmt_num(wos$Gamma_closed, 2))[order(wos$h), ],
  "tabOA_bis_within_order.tex",
  caption = "International banking network: within-order structural decomposition",
  label = "tab:bis-within", align = "lcccc",
  header = c("Order", "Quarter", "Closed share of $\\psi^{(h)}_T$", "$\\Gamma$ open", "$\\Gamma$ closed"))

# ---- permutation band for Gamma of reciprocal exposures -------------------------
set.seed(12)
ts$x <- dec_year(ts$quarter)
band <- do.call(rbind, lapply(seq_len(nrow(ts)), function(k) {
  m <- mats[[qlab(ts$quarter[k])]]
  pa <- perm_alignment(m$w, m$F1, B = 300)
  data.frame(x = ts$x[k], lo = quantile(pa$null, .05), hi = quantile(pa$null, .95))
}))

# ---- Figure 3 --------------------------------------------------------------------
pa <- fig_levels(ts, "x", c("hhi", "density", "psiD", "psiF1"),
  c("HHI of w", "Density", "Baseline", "Reciprocal"),
  c(COL_HHI, MUTED, INK, COL_MECH[["Reciprocal"]]), c("22", "solid", "solid", "solid"),
  title = "(a) Levels: stable baseline, thinning network",
  subtitle = expression("Baseline "*psi[D]*", directed density, reciprocal "*psi[F]^{(1)}),
  nudge = 1, crisis = TRUE, ylim = c(0, 1))
pb <- fig_alignment(ts, "x", band, "(b) Alignment with banking-system size",
  "Shaded: 5–95% of reciprocal Γ under random weights",
  nudge = 1, crisis = TRUE)
dec <- rbind(dec_between("2007-Q1", "2013-Q4"), dec_between("2013-Q4", "2025-Q4"))
dec$period <- sub("/", "–", gsub("20(\\d\\d)-Q(\\d)", "\\1Q\\2", dec$period))
pc <- fig_decomp(dec, "(c) Sources of change (Proposition 5)")
ld <- to_long(ts, "x", names(DEST), DEST); ld$s <- factor(ld$s, DEST)
COL_DEST <- setNames(CAT[1:5], DEST)
pd <- ggplot(ld, aes(x, v, colour = s)) + geom_crisis() + geom_line(linewidth = .6) +
  end_labels(ld, "x", "v", "s", COL_DEST, nudge = 1) +
  scale_colour_manual(values = COL_DEST) +
  scale_x_continuous(breaks = data_breaks(ld$x), expand = expansion(mult = c(.02, .32))) +
  scale_y_continuous(labels = scales::percent) + coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL, title = "(d) Where size-weighted exposures land",
       subtitle = "Share of one-step propagated exposures by counterparty country") +
  theme_paper() + theme(legend.position = "none")
save_fig((pa | pb) / (pc | pd), "fig3_bis", width = 7.5, height = 7.4)

# ---- Figure 5b: pairs driving the alignment of reciprocal exposures -------------
PQ <- c("2007-Q1", "2012-Q4", "2025-Q4")
pcs <- do.call(rbind, lapply(PQ, function(q) {
  m <- mats[[q]]; p <- pair_contributions(m$w, m$F1, labels = REPORTERS, top = 8)
  G <- attr(p, "Gamma")
  cbind(p, date = sprintf("%s (Γ = %.2f)", q, G), share = p$contribution / G)
}))
write.csv(pcs, file.path(DIR_TAB, "fig5_pairs_bis.csv"), row.names = FALSE)
save_fig(fig_pairs(pcs, "Pairs driving the alignment of reciprocal exposures",
  "Top-8 pairs of banking systems; labels: share of Γ"), "fig5b_pairs_bis",
  width = 7.5, height = 3.2)

# ---- Figure OA: threshold profile ---------------------------------------------
THG <- c(0.01, 0.02, 0.05)
tpl <- do.call(rbind, lapply(seq_len(nrow(ts)), function(k) {
  m <- mats[[qlab(ts$quarter[k])]]
  rbind(data.frame(x = ts$x[k], v = gamma_threshold_profile(m$w, m$S, THG),
                   s = sprintf("Binary, %g%%", 100 * THG)),
        data.frame(x = ts$x[k], v = ts$Gw_rec[k], s = "Weighted"))
}))
save_fig(fig_threshold(tpl, "x", "Threshold dependence of the reciprocal alignment",
  "Binary networks at three thresholds (blue, light to dark) and weighted network (orange)",
  crisis = TRUE, nudge = 1), "figOA_threshold_bis", width = 6, height = 3.6)
log_msg("Table 3 and Figures 3, 5b, OA written")
