# =============================================================================
# 02_wiod.R -- Section 7.1: world production network, Table 2, Figure 2,
#              Online Appendix OA.1 (threshold / weighted / size-neutral)
# =============================================================================
source("R/utils.R"); source("R/core_indices.R"); source("R/data_loaders.R")
source("R/null_models.R"); source("R/fixes.R"); source("R/figures.R")

# years available for the chosen source (NCI_PROD = "wiod" or "icio")
YEARS  <- sort(as.integer(sub(paste0("^", PROD_SOURCE, "_([0-9]{4})\\.rds$"), "\\1",
  list.files(DIR_PROC, paste0("^", PROD_SOURCE, "_[0-9]{4}\\.rds$")))))
stopifnot(length(YEARS) > 1)
SRC <- PROD_SOURCE                      # suffix of every output file
THETA  <- 0.005                     # baseline threshold
THETAS <- c(0.002, 0.005, 0.01)     # robustness
H <- 3; ALPHA <- rep(1 / 3, 3)
KMAX <- 3                           # upstream steps for propagation
TABLE_YEARS <- intersect(if (SRC == "icio") c(1995, 2001, 2007, 2014, 2019, 2022)
                         else c(2000, 2007, 2014), YEARS)
if (length(TABLE_YEARS) < 2) TABLE_YEARS <- YEARS[unique(round(seq(1, length(YEARS), length.out = 3)))]
ALIGN <- function(lead) paste0(lead, strrep("c", length(TABLE_YEARS)))

gam <- function(tb, comp, h) tb$Gamma[tb$comp == comp & tb$h == h]
psv <- function(tb, comp, h) tb$psi[tb$comp == comp & tb$h == h]

country_measures <- function(dc, w, theta) {
  S <- wiod_intensity(dc, "share")
  A <- (S >= theta) * 1
  b  <- baseline_index(w, A)
  tp <- topological_profile(w, A)
  sd <- structural_decomposition(w, A, H, ALPHA)
  sw <- structural_decomposition(w, pmin(S, 1), H, ALPHA)     # weighted
  # size-neutral binary network with the same number of links
  Sn <- wiod_intensity(dc, "neutral")
  thr_n <- sort(Sn[row(Sn) != col(Sn)], decreasing = TRUE)[sum(A)]
  An <- (Sn >= thr_n) * 1
  sn <- structural_decomposition(w, An, H, ALPHA)
  tb <- sd$table
  list(
    row = data.frame(
      hhi = hhi(w), density = sum(A) / (nrow(A) * (nrow(A) - 1)),
      psiD = b$psiD, Gamma1 = b$Gamma1,
      psiG2 = tp$profile$psi[2], Gamma2 = tp$profile$Gamma[2],
      dw = tp$dw, dbar = tp$dbar,
      psiT1 = psv(tb, "R", 1) + psv(tb, "F", 1),
      psiR1 = psv(tb, "R", 1), psiF1 = psv(tb, "F", 1),
      share_rec = sd$shares[["reciprocal"]],
      G_rec = gam(tb, "F", 1), G_c3 = gam(tb, "F", 2), G_c4 = gam(tb, "F", 3),
      Gw_rec = gam(sw$table, "F", 1), Gw_c3 = gam(sw$table, "F", 2),
      Gw_c4 = gam(sw$table, "F", 3),
      Gn_rec = gam(sn$table, "F", 1), Gn_c3 = gam(sn$table, "F", 2),
      n_in_CHN_neutral = if ("CHN" %in% dc$country)
        sum(An[match("CHN", dc$country), ]) else NA),
    M = list(R1 = sd$mats$R1, F1 = sd$mats$F1, S = S, A = A,
             within = within_order_shares(sd)))
}

# Propagation on the country-industry network of foreign intermediate inputs.
# Column b of P = shares of b's foreign inputs supplied by each foreign
# country-industry (upstream direction).  Nodes: positive VA and positive
# foreign inputs; columns are renormalized on this set.
# NOTE: the paper does not say whether shares are renormalized after dropping
# nodes; check against the original code.
ci_propagation <- function(d, kmax = KMAX) {
  same <- outer(d$country, d$country, "==")
  Zf <- d$Z; Zf[same] <- 0
  keep <- d$va > 0 & colSums(Zf) > 0
  repeat {                                   # renormalized set must be closed
    k2 <- keep & colSums(Zf[keep, , drop = FALSE]) > 0
    if (all(k2 == keep)) break; keep <- k2
  }
  Zk <- Zf[keep, keep]
  P  <- sweep(Zk, 2, colSums(Zk), "/")
  w  <- d$va[keep] / sum(d$va[keep])
  cc <- d$country[keep]; lab <- paste(cc, d$industry[keep], sep = "_")
  Pk <- diag(length(w)); out <- list()
  for (k in seq_len(kmax)) {
    Pk <- if (k == 1) P else Pk %*% P      # the only O(N^3) step
    wk <- as.numeric(Pk %*% w)
    Y  <- colSums(Pk^2)
    own <- sum(w^2 * Y)
    # Proposition 3(i) gives psi_O without forming O = Pk'Pk (O(N^3)):
    #   (1 - HHI) psi_O = HHI(w^(k)) - sum_j w_j^2 Y_j
    # and sum of off-diagonal O = ||Pk 1||^2 - sum_j Y_j gives mbar_O.
    N <- length(w); psiO <- (sum(wk^2) - own) / (1 - sum(w^2))
    mbarO <- (sum(rowSums(Pk)^2) - sum(Y)) / (N * (N - 1))
    out[[k]] <- data.frame(k = k, n_nodes = length(w), psiO = psiO,
      GammaO = psiO / mbarO, hhi_wk = sum(wk^2),
      share_conv = 1 - own / sum(wk^2),
      share_CHN = sum(wk[cc == "CHN"]), share_USA = sum(wk[cc == "USA"]),
      top_dest = lab[which.max(wk)], top_share = max(wk))
  }
  do.call(rbind, out)
}

# ---- run over years ------------------------------------------------------------
# NCI_FIG_ONLY=1 redraws tables and figures from the saved results
RES_FILE <- file.path(DIR_PROC, paste0("results_", SRC, ".rds"))
FIG_ONLY <- Sys.getenv("NCI_FIG_ONLY") == "1" && file.exists(RES_FILE)
rows <- list(); mats <- list(); prop <- list(); robust <- list()
for (yr in if (FIG_ONLY) integer(0) else YEARS) {
  d  <- load_wiod_year(yr)
  dc <- wiod_country(d)
  w  <- dc$VA / sum(dc$VA)
  cm <- country_measures(dc, w, THETA)
  rows[[as.character(yr)]] <- cbind(year = yr, cm$row)
  mats[[as.character(yr)]] <- c(cm$M, list(w = w, country = dc$country))
  prop[[as.character(yr)]] <- cbind(year = yr,
    ci_propagation(d, kmax = if (yr %in% TABLE_YEARS) KMAX else 1))
  for (th in THETAS) robust[[paste(yr, th)]] <-
    cbind(year = yr, theta = th, country_measures(dc, w, th)$row)
  log_msg("WIOD ", yr, " done")
}
if (FIG_ONLY) {
  R0 <- readRDS(RES_FILE); ts <- R0$ts; prop <- R0$prop; rob <- R0$rob; mats <- R0$mats
} else {
  ts   <- do.call(rbind, rows)
  prop <- do.call(rbind, prop)
  rob  <- do.call(rbind, robust)
  saveRDS(list(ts = ts, prop = prop, rob = rob, mats = mats), RES_FILE)
}

# ---- Table 2 ---------------------------------------------------------------------
ty <- ts[ts$year %in% TABLE_YEARS, ]
py <- prop[prop$year %in% TABLE_YEARS & prop$k == 1, ]
f3 <- function(x) fmt_num(x, 3); f2 <- function(x) fmt_num(x, 2)
pair <- function(a, b, f = f2) paste(f(a), f(b), sep = "; ")
mk <- function(dim, lab, vals) setNames(data.frame(dim, lab, t(vals)),
                                        c("Dimension", "Measure", TABLE_YEARS))
tab2 <- rbind(
  mk("Marginal", "Herfindahl index of $w$", f3(ty$hhi)),
  mk("", "Density of the network", f3(ty$density)),
  mk("Direct", "$\\psi_D(w)$; $\\Gamma_1(w)$", pair(ty$psiD, ty$Gamma1, f3)),
  mk("Topological", "$\\psi^{(2)}_G(w)$; $\\Gamma_2(w)$", pair(ty$psiG2, ty$Gamma2, f3)),
  mk("", "Economic mean distance $\\bar d_w$; average path length $\\bar d$",
     pair(ty$dw, ty$dbar)),
  mk("Structural", "$\\psi^{(1)}_T(w)$: one-way + reciprocal", f3(ty$psiT1)),
  mk("", "Share of reciprocal links in $\\psi_{TOT}$", f2(ty$share_rec)),
  mk("", "$\\Gamma$ of reciprocal links", f2(ty$G_rec)),
  mk("", "$\\Gamma$ of circuits of length 3; of length 4", pair(ty$G_c3, ty$G_c4)),
  mk("", "Weighted: $\\Gamma$ of reciprocal links", f2(ty$Gw_rec)),
  mk("", "Weighted: $\\Gamma$ of circuits of length 3; of length 4",
     pair(ty$Gw_c3, ty$Gw_c4)),
  mk("Propagation", "$\\psi^{(1)}_O(w)$; $\\Gamma_{O^{(1)}}(w)$",
     paste(f3(py$psiO), f2(py$GammaO), sep = "; ")),
  mk("", "Share of convergence in the Herfindahl index of $w^{(1)}$",
     f2(py$share_conv)),
  mk("", "Share of $w^{(1)}$ landing in China; in the United States",
     pair(py$share_CHN, py$share_USA)))
write_tex_table(tab2, paste0("tab2_", SRC, ".tex"),
  caption = paste0("Concentration of the world production network, selected years",
                   if (SRC == "icio") " (OECD ICIO)" else ""),
  label = paste0("tab:", SRC), align = ALIGN("ll"),
  header = c("Dimension", "Measure", TABLE_YEARS),
  notes = paste0(if (SRC == "icio") "OECD Inter-Country Input--Output tables, 2025 edition (China and Mexico aggregated, rest of the world excluded); "
                 else "World Input--Output Database, 2016 release; ", nrow(mats[[1]]$R1),
    " countries weighted by value added. Country network with threshold $\\theta=",
    100 * THETA, "\\%$; weighted rows use the input shares $Z_{ji}/X_i$ as link",
    " intensities, without threshold. Structural indices with $H=3$ and",
    " $\\alpha_h=1/3$; circuits of length 3 and 4 correspond to $h=2$ and $h=3$.",
    " Propagation measures on the network of foreign intermediate inputs among ",
    py$n_nodes[1], " country--industry pairs, one step upstream."))

# ---- Online Appendix OA.1: robustness of the reciprocal/circuit alignment --------
oa <- rob[rob$year %in% TABLE_YEARS, ]
oa_tab <- do.call(rbind, lapply(THETAS, function(th) {
  o <- oa[oa$theta == th, ]
  data.frame(Spec = sprintf("Binary, $\\theta=%s\\%%$", 100 * th),
             Measure = c("$\\Gamma$ reciprocal", "$\\Gamma$ circuits (3)",
                         "Share of reciprocal in $\\psi_{TOT}$"),
             rbind(f2(o$G_rec), f2(o$G_c3), f2(o$share_rec)))
}))
o <- oa[oa$theta == THETA, ]
oa_tab <- rbind(oa_tab,
  data.frame(Spec = "Weighted, $Z_{ji}/X_i$", Measure = c("$\\Gamma$ reciprocal",
             "$\\Gamma$ circuits (3)", "$\\Gamma$ circuits (4)"),
             rbind(f2(o$Gw_rec), f2(o$Gw_c3), f2(o$Gw_c4))),
  data.frame(Spec = "Size-neutral binary", Measure = c("$\\Gamma$ reciprocal",
             "$\\Gamma$ circuits (3)", "Incoming links of China"),
             rbind(f2(o$Gn_rec), f2(o$Gn_c3), as.character(o$n_in_CHN_neutral))))
names(oa_tab)[-(1:2)] <- TABLE_YEARS
write_tex_table(oa_tab, paste0("tabOA1_", SRC, "_robustness.tex"),
  caption = "World production network: threshold, weighted and size-neutral specifications",
  label = paste0("tab:oa-", SRC), align = ALIGN("ll"), header = c("Specification", "Measure",
  TABLE_YEARS))

# ---- κ-free within-order shares (fix for the cross-order share problem) ------
wos <- do.call(rbind, lapply(TABLE_YEARS, function(y)
  cbind(year = y, mats[[as.character(y)]]$within)))
tab_wos <- data.frame(Order = sprintf("$h=%d$ (%s)", wos$h,
                        c("reciprocal vs one-way", "circuits of length 3 vs open chains",
                          "circuits of length 4 vs open chains")[wos$h]),
                      Year = wos$year, Closed = fmt_num(wos$closed_share, 3),
                      Gopen = fmt_num(wos$Gamma_open, 2), Gclosed = fmt_num(wos$Gamma_closed, 2))
write_tex_table(tab_wos[order(wos$h, wos$year), ], paste0("tabOA_", SRC, "_within_order.tex"),
  caption = "World production network: within-order structural decomposition",
  label = paste0("tab:", SRC, "-within"), align = "lcccc",
  header = c("Order", "Year", "Closed share of $\\psi^{(h)}_T$", "$\\Gamma$ open", "$\\Gamma$ closed"),
  notes = paste("Within an order, open and closed components share the normalising",
    "constant $\\kappa_{N,h}$, so the closed share does not depend on it."))

# ---- permutation band for Gamma of reciprocal links ---------------------------
set.seed(11)
band <- do.call(rbind, lapply(names(mats), function(y) {
  pa <- perm_alignment(mats[[y]]$w, mats[[y]]$F1, B = 500)
  data.frame(x = as.numeric(y), lo = quantile(pa$null, .05), hi = quantile(pa$null, .95),
             p = pa$p_plus1)
}))

# ---- Figure 2 --------------------------------------------------------------------
ts$x <- ts$year
pa <- fig_levels(ts, "year", c("hhi", "psiT1", "psiR1", "psiF1"),
  c("HHI of w", "Direct", "One-way", "Reciprocal"),
  c(COL_HHI, INK, COL_MECH[["One-way"]], COL_MECH[["Reciprocal"]]),
  c("22", "solid", "solid", "solid"),
  title = "(a) Size concentration and first-order structure",
  subtitle = expression("Direct "*psi[T]^{(1)}*" = one-way "*psi[R]^{(1)}*" + reciprocal "*psi[F]^{(1)}),
  nudge = .4)
pb <- fig_alignment(ts, "year", transform(band, lo = pmax(lo, .2)),
  "(b) Alignment with economic size", log_y = SRC == "icio",
  "Shaded: 5–95% of reciprocal Γ under random weights", nudge = .4)
decomp_periods <- Filter(function(p) all(p %in% YEARS),
  if (SRC == "icio") list(c(1995, 2007), c(2007, 2014), c(2008, 2009), c(2014, 2019),
                          c(2019, 2020), c(2020, 2022))
  else list(c(2000, 2007), c(2007, 2014), c(2008, 2009)))
dec <- do.call(rbind, lapply(decomp_periods, function(p) {
  a <- mats[[as.character(p[1])]]; b <- mats[[as.character(p[2])]]
  lab <- sprintf("%d-%02d", p[1], p[2] %% 100)
  rbind(cbind(comp = "OW", period = lab, dynamic_decomposition(a$w, a$R1, b$w, b$R1)),
        cbind(comp = "RE", period = lab, dynamic_decomposition(a$w, a$F1, b$w, b$F1)))
}))
write.csv(dec, file.path(DIR_TAB, paste0("fig2c_decomposition_", SRC, ".csv")), row.names = FALSE)
pc <- fig_decomp(dec, "(c) Sources of change (Proposition 5)")
D_YEARS <- TABLE_YEARS[unique(round(seq(1, length(TABLE_YEARS), length.out = 3)))]
pdd <- prop[prop$year %in% D_YEARS, ]
ld <- rbind(data.frame(k = pdd$k, v = pdd$share_CHN, cty = "China", year = pdd$year),
            data.frame(k = pdd$k, v = pdd$share_USA, cty = "United States", year = pdd$year))
ld$year <- factor(ld$year)
COL_CTY2 <- c("China" = CAT[8], "United States" = CAT[1])
lab_d <- ld[ld$k == KMAX, ]; lab_d$lab <- paste(lab_d$cty, lab_d$year)
pd <- ggplot(ld, aes(k, v, colour = cty, group = interaction(cty, year), alpha = year)) +
  geom_line(linewidth = .7) + geom_point(size = 2) +
  geom_text_repel(data = lab_d, aes(label = lab), nudge_x = .25, hjust = 0,
                  direction = "y", size = 2.3, segment.colour = GRID, alpha = 1,
                  show.legend = FALSE, seed = 1) +
  scale_colour_manual(values = COL_CTY2) +
  scale_alpha_manual(values = setNames(seq(.3, 1, length.out = nlevels(ld$year)),
                                       levels(ld$year))) +
  scale_x_continuous(breaks = 1:KMAX, expand = expansion(mult = c(.05, .45))) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Upstream step k", y = NULL, title = "(d) Where upstream dependence lands",
       subtitle = "Share of value-added-weighted foreign dependence") +
  theme_paper() + theme(legend.position = "none")
bottom <- if (length(decomp_periods) > 3) (pc | pd) + plot_layout(widths = c(1.7, 1)) else (pc | pd)
save_fig((pa | pb) / bottom, paste0("fig2_", SRC), width = if (length(decomp_periods) > 3) 8.5 else 7.5,
         height = 7.4)

# ---- Figure 5a: who drives the alignment of reciprocal links ------------------
PAIR_YEARS <- TABLE_YEARS[unique(round(seq(1, length(TABLE_YEARS), length.out = 3)))]
pcs <- do.call(rbind, lapply(PAIR_YEARS, function(y) {
  m <- mats[[as.character(y)]]
  p <- pair_contributions(m$w, m$F1, labels = m$country, top = 8)
  G <- attr(p, "Gamma")
  cbind(p, date = sprintf("%d (Γ = %.2f)", y, G), share = p$contribution / G)
}))
write.csv(pcs, file.path(DIR_TAB, paste0("fig5_pairs_", SRC, ".csv")), row.names = FALSE)
p5 <- fig_pairs(pcs, "Pairs driving the alignment of reciprocal links",
  "Top-8 country pairs; labels: share of Γ; Γ is the sum of the contributions of all pairs")
save_fig(p5, paste0("fig5a_pairs_", SRC), width = 7.5, height = 3.2)

# ---- Figure OA: threshold profile of Gamma (binary) vs weighted ----------------
THG <- c(0.002, 0.005, 0.01, 0.02)
tpl <- do.call(rbind, lapply(names(mats), function(y) {
  m <- mats[[y]]
  gb <- gamma_threshold_profile(m$w, m$S, THG)
  rbind(data.frame(x = as.numeric(y), v = gb, s = sprintf("Binary, %.1f%%", 100 * THG)),
        data.frame(x = as.numeric(y), v = ts$Gw_rec[ts$year == as.numeric(y)],
                   s = "Weighted"))
}))
# thresholds with no reciprocal links in any year (2% in the production network) are dropped
tpl <- tpl[!is.na(tpl$v), ]
n_th <- length(unique(tpl$s)) - 1
pt <- fig_threshold(tpl, "x", "Threshold dependence of the reciprocal alignment",
  sprintf("Binary networks at %s thresholds (blue, light to dark) and weighted network (orange)",
          c("one", "two", "three", "four")[n_th]),
  nudge = .4)
save_fig(pt, paste0("figOA_threshold_", SRC), width = 6, height = 3.6)
log_msg("Table 2 and Figures 2, 5a, OA written")
