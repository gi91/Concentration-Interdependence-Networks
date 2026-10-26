# =============================================================================
# 09_sim_design.R -- simulation design on network size and topology
#
# Questions:
#  Q1 (inference)   Do the tests keep their size when node importance is
#                   correlated with connectivity, as in real networks?  The
#                   permutation null breaks the weight-degree link, the
#                   degree-conditional null (DBCM; exact rewiring at N = 19)
#                   keeps it.  Which one answers "is reciprocity aligned with
#                   size beyond what degrees imply"?
#  Q2 (power)       How much power do the tests have at the N and the weight
#                   concentration of the applications?  (Needed to read the
#                   2014 non-rejection in the production network.)
#  Q3 (comparability) Which quantities are comparable across N and topology:
#                   Gamma and the within-order closed share (kappa-free) or the
#                   cross-order shares of psi_TOT?
# Factorial design: topology (4) x N (19, 43, 100) x rho (0, 0.8) x delta (0, 0.5)
# Output: figS4, figS5, figS6, tabS3 (+ csv)
# =============================================================================
source("R/utils.R"); source("R/core_indices.R"); source("R/null_models.R")
source("R/fixes.R")
REPS <- as.integer(Sys.getenv("NCI_REPS_DESIGN", 200))
BT   <- as.integer(Sys.getenv("NCI_B_DESIGN", 199))
KBAR <- 6                                   # expected out-degree
set.seed(2026)

# ---- generators: directed binary A[to, from] with E(out-degree) ~ KBAR -------
# f: latent node fitness (drives degree where the topology allows it)
gen_network <- function(topo, N, f) {
  p0 <- KBAR / (N - 1)
  P <- switch(topo,
    "Erdos-Renyi" = matrix(p0, N, N),
    "Chung-Lu" = { P <- outer(f, f); P * (p0 * N * (N - 1) / sum(P[row(P) != col(P)])) },
    "Core-periphery" = { core <- f >= quantile(f, .8)
      P <- outer(core, core, function(a, b) ifelse(a & b, 8, ifelse(a | b, 3, 1)))
      P * (p0 * N * (N - 1) / sum(P[row(P) != col(P)])) },
    "Gravity" = { Z <- outer(f^.6, f^.6) * exp(matrix(rnorm(N * N, 0, 1), N, N))
      diag(Z) <- 0; S <- Z / rowSums(Z)
      thr <- threshold_for_density(S, p0); return((S >= thr) * 1) })
  P <- pmin(P, 1); diag(P) <- 0
  A <- (matrix(runif(N * N), N, N) < P) * 1; diag(A) <- 0; A
}

draw_case <- function(topo, N, rho, delta) {
  zf <- rnorm(N); f <- exp(1.0 * zf)                     # heavy-tailed fitness
  zw <- rho * zf + sqrt(1 - rho^2) * rnorm(N)
  wv <- sort(rlnorm(N, 0, 1.2)); w <- wv[rank(zw)]; w <- w / sum(w)
  A <- gen_network(topo, N, f)
  # planted alignment: reciprocal closure among the 4 heaviest nodes
  top <- order(w, decreasing = TRUE)[1:4]
  for (i in top) for (j in top) if (i != j && A[i, j] == 1 && runif(1) < delta) A[j, i] <- 1
  list(A = A, w = w)
}

one_rep <- function(topo, N, rho, delta) {
  cs <- draw_case(topo, N, rho, delta); A <- cs$A; w <- cs$w
  F <- A * t(A)
  if (sum(F) == 0) return(NULL)
  pp <- perm_alignment(w, F, B = BT)
  psiF <- function(X) psi_index(w, X * t(X))
  fit <- dbcm_fit(A)
  db <- dbcm_test(A, psiF, B = BT, fit = fit)
  sd <- structural_decomposition(w, A, H = 3)
  wo <- within_order_shares(sd)
  # degree-conditional benchmark of the closed share (invariance cells only)
  cs_db <- if (delta == 0 && rho == 0) mean(replicate(10, {
    Xs <- dbcm_sample(fit)
    within_order_shares(structural_decomposition(w, Xs, H = 2, alpha = c(.5, .5)))$closed_share[2]
  }), na.rm = TRUE) else NA
  data.frame(topo = topo, N = N, rho = rho, delta = delta, hhi = hhi(w),
             cor_w_deg = cor(w, colSums(A) + rowSums(A), method = "spearman"),
             GammaF = pp$Gamma, p_perm = pp$p_plus1, p_dbcm = db$p,
             closed_share_h2 = wo$closed_share[2],
             closure_ratio_h2 = wo$closed_share[2] / (sum(A) / (N * (N - 1))),
             closure_ratio_dbcm = wo$closed_share[2] / cs_db,
             circuit_share_tot = sd$shares[["circuits"]],
             recip_share_tot = sd$shares[["reciprocal"]])
}

TOPOS <- c("Erdos-Renyi", "Chung-Lu", "Core-periphery", "Gravity")
grid <- expand.grid(topo = TOPOS, N = c(19, 43, 100), rho = c(0, 0.8),
                    delta = c(0, 0.5), stringsAsFactors = FALSE)
# NCI_SIM_REUSE=1 redraws figures and tables from the saved simulation
REUSE <- Sys.getenv("NCI_SIM_REUSE") == "1" &&
  file.exists(file.path(DIR_PROC, "results_sim_design.rds"))
if (REUSE) {
  res <- readRDS(file.path(DIR_PROC, "results_sim_design.rds"))
} else {
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(g) {
  r <- do.call(rbind, lapply(seq_len(REPS), function(k)
    one_rep(grid$topo[g], grid$N[g], grid$rho[g], grid$delta[g])))
  if (g %% 8 == 0) log_msg("design cell ", g, "/", nrow(grid))
  r
}))
saveRDS(res, file.path(DIR_PROC, "results_sim_design.rds"))
}

agg <- aggregate(cbind(rej_perm = p_perm <= .05, rej_dbcm = p_dbcm <= .05, GammaF,
                       closed_share_h2, closure_ratio_h2, closure_ratio_dbcm,
                       circuit_share_tot,
                       recip_share_tot, hhi,
                       cor_w_deg) ~ topo + N + rho + delta,
                 data = res, FUN = function(x) mean(x, na.rm = TRUE), na.action = na.pass)
write.csv(agg, file.path(DIR_TAB, "simS3_design.csv"), row.names = FALSE)

# ---- exact rewiring vs DBCM at N = 19 (size, rho = 0.8) -------------------------
if (!REUSE) {
chk <- do.call(rbind, lapply(TOPOS, function(topo) {   # Section 3.2: exact rewiring vs DBCM
  do.call(rbind, lapply(seq_len(min(REPS, 100)), function(k) {
    cs <- draw_case(topo, 19, 0.8, 0); A <- cs$A
    if (sum(A * t(A)) == 0 || any(colSums(A) + rowSums(A) == 0)) return(NULL)
    st <- function(X) c(psi_F1 = psi_index(cs$w, X * t(X)))
    rw <- rewiring_test(A, st, B = 99)
    db <- dbcm_test(A, function(X) st(X)[[1]], B = 99)
    data.frame(topo = topo, p_rewire = rw$p_plus1[[1]], p_dbcm = db$p)
  }))
}))
chk_s <- aggregate(cbind(rej_rewire = p_rewire <= .05, rej_dbcm = p_dbcm <= .05,
                         cor_p = p_rewire) ~ topo, data = chk, FUN = mean)
chk_s$cor_p <- sapply(chk_s$topo, function(t)
  cor(chk$p_rewire[chk$topo == t], chk$p_dbcm[chk$topo == t]))
write.csv(chk_s, file.path(DIR_TAB, "simS3_rewire_vs_dbcm.csv"), row.names = FALSE)
log_msg("rewiring vs DBCM check done")
}

# ---- figures --------------------------------------------------------------------
COL_TEST <- c("Permutation of weights" = INK, "Degree-conditional (DBCM)" = CAT[1])
long_rej <- function(a) rbind(
  data.frame(a[, c("topo", "N", "rho", "delta")], rej = a$rej_perm,
             test = "Permutation of weights"),
  data.frame(a[, c("topo", "N", "rho", "delta")], rej = a$rej_dbcm,
             test = "Degree-conditional (DBCM)"))
plot_rej <- function(L, title, subtitle, band) {
  L$rho_l <- factor(ifelse(L$rho == 0, "weights unrelated to degree",
                           "weights correlated with degree (ρ = 0.8)"),
                    c("weights unrelated to degree", "weights correlated with degree (ρ = 0.8)"))
  L$topo <- factor(L$topo, TOPOS)
  ggplot(L, aes(N, rej, colour = test, linetype = rho_l, shape = test)) +
    (if (band) annotate("rect", xmin = 17, xmax = 112, ymin = .025, ymax = .075,
                        fill = SHADE) else NULL) +
    geom_hline(yintercept = .05, colour = MUTED, linewidth = .3, linetype = "22") +
    geom_line(linewidth = .6) + geom_point(size = 1.8) +
    facet_wrap(~ topo, nrow = 1) +
    scale_colour_manual(values = COL_TEST) +
    scale_linetype_manual(values = c("solid", "42")) +
    scale_x_continuous(breaks = c(19, 43, 100), trans = "log10") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = "Number of nodes N (log scale)", y = "Rejection rate at 5%",
         title = title, subtitle = subtitle) +
    guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2),
           shape = guide_legend(nrow = 2)) + theme_paper()
}
save_fig(plot_rej(long_rej(agg[agg$delta == 0, ]),
  "Size: no alignment of reciprocity beyond what degrees imply",
  "With size correlated with degree (dashed) the permutation test rejects; the DBCM test holds, except with heavy-tailed degrees",
  band = TRUE), "figS4_size_topology", width = 7.5, height = 3.8)
save_fig(plot_rej(long_rej(agg[agg$delta > 0, ]),
  "Power: reciprocal closure among the four heaviest nodes (delta = 0.5)",
  "Power is moderate at N = 19-43; with rho = 0.8 the permutation test mixes the effect with degrees",
  band = FALSE), "figS5_power_topology", width = 7.5, height = 3.8)

inv <- agg[agg$delta == 0 & agg$rho == 0, ]
L6 <- rbind(
  data.frame(inv[, c("topo", "N")], v = inv$GammaF, q = "Γ of reciprocal links"),
  data.frame(inv[, c("topo", "N")], v = inv$closed_share_h2,
             q = "Within-order closed share, h = 2 (κ-free)"),
  data.frame(inv[, c("topo", "N")], v = inv$closure_ratio_h2,
             q = "Closed share / directed density"),
  data.frame(inv[, c("topo", "N")], v = inv$closure_ratio_dbcm,
             q = "Closed share / degree-conditional expectation"),
  data.frame(inv[, c("topo", "N")], v = inv$circuit_share_tot,
             q = "Cross-order share of circuits in ψ_TOT"))
L6$q <- factor(L6$q, unique(L6$q)); L6$topo <- factor(L6$topo, TOPOS)
COL_TOPO <- setNames(c(INK, CAT[8], CAT[6], CAT[7]), TOPOS)
p6 <- ggplot(L6, aes(N, v, colour = topo, shape = topo)) +
  geom_line(linewidth = .6) + geom_point(size = 1.8) +
  facet_wrap(~ q, scales = "free_y", nrow = 1, labeller = label_wrap_gen(24)) +
  scale_colour_manual(values = COL_TOPO) +
  scale_x_continuous(breaks = c(19, 43, 100), trans = "log10") +
  expand_limits(y = 0) +
  labs(x = "Number of nodes N (log scale)", y = NULL,
       title = "Comparability across network size and topology (no planted structure)",
       subtitle = "Comparable across N and topology: Γ and the closed share relative to its degree-conditional expectation") +
  guides(colour = guide_legend(nrow = 1)) + theme_paper()
save_fig(p6, "figS6_comparability", width = 9, height = 3.8)

# ---- summary table --------------------------------------------------------------
pf <- function(x) formatC(100 * x, format = "f", digits = 0)
tab <- agg[order(agg$delta, agg$rho, match(agg$topo, TOPOS), agg$N), ]
tabS3 <- data.frame(Topology = tab$topo, N = tab$N, rho = tab$rho,
  HHI = fmt_num(tab$hhi, 3), cor = fmt_num(tab$cor_w_deg, 2),
  perm = pf(tab$rej_perm), dbcm = pf(tab$rej_dbcm), G = fmt_num(tab$GammaF, 2),
  ws = fmt_num(tab$closed_share_h2, 3), wr = fmt_num(tab$closure_ratio_h2, 2),
  cs = fmt_num(tab$circuit_share_tot, 4))
write_tex_table(tabS3, "tabS3_design.tex",
  panel = ifelse(tab$delta == 0, "Size: no planted alignment ($\\delta=0$)",
                 "Power: reciprocal closure among the 4 heaviest nodes ($\\delta=0.5$)"),
  caption = "Simulation design: network size, topology and weight--degree correlation",
  label = "tab:sim-design", align = "lcccccccccc",
  header = c("Topology", "$N$", "$\\rho$", "HHI", "Corr($w$, deg)",
             "Perm.\\ (\\%)", "DBCM (\\%)", "$\\Gamma_F$", "Closed share $h=2$",
             "Closed share / density", "Circuit share $\\psi_{TOT}$"),
  notes = paste0(REPS, " replications per cell, ", BT, " null draws per test,",
    " expected out-degree ", KBAR, ". Rejection rates at the 5\\% level."))
log_msg("design simulation done")
