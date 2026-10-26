# =============================================================================
# 08_sim_robustness.R -- Monte Carlo experiments on the assumptions behind
# the empirical results.  Each block writes a figure and/or table (prefix S).
#   B1 power of the permutation test vs. concentration of the weights
#   B2 threshold dependence of binary alignment factors vs. weighted version
#   B3 dependence of the psi_TOT shares on the kappa normalisation and on N
#   B4 Proposition 5 (nested) vs. two-player Shapley decomposition
#   B5 small-sample inference in the SRISK regressions: HC1/HC3/Freedman-Lane
#      size, and spurious significance of reciprocal embeddedness when only
#      importance as a counterparty matters
#   B6 contagion experiment on synthetic topologies (ER, core-periphery)
# =============================================================================
source("R/utils.R"); source("R/core_indices.R"); source("R/null_models.R")
source("R/validation.R"); source("R/contagion.R")
REPS <- as.integer(Sys.getenv("NCI_REPS", 300))
set.seed(99)

lognorm_w <- function(N, s) { x <- rlnorm(N, 0, s); x / sum(x) }
gamF <- function(w, A) { F <- A * t(A); alignment(w, F) }

# ---- B1: power of the permutation test -------------------------------------------
# DGP: directed ER; among the k heaviest nodes each pair is made reciprocal with
# probability delta (delta = 0: null).  Weights log-normal with dispersion s.
B1 <- expand.grid(N = c(19, 43), s = c(0.5, 1, 1.5, 2), delta = c(0, 0.3, 0.6))
B1$rej_perm <- B1$rej_mantel <- B1$hhi <- NA
for (r in seq_len(nrow(B1))) {
  N <- B1$N[r]; p <- if (N == 19) 0.3 else 0.12; k <- 4
  out <- replicate(REPS, {
    w <- lognorm_w(N, B1$s[r]); A <- (matrix(runif(N^2), N, N) < p) * 1
    diag(A) <- 0
    top <- order(w, decreasing = TRUE)[1:k]
    for (i in top) for (j in top) if (i < j && runif(1) < B1$delta[r]) {
      A[i, j] <- 1; A[j, i] <- 1 }
    F <- A * t(A)
    if (sum(F) == 0) return(c(NA, NA, hhi(w)))
    pa <- perm_alignment(w, F, B = 199); mm <- mantel_moments(w, F)
    c(pa$p_plus1 <= 0.05, mm$z > qnorm(0.95), hhi(w)) })
  B1$rej_perm[r] <- mean(out[1, ], na.rm = TRUE)
  B1$rej_mantel[r] <- mean(out[2, ], na.rm = TRUE)
  B1$hhi[r] <- mean(out[3, ])
}
write.csv(B1, file.path(DIR_TAB, "simB1_power.csv"), row.names = FALSE)
l1 <- rbind(data.frame(B1[, c("N", "s", "delta", "hhi")], rej = B1$rej_perm,
                       test = "Permutation"),
            data.frame(B1[, c("N", "s", "delta", "hhi")], rej = B1$rej_mantel,
                       test = "Mantel normal approx."))
l1$panel <- paste0("N = ", l1$N); l1$dl <- paste0("delta = ", l1$delta)
l1$dl <- factor(ifelse(l1$delta == 0, "No planted effect (size)",
                       paste0("Planted reciprocity, delta = ", l1$delta)))
COL_D <- setNames(c(MUTED, SEQ3[1], SEQ3[3]), levels(l1$dl))
pB1 <- ggplot(l1, aes(hhi, rej, colour = dl, linetype = test, shape = test,
                      group = interaction(dl, test))) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.025, ymax = 0.075,
           fill = SHADE) +
  geom_hline(yintercept = 0.05, colour = MUTED, linewidth = .3, linetype = "22") +
  geom_line(linewidth = .6) + geom_point(size = 1.8) +
  facet_wrap(~ panel, scales = "free_x") +
  scale_colour_manual(values = COL_D) +
  scale_linetype_manual(values = c("Permutation" = "solid", "Mantel normal approx." = "22")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Herfindahl index of the weights", y = "Rejection rate at 5%",
       title = "Size and power of the alignment test for reciprocal links",
       subtitle = "Grey band: acceptable size range; planted effect among the 4 heaviest nodes") +
  guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2),
         shape = guide_legend(nrow = 2)) + theme_paper()
save_fig(pB1, "figS1_power_permutation", 7, 4)
log_msg("B1 done")

# ---- B2: threshold dependence ------------------------------------------------------
# Gravity flows with a symmetric pair component (reciprocity) of strength rho.
gravity <- function(N, rho, s) {
  size <- rlnorm(N, 0, s); u <- matrix(rnorm(N^2), N, N); v <- matrix(rnorm(N^2), N, N)
  v <- (v + t(v)) / sqrt(2)
  Z <- outer(size^0.6, size^0.6) * exp(sqrt(1 - rho) * u + sqrt(rho) * v)
  diag(Z) <- 0; X <- rowSums(Z) + colSums(Z) + size * 5
  list(S = t(Z) / X, w = size / sum(size))       # S[i, j] = Z_ji / X_i
}
thetas <- c(0.002, 0.005, 0.01, 0.02)
B2 <- do.call(rbind, lapply(1:min(REPS, 200), function(r) {
  g <- gravity(43, rho = 0.5, s = 1.5)
  gb <- sapply(thetas, function(th) { A <- (g$S >= th) * 1; gamF(g$w, A) })
  gw <- gamF(g$w, pmin(g$S / max(g$S), 1))
  data.frame(rep = r, theta = c(thetas, 0), spec = c(rep("binary", 4), "weighted"),
             Gamma = c(gb, gw))
}))
B2s <- aggregate(Gamma ~ spec + theta, data = B2,
                 FUN = function(x) c(mean = mean(x), sd = sd(x),
                                     cv = sd(x) / mean(x)))
B2s <- do.call(data.frame, B2s)
write.csv(B2s, file.path(DIR_TAB, "simB2_threshold.csv"), row.names = FALSE)
# rank correlation across thresholds within replications: does the ordering of
# networks by Gamma survive a change in theta?
wide <- reshape(B2[, c("rep", "theta", "Gamma")], idvar = "rep", timevar = "theta",
                direction = "wide")
cors <- cor(wide[, -1], method = "spearman", use = "pairwise.complete.obs")
write.csv(round(cors, 3), file.path(DIR_TAB, "simB2_rank_correlations.csv"))
log_msg("B2 done")

# ---- B3: shares of psi_TOT depend on kappa and N ------------------------------
B3 <- do.call(rbind, lapply(c(8, 19, 43, 80), function(N) {
  do.call(rbind, lapply(1:min(REPS, 50), function(r) {
    w <- lognorm_w(N, 1); A <- (matrix(runif(N^2), N, N) < 4 / N) * 1; diag(A) <- 0
    sd <- structural_decomposition(w, A, H = 3)
    tb <- sd$table
    data.frame(N = N, t(sd$shares),
               G_open2 = tb$Gamma[tb$comp == "R" & tb$h == 2],
               G_circ2 = tb$Gamma[tb$comp == "F" & tb$h == 2])
  }))
}))
B3s <- aggregate(. ~ N, data = B3, FUN = function(x) mean(x, na.rm = TRUE),
                 na.action = na.pass)
write.csv(B3s, file.path(DIR_TAB, "simB3_shares_vs_N.csv"), row.names = FALSE)
l3 <- reshape(B3s[, c("N", "one_way", "reciprocal", "open_chains", "circuits")],
              direction = "long", varying = list(2:5), v.names = "share",
              timevar = "mech", times = c("One-way", "Reciprocal", "Open chains",
                                          "Circuits"))
l3$mech <- factor(l3$mech, names(COL_MECH))
pB3 <- ggplot(l3, aes(N, share, colour = mech)) + geom_line(linewidth = .7) +
  geom_point(size = 2) + end_labels(l3, "N", "share", "mech", COL_MECH, nudge = 3) +
  scale_colour_manual(values = COL_MECH) +
  scale_x_continuous(breaks = c(8, 19, 43, 80), expand = expansion(mult = c(.03, .3))) +
  scale_y_continuous(labels = scales::percent) + coord_cartesian(clip = "off") +
  labs(x = "Number of nodes N (mean degree held at 4, no planted structure)",
       y = expression("Share of "*psi[TOT]),
       title = "Cross-order shares of the total structural index depend on N",
       subtitle = "Higher orders shrink with the normalisation by the complete network") +
  theme_paper() + theme(legend.position = "none")
save_fig(pB3, "figS3_shares_vs_N", 6, 3.8)
log_msg("B3 done")

# ---- B4: nested (Prop. 5) vs two-player Shapley decomposition ----------------------
B4 <- do.call(rbind, lapply(seq_len(REPS), function(r) {
  N <- sample(c(19, 43), 1); p0 <- runif(1, .1, .4)
  A0 <- (matrix(runif(N^2), N, N) < p0) * 1; diag(A0) <- 0
  A1 <- (matrix(runif(N^2), N, N) < runif(1, .1, .4)) * 1; diag(A1) <- 0
  w0 <- lognorm_w(N, 1); w1 <- w0 * exp(rnorm(N, 0, .5)); w1 <- w1 / sum(w1)
  d <- dynamic_decomposition(w0, A0, w1, A1)
  data.frame(delta = d$delta, nested_w = d$weights,
             shapley_w = d$shapley_weights,
             nested_net = d$availability + d$network, shapley_net = d$shapley_network)
}))
B4$gap_rel <- abs(B4$nested_w - B4$shapley_w) / abs(B4$delta)
B4$sign_disagree <- sign(B4$nested_w) != sign(B4$shapley_w)
s4 <- data.frame(median_rel_gap = median(B4$gap_rel),
                 p95_rel_gap = quantile(B4$gap_rel, .95),
                 share_sign_disagreement = mean(B4$sign_disagree))
write.csv(s4, file.path(DIR_TAB, "simB4_nested_vs_shapley.csv"), row.names = FALSE)
log_msg("B4 done")

# ---- B5: small-sample inference in the SRISK regressions --------------------------
vfile <- file.path(DIR_PROC, "validation_data_gfc.rds")
if (file.exists(vfile)) {
  d <- readRDS(vfile)
  base <- lm(y ~ lfc + lpre, data = d)
  sig <- summary(base)$sigma; yhat <- fitted(base)
  # (a) size under the null: no network effect at all
  sizes <- t(replicate(REPS, {
    dd <- d; dd$y <- yhat + rnorm(nrow(d), 0, sig)
    f <- y ~ lfc + lpre + e_rec
    c(HC1 = abs(robust_ols(f, dd, "HC1")$t["e_rec"]) > qt(.975, nrow(d) - 4),
      HC3 = abs(robust_ols(f, dd, "HC3")$t["e_rec"]) > qt(.975, nrow(d) - 4),
      FL  = freedman_lane(f, "e_rec", dd, B = 199)$p <= 0.05) }))
  # (b) confounding: the outcome depends on importance as a counterparty only,
  # with the expected t-statistic of size_cpty set to the paper's value (7.28)
  T_TARGET <- as.numeric(Sys.getenv("NCI_T_CPTY", 7.28))
  rz <- resid(lm(size_cpty ~ lfc + lpre, data = d))
  bc <- T_TARGET * sig / sqrt(sum(rz^2))
  spur <- replicate(REPS, {
    dd <- d; dd$y <- yhat + bc * (d$size_cpty - mean(d$size_cpty)) +
      rnorm(nrow(d), 0, sig)
    abs(robust_ols(y ~ lfc + lpre + e_rec, dd, "HC1")$t["e_rec"]) > qt(.975, nrow(d) - 4) })
  s5 <- data.frame(test = c("HC1", "HC3", "Freedman-Lane", "Spurious e_rec (HC1)"),
                   rejection = c(colMeans(sizes), mean(spur)),
                   nominal = 0.05,
                   corr_e_rec_size_cpty = cor(d$e_rec, d$size_cpty),
                   partial_corr = cor(resid(lm(e_rec ~ lfc + lpre, d)), rz))
  write.csv(s5, file.path(DIR_TAB, "simB5_validation_inference.csv"), row.names = FALSE)
  log_msg("B5 done")
}

# ---- B6: contagion on synthetic topologies --------------------------------------
DR <- as.integer(Sys.getenv("NCI_DRAWS_SIM", 30))
topo <- function(type, N = 19) {
  repeat {
    if (type == "ER") A <- (matrix(runif(N^2), N, N) < .35) * 1
    else { core <- 1:5; P <- matrix(.1, N, N); P[core, ] <- .5; P[, core] <- .5
           P[core, core] <- .9; A <- (matrix(runif(N^2), N, N) < P) * 1 }
    diag(A) <- 0; if (all(rowSums(A) > 0)) return(A)
  }
}
B6 <- do.call(rbind, lapply(c("ER", "Core-periphery"), function(tp) {
  A <- topo(tp); w <- sort(lognorm_w(19, 1), decreasing = TRUE)
  do.call(rbind, lapply(c(-20, 0, 20), function(b) {
    dr <- reciprocity_mcmc(A, w, b, n_draws = DR, burn = 1000, thin = 100)
    data.frame(topology = tp, beta = b,
               Gamma = mean(sapply(dr, function(X) gamF(w, X))),
               loss = mean(sapply(dr, function(X) debtrank_loss(X, w, .8)$total)))
  }))
}))
B6$pct_vs_random <- ave(B6$loss, B6$topology, FUN = function(x) 100 * (x / x[2] - 1))
write.csv(B6, file.path(DIR_TAB, "simB6_contagion_topologies.csv"), row.names = FALSE)
log_msg("B6 done")

# ---- summary table (LaTeX) ---------------------------------------------------------
f2 <- function(x) fmt_num(x, 2); f3 <- function(x) fmt_num(x, 3)
b1n <- B1[B1$delta == 0, ]
rows <- list(
  c("E1", "Permutation test, size at 5\\% (range over $N$, dispersion)",
    sprintf("%s--%s", f3(min(b1n$rej_perm)), f3(max(b1n$rej_perm)))),
  c("E1", "Mantel normal approximation, size at 5\\% (range)",
    sprintf("%s--%s", f3(min(b1n$rej_mantel)), f3(max(b1n$rej_mantel)))),
  c("E2", "Binary $\\Gamma_F$, mean at $\\theta=0.2\\%$ / $2\\%$",
    paste(f2(B2s$Gamma.mean[B2s$theta == 0.002 & B2s$spec == "binary"]),
          f2(B2s$Gamma.mean[B2s$theta == 0.02]), sep = " / ")),
  c("E2", "Spearman correlation of $\\Gamma_F$ across $\\theta=0.2\\%$ and $2\\%$",
    f2(cors["Gamma.0.002", "Gamma.0.02"])),
  c("E3", "Share of circuits in $\\psi_{TOT}$ at $N=8$ / $N=80$ (no structure)",
    paste(f3(B3s$circuits[B3s$N == 8]), f3(B3s$circuits[B3s$N == 80]), sep = " / ")),
  c("E4", "Nested vs Shapley weight effect: median / 95th pct gap ($\\%$ of $|\\Delta\\psi|$)",
    paste(f2(100 * s4$median_rel_gap), f2(100 * s4$p95_rel_gap), sep = " / ")),
  c("E4", "Share of sign disagreements in the weight effect", f2(s4$share_sign_disagreement)))
if (exists("s5")) rows <- c(rows, list(
  c("E5", "Rejection under the null: HC1 / HC3 / Freedman--Lane",
    paste(f3(s5$rejection[1]), f3(s5$rejection[2]), f3(s5$rejection[3]), sep = " / ")),
  c("E5", "Spurious rejection of reciprocal emb.\\ when only counterparty size matters",
    f3(s5$rejection[4]))))
rows <- c(rows, list(
  c("E6", "Loss change at highest tilt: ER / core--periphery ($\\%$)",
    paste(f2(B6$pct_vs_random[B6$topology == "ER" & B6$beta == 20]),
          f2(B6$pct_vs_random[B6$topology == "Core-periphery" & B6$beta == 20]),
          sep = " / "))))
tabS2 <- as.data.frame(do.call(rbind, rows)); names(tabS2) <- c("Experiment", "Quantity", "Value")
write.csv(tabS2, file.path(DIR_TAB, "tabS2_robustness_summary.csv"), row.names = FALSE)  # E7-E8 appended by script 18
write_tex_table(tabS2, "tabS2_robustness_summary.tex",
  caption = "Monte Carlo robustness experiments: summary",
  label = "tab:sim-robust", align = "llc", header = names(tabS2),
  notes = paste0(REPS, " replications per cell; see scripts/08\\_sim\\_robustness.R."))
log_msg("robustness summary written")
