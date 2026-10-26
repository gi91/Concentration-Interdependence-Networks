# =============================================================================
# 06_contagion.R -- Section 7.4: contagion experiment, Figure 4
# =============================================================================
source("R/utils.R"); source("R/core_indices.R"); source("R/null_models.R")
source("R/contagion.R")

DRAWS <- as.integer(Sys.getenv("NCI_DRAWS", 150))     # networks per tilt
BETAS <- c(-40, -20, -10, -5, 0, 5, 10, 20, 40)        # tilt on Gamma_F (tune: see acceptance)
LEVS  <- c(0.5, 0.8, 0.95)                            # interbank exposure / equity
DATES <- c("2007-Q2", "2012-Q4")
B_REW <- as.integer(Sys.getenv("NCI_B", 1000))
bis <- readRDS(file.path(DIR_PROC, "results_bis.rds"))

# proportional allocation variant: claims spread in proportion to the size of
# the counterparties (robustness mentioned in the paper)
loss_prop <- function(X, w, c) {
  Xp <- sweep(X, 2, w, "*"); Xp <- Xp / rowSums(Xp)
  Lam <- c * Xp; N <- nrow(X); H0 <- diag(N); H <- H0
  repeat { Hn <- pmin(H0 + Lam %*% H, 1); if (max(abs(Hn - H)) < 1e-12) break
    H <- Hn }
  L <- vapply(seq_len(N), function(j) sum(w[-j] * Hn[-j, j]) / sum(w[-j]), 0)
  sum(w * L)
}

net_stats <- function(X, w) {
  sd <- structural_decomposition(w, X, H = 2, alpha = c(.5, .5)); tb <- sd$table
  top5 <- order(w, decreasing = TRUE)[1:5]
  c(GammaF = tb$Gamma[tb$comp == "F" & tb$h == 1],
    psiR1 = tb$psi[tb$comp == "R" & tb$h == 1],
    psiF1 = tb$psi[tb$comp == "F" & tb$h == 1],
    psiF2 = tb$psi[tb$comp == "F" & tb$h == 2],
    setNames(sapply(LEVS, function(c) debtrank_loss(X, w, c)$total),
             paste0("loss_", LEVS)),
    loss_top5_0.8 = { L <- debtrank_loss(X, w, 0.8)$L
                      sum(w[top5] * L[top5]) / sum(w[top5]) },
    loss_prop_0.8 = loss_prop(X, w, 0.8))
}

all <- list(); obs <- list(); pct <- list()
set.seed(2007)
for (q in DATES) {
  m <- bis$mats[[q]]; X0 <- m$A; w <- m$w
  obs[[q]] <- data.frame(date = q, t(net_stats(X0, w)))
  for (b in BETAS) {
    draws <- reciprocity_mcmc(X0, w, beta = b, n_draws = DRAWS)
    st <- t(vapply(draws, net_stats, numeric(ncol(obs[[q]]) - 1), w = w))
    all[[paste(q, b)]] <- data.frame(date = q, beta = b, st,
                                     acc = attr(draws, "acceptance"))
    log_msg("contagion ", q, " beta = ", b)
  }
  # observed network within a degree-preserving ensemble (in/out degrees fixed)
  rew <- vapply(seq_len(B_REW), function(i)
    debtrank_loss(rewire_directed(X0), w, 0.8)$total, numeric(1))
  pct[[q]] <- data.frame(date = q,
                         percentile = mean(rew <= obs[[q]]$loss_0.8))
}
sims <- do.call(rbind, all); obs <- do.call(rbind, obs); pct <- do.call(rbind, pct)
saveRDS(list(sims = sims, obs = obs, pct = pct),
        file.path(DIR_PROC, "results_contagion.rds"))

# sanity checks: invariants of the chain
chk <- sapply(DATES, function(q) {
  X0 <- bis$mats[[q]]$A
  d <- reciprocity_mcmc(X0, bis$mats[[q]]$w, beta = 2, n_draws = 20, burn = 200,
                        thin = 50)
  c(same_support = all(sapply(d, function(X)
      identical(support_undirected(X), support_undirected(X0)))),
    same_reciprocal_pairs = all(sapply(d, function(X)
      sum(X * t(X)) == sum(X0 * t(X0)))),
    spectral_radius_c = all(sapply(d, function(X)
      abs(max(Mod(eigen(debtrank_lambda(X, 0.8))$values)) - 0.8) < 1e-8)))
})
print(chk)

# ---- summaries and Figure 4 -----------------------------------------------------
summ <- do.call(rbind, lapply(split(sims, list(sims$date, sims$beta)), function(s)
  data.frame(date = s$date[1], beta = s$beta[1], GammaF = mean(s$GammaF),
             t(colMeans(s[, grep("^loss", names(s))])))))
base <- summ[summ$beta == 0, ]
# per-draw % change w.r.t. the untilted mean, then mean and 95% band per tilt
fig <- do.call(rbind, lapply(LEVS, function(c) {
  v <- paste0("loss_", c)
  d <- sims; d$pct <- 100 * (d[[v]] / base[[v]][match(d$date, base$date)] - 1)
  do.call(rbind, lapply(split(d, list(d$date, d$beta)), function(z)
    data.frame(date = z$date[1], beta = z$beta[1], GammaF = mean(z$GammaF), c = c,
               pct = mean(z$pct), se = sd(z$pct) / sqrt(nrow(z)))))
}))
fo <- do.call(rbind, lapply(LEVS, function(c) {
  v <- paste0("loss_", c)
  data.frame(date = obs$date, GammaF = obs$GammaF, c = c,
             pct = 100 * (obs[[v]] / base[[v]][match(obs$date, base$date)] - 1))
}))
lab_c <- function(c) paste0("c = ", c)
fig$lev <- factor(lab_c(fig$c), lab_c(LEVS)); fo$lev <- factor(lab_c(fo$c), lab_c(LEVS))
fig <- fig[order(fig$date, fig$c, fig$GammaF), ]
COL_LEV <- setNames(SEQ3, lab_c(LEVS))
ann <- merge(pct, obs[, c("date", "GammaF")], by = "date")
ann$lab <- sprintf("Observed network: %.0fth percentile\nof degree-preserving ensemble (c = 0.8)",
                   100 * ann$percentile)
p4 <- ggplot(fig, aes(GammaF, pct, colour = lev, fill = lev)) +
  geom_hline(yintercept = 0, colour = MUTED, linewidth = .3, linetype = "22") +
  geom_ribbon(aes(ymin = pct - 1.96 * se, ymax = pct + 1.96 * se), colour = NA,
              alpha = .18) +
  geom_line(linewidth = .7) + geom_point(size = 1.8) +
  geom_point(data = fo, shape = 21, size = 3.4, stroke = 1, fill = "white") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = -0.03, vjust = 1.2, size = 2.2, colour = INK2, lineheight = .9) +
  facet_wrap(~ paste("Network of", date), scales = "free_x") +
  scale_colour_manual(values = COL_LEV, name = "Interbank exposure / equity") +
  scale_fill_manual(values = COL_LEV, name = "Interbank exposure / equity") +
  labs(x = expression("Alignment of reciprocal exposures with size, "*Gamma),
       y = "Contagion loss, % change vs. random reciprocity",
       title = "Same baseline index, density and reciprocity: losses rise with alignment",
       subtitle = "Mean over networks at each tilt of the chain, 95% band; open circles: observed networks",
       caption = "Linear DebtRank, failure of each system in turn, losses weighted by size") +
  theme_paper() + theme(legend.title = element_text(size = 7, colour = INK2))
save_fig(p4, "fig4_contagion", width = 7, height = 4)

# regression of losses on the structural components (c = 0.8).
# NB: with the undirected support fixed, psi_R^(1) = (psi_D - psi_F^(1)) / 2
# for binary networks, so one-way and reciprocal components are EXACTLY
# collinear along the chain and cannot enter the same regression.  We report
# the correlation and fit reciprocal + circuits only.
reg <- do.call(rbind, lapply(DATES, function(q) {
  s <- sims[sims$date == q, ]; z <- function(x) (x - mean(x)) / sd(x)
  f <- lm(log(loss_0.8) ~ z(psiF1) + z(psiF2), data = s)
  data.frame(date = q, cor_R1_F1 = cor(s$psiR1, s$psiF1),
             pct_per_sd_F1 = 100 * coef(f)[2], pct_per_sd_F2 = 100 * coef(f)[3])
}))
summ_out <- merge(summ, pct, by = "date")
write.csv(summ_out, file.path(DIR_TAB, "fig4_contagion_summary.csv"), row.names = FALSE)
write.csv(reg, file.path(DIR_TAB, "contagion_regression.csv"), row.names = FALSE)
log_msg("Figure 4 written; acceptance rates: ",
        paste(round(range(sims$acc), 2), collapse = "-"))
