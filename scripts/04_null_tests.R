# =============================================================================
# 04_null_tests.R -- Table 4: permutation and degree-preserving null models
# =============================================================================
source("R/utils.R"); source("R/core_indices.R"); source("R/data_loaders.R")
source("R/null_models.R"); source("R/fixes.R"); source("R/contagion.R")

B <- as.integer(Sys.getenv("NCI_B", 1000))       # draws for both nulls
set.seed(123)

null_row <- function(w, A, label) {
  sd <- structural_decomposition(w, A, H = 2, alpha = c(.5, .5))
  pF1 <- perm_alignment(w, sd$mats$F1, B)
  pF2 <- perm_alignment(w, sd$mats$F2, B)
  mF1 <- mantel_moments(w, sd$mats$F1)
  rw  <- rewiring_test(A, structural_stats(w), B)
  # soft-constraint alternative: directed binary configuration model
  db  <- dbcm_test(A, function(X) structural_stats(w)(X)[["psi_F1"]], B)
  data.frame(network = label,
    G_rec = pF1$Gamma, p_perm_rec = pF1$p_share, z_mantel_rec = mF1$z,
    z_avail_rec = rw$z[["mbar_F1"]], p_rew_rec = rw$p_share[["psi_F1"]],
    G_c3 = pF2$Gamma, p_perm_c3 = pF2$p_share, p_rew_c3 = rw$p_share[["psi_F2"]],
    p_rew_oneway = rw$p_share[["psi_R1"]],
    p_dbcm_rec = db$p, z_dbcm_rec = db$z)
}

# ---- production network (binary, theta = 0.5%) ---------------------------------
avail <- sort(as.integer(sub(paste0("^", PROD_SOURCE, "_([0-9]{4})\\.rds$"), "\\1",
  list.files(DIR_PROC, paste0("^", PROD_SOURCE, "_[0-9]{4}\\.rds$")))))
NULL_YEARS <- intersect(if (PROD_SOURCE == "icio") c(1995, 2007, 2014, 2022)
                        else c(2000, 2007, 2014), avail)
out <- list()
for (yr in NULL_YEARS) {
  dc <- wiod_country(load_wiod_year(yr)); w <- dc$VA / sum(dc$VA)
  A <- (wiod_intensity(dc, "share") >= 0.005) * 1
  out[[length(out) + 1]] <- cbind(panel = "World production network",
                                  null_row(w, A, as.character(yr)))
  log_msg("null models WIOD ", yr)
}
# ---- banking network (binary, theta = 2%) ---------------------------------------
bis <- readRDS(file.path(DIR_PROC, "results_bis.rds"))
for (q in c("2007-Q1", "2009-Q4", "2012-Q4", "2019-Q4", "2025-Q4")) {
  m <- bis$mats[[q]]
  out[[length(out) + 1]] <- cbind(panel = "International banking network",
                                  null_row(m$w, m$A, q))
  log_msg("null models BIS ", q)
}
res <- do.call(rbind, out)
saveRDS(res, file.path(DIR_PROC, "results_null.rds"))

pfmt <- function(p) ifelse(p < 0.001, "$<0.001$", formatC(p, format = "f", digits = 3))
tab4 <- data.frame(Date = res$network, G = fmt_num(res$G_rec, 2),
  pp = pfmt(res$p_perm_rec), z = fmt_num(res$z_avail_rec, 1),
  pr = pfmt(res$p_rew_rec), G3 = fmt_num(res$G_c3, 2),
  pp3 = pfmt(res$p_perm_c3), pr3 = pfmt(res$p_rew_c3))
write_tex_table(tab4, "tab4_null_models.tex", panel = res$panel,
  caption = "Structural components against permutation and degree-preserving null models",
  label = "tab:null", align = "lccccccc",
  header = c("", "$\\Gamma$", "$p$, perm.", "$z$, avail.", "$p$, rewiring",
             "$\\Gamma$", "$p$, perm.", "$p$, rewiring"),
  notes = paste0("Columns 2--5: reciprocal component; columns 6--8: circuits of",
    " length 3. Binary networks with $\\theta=0.5\\%$ (production) and",
    " $\\theta=2\\%$ (banking). ``$p$, perm.'': share of ", B, " random",
    " permutations of the weights with alignment factor at least as large as",
    " the observed one; under this null $E(\\Gamma)=1$ (Proposition 1).",
    " ``$z$, avail.'': standardized difference between the observed availability",
    " of reciprocal links and its mean over ", B, " directed rewirings preserving",
    " in- and out-degrees (double-edge swaps and 3-cycle reversals).",
    " ``$p$, rewiring'': share of rewired networks whose component $\\psi$ is",
    " at least as large as the observed one, with the weights held fixed."))
# full results (Online Appendix OA.2), including analytic Mantel z-scores
write.csv(res, file.path(DIR_TAB, "tabOA2_null_models_full.csv"), row.names = FALSE)
log_msg("Table 4 written")

# ---- power of the permutation test at the observed weights (WIOD 2014) ---------
# Alternatives: networks with the 2014 undirected support and number of
# reciprocal pairs, tilted towards higher alignment of reciprocal links.
# Question: would the test have detected the 2000 level of alignment in 2014?
Y_LAST <- max(NULL_YEARS); Y_FIRST <- min(NULL_YEARS)
dc14 <- wiod_country(load_wiod_year(Y_LAST)); w14 <- dc14$VA / sum(dc14$VA)
A14 <- (wiod_intensity(dc14, "share") >= 0.005) * 1
G2000 <- res$G_rec[res$network == as.character(Y_FIRST)]
pw <- do.call(rbind, lapply(c(0, 5, 10, 20, 40), function(b) {
  dr <- reciprocity_mcmc(A14, w14, b, n_draws = 40, burn = 1000, thin = 100,
                         require_rows = FALSE)
  do.call(rbind, lapply(dr, function(X) { pa <- perm_alignment(w14, X * t(X), B = 199)
    data.frame(beta = b, Gamma = pa$Gamma, reject = pa$p_plus1 <= 0.05) }))
}))
pw$bin <- cut(pw$Gamma, breaks = unique(quantile(pw$Gamma, seq(0, 1, .1))),
              include.lowest = TRUE)
pc_ <- aggregate(cbind(Gamma, reject) ~ bin, data = pw, FUN = mean)
pc_ <- pc_[order(pc_$Gamma), ]
near <- pw[abs(pw$Gamma - G2000) <= 0.25 * G2000, ]
pow_at <- if (nrow(near) > 0) mean(near$reject) else NA
write.csv(pc_, file.path(DIR_TAB, sprintf("power_%s%d_curve.csv", PROD_SOURCE, Y_LAST)),
          row.names = FALSE)
writeLines(sprintf("Power at Gamma(%d) = %.2f with %d weights: %s (%d networks within 25%%)",
                   Y_FIRST, G2000, Y_LAST, fmt_num(pow_at, 2), nrow(near)),
           file.path(DIR_TAB, sprintf("power_%s%d.txt", PROD_SOURCE, Y_LAST)))
pP <- ggplot(pc_, aes(Gamma, reject)) +
  geom_hline(yintercept = .05, colour = MUTED, linewidth = .3, linetype = "22") +
  geom_vline(xintercept = G2000, colour = COL_MECH[["Reciprocal"]], linewidth = .4) +
  annotate("text", x = G2000, y = .55, hjust = -.05, size = 2.3,
           colour = COL_MECH[["Reciprocal"]], label = sprintf("Γ in %d = %.2f", Y_FIRST, G2000)) +
  geom_line(colour = INK, linewidth = .7) + geom_point(colour = INK, size = 2) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Alignment of reciprocal links in the simulated network, Γ",
       y = "Rejection rate at 5%",
       title = sprintf("Power of the permutation test with the %d value-added weights", Y_LAST),
       subtitle = sprintf("Networks with the %d support and number of reciprocal pairs, tilted towards alignment", Y_LAST)) +
  theme_paper()
save_fig(pP, sprintf("figOA_power_%s%d", PROD_SOURCE, Y_LAST), width = 6, height = 3.6)
log_msg("power analysis written")
