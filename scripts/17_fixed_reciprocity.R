# =============================================================================
# 17_fixed_reciprocity.R -- tests with fixed reciprocity (Definition S1(c), (d))
#
# (c) networks with the observed in- and out-degrees and number of reciprocal
#     pairs; (d) the same degrees and, for every node, the same number of
#     reciprocal partners.  Sampler: R/constrained_nulls.R (C++ via Rcpp).
#
#   1. check of the sampler by complete enumeration of 5-node networks
#   2. Main-text Table 4 and Table S11: "p, fixed recip." and "p, fixed partners"
#      (4 chains x 1,000 draws), production (theta = 0.5%) and banking (2%)
#   3. Table S12, lower panel: banking, 19/22 countries x theta in {1, 2, 5}%,
#      fixed reciprocal pairs (4 chains x 500 draws; the 19 / 2% row is the test of Table 4)
# Rows (12c)-(12d) of Table S4 are computed in 05_srisk_validation.R.
# Run after 03, 04 and 16.  NCI_CORES sets the number of parallel chains
# (Linux/macOS; on Windows the chains run one after the other).
# Output: fixed_recip_check_enumeration.csv, fixed_recip_tests.csv,
#         tab4_null_models.tex (adds the column), tabOA2_null_models_full.csv
#         (adds four columns), tabOA_bank_null_fixed.csv
# =============================================================================
source("R/utils.R"); source("R/data_loaders.R"); source("R/core_indices.R")
source("R/null_models.R"); source("R/fixes.R"); source("R/constrained_nulls.R")
suppressPackageStartupMessages(library(data.table))
cn_compile()
N_DRAWS <- as.integer(Sys.getenv("NCI_B_FIXED", 1000))   # draws per chain (x 4 chains)

# ---- 1. complete enumeration on 5-node networks --------------------------------
# All 2^20 directed networks on 5 nodes; for random starting networks the chain
# must visit only the constrained set and, on the part it reaches, every network
# with the same frequency.
N5 <- 5
pos <- expand.grid(j = 1:N5, i = 1:N5)[, c("i", "j")]; pos <- pos[pos$i != pos$j, ]
codes <- 0:(2^20 - 1)
bits <- sapply(0:19, function(b) bitwAnd(bitwShiftR(codes, b), 1L))
kout5 <- sapply(1:N5, function(v) rowSums(bits[, pos$i == v, drop = FALSE]))
kin5  <- sapply(1:N5, function(v) rowSums(bits[, pos$j == v, drop = FALSE]))
recm <- bits * bits[, match(paste(pos$j, pos$i), paste(pos$i, pos$j))]
part5 <- sapply(1:N5, function(v) rowSums(recm[, pos$i == v, drop = FALSE])); np5 <- rowSums(recm) / 2
rm(recm)
set.seed(11); enum <- list()
for (r in 1:6) {
  repeat { A <- (matrix(runif(25), 5, 5) < .45) * 1; diag(A) <- 0; if (sum(A) >= 6) break }
  deg_ok <- rowSums(sweep(kout5, 2, rowSums(A)) != 0) == 0 & rowSums(sweep(kin5, 2, colSums(A)) != 0) == 0
  for (nl in c("pairs", "partners")) {
    S <- codes[deg_ok & (if (nl == "pairs") np5 == sum(A * t(A)) / 2
                         else rowSums(sweep(part5, 2, rowSums(A * t(A))) != 0) == 0)]
    x <- constrained_chain(A, rep(.2, 5), nl, n_draws = 20000, burn_mult = 200,
                           thin_mult = 200, keep_codes = TRUE)
    tab <- table(x$codes)
    enum[[length(enum) + 1]] <- data.frame(network = r, null = nl, constrained_set = length(S),
      reached = length(tab), outside_set = sum(!as.numeric(names(tab)) %in% S),
      chi2_p = if (length(tab) > 1) suppressWarnings(chisq.test(as.numeric(tab))$p.value) else NA)
  }
}
enum <- do.call(rbind, enum)
write.csv(enum, file.path(DIR_TAB, "fixed_recip_check_enumeration.csv"), row.names = FALSE)
print(enum, row.names = FALSE, digits = 3)
stopifnot(all(enum$outside_set == 0))
log_msg(sprintf("enumeration check: no draw outside the constrained set; chi-square p on the reached set between %.2f and %.2f",
                min(enum$chi2_p, na.rm = TRUE), max(enum$chi2_p, na.rm = TRUE)))

# ---- 2. Table 4 and Table S11 ----------------------------------------------------
nets <- list()
for (yr in c(1995, 2007, 2014, 2022)) {
  dc <- wiod_country(load_wiod_year(yr))
  nets[[as.character(yr)]] <- list(panel = "World production network", w = dc$VA / sum(dc$VA),
                                   A = (wiod_intensity(dc, "share") >= 0.005) * 1)
}
bis <- readRDS(file.path(DIR_PROC, "results_bis.rds"))
for (q in c("2007-Q1", "2009-Q4", "2012-Q4", "2019-Q4", "2025-Q4"))
  nets[[q]] <- list(panel = "International banking network", w = bis$mats[[q]]$w, A = bis$mats[[q]]$A)
fx <- do.call(rbind, lapply(seq_along(nets), function(k) {
  n <- nets[[k]]
  tc <- constrained_test(n$A, n$w, "pairs", n_draws = N_DRAWS, seed = 100 + k)
  td <- constrained_test(n$A, n$w, "partners", n_draws = N_DRAWS, seed = 200 + k)
  log_msg("fixed reciprocity: ", names(nets)[k])
  data.frame(panel = n$panel, network = names(nets)[k], links = sum(n$A),
             reciprocal_pairs = sum(n$A * t(n$A)) / 2,
             p_fixed_pairs = tc$p, p_fixed_pairs_min = tc$p_range[1], p_fixed_pairs_max = tc$p_range[2],
             p_fixed_partners = td$p, p_fixed_partners_min = td$p_range[1],
             p_fixed_partners_max = td$p_range[2], draws = tc$B)
}))
write.csv(fx, file.path(DIR_TAB, "fixed_recip_tests.csv"), row.names = FALSE)
print(fx[, c("network", "p_fixed_pairs", "p_fixed_pairs_min", "p_fixed_pairs_max",
             "p_fixed_partners", "p_fixed_partners_min", "p_fixed_partners_max")], digits = 3,
      row.names = FALSE)

# Table 4 with the new column (results of 04_null_tests.R + fixed pairs)
res <- readRDS(file.path(DIR_PROC, "results_null.rds"))
res <- merge(res, fx[, c("network", "p_fixed_pairs", "p_fixed_partners", "p_fixed_pairs_min",
                         "p_fixed_pairs_max", "p_fixed_partners_min", "p_fixed_partners_max")],
             by = "network", sort = FALSE)
res <- res[match(names(nets), res$network), ]
write.csv(res, file.path(DIR_TAB, "tabOA2_null_models_full.csv"), row.names = FALSE)
pfmt <- function(p) ifelse(p < 0.001, "$<0.001$", formatC(p, format = "f", digits = 3))
tab4 <- data.frame(Date = res$network, G = fmt_num(res$G_rec, 2), pp = pfmt(res$p_perm_rec),
  z = fmt_num(res$z_avail_rec, 1), pr = pfmt(res$p_rew_rec), pf = pfmt(res$p_fixed_pairs),
  G3 = fmt_num(res$G_c3, 2), pp3 = pfmt(res$p_perm_c3), pr3 = pfmt(res$p_rew_c3))
write_tex_table(tab4, "tab4_null_models.tex", panel = res$panel,
  caption = "Structural components against permutation and degree-preserving null models",
  label = "tab:null", align = "lcccccccc",
  header = c("", "$\\Gamma$", "$p$, perm.", "$z$, avail.", "$p$, rewiring", "$p$, fixed recip.",
             "$\\Gamma$", "$p$, perm.", "$p$, rewiring"),
  notes = paste0("``$p$, fixed recip.'': test of $\\psi_F^{(1)}$ against ", 4 * N_DRAWS,
    " networks, in four chains, with the observed in- and out-degrees and number of reciprocal",
    " pairs (Definition S1(c)). Other columns as in 04_null_tests.R."))

# ---- 3. Table S12, lower panel ---------------------------------------------------
cbs_zip <- list.files(file.path(DIR_RAW, "bis"), "^WS_CBS_PUB.*\\.zip$", full.names = TRUE)
if (length(cbs_zip)) {
  BIS_FILTER <- list(L_MEASURE = "S", CBS_BANK_TYPE = "4R", CBS_BASIS = "U", L_POSITION = "C",
                     L_INSTR = "A", REM_MATURITY = "A", CURR_TYPE_BOOK = "TO1", L_CP_SECTOR = "A")
  long <- suppressWarnings(read_bis_wide(cbs_zip[1], BIS_FILTER, from_year = 2006))
  R19 <- c("AU","AT","BE","CL","FR","DE","GR","IE","IT","JP","NL","PT","ES","SE","CH","TW","TR","GB","US")
  R22 <- sort(c(R19, "CA", "FI", "IN"))
  QS  <- c("2007-01-01", "2009-10-01", "2012-10-01", "2019-10-01", "2025-10-01")
  qlab <- function(d) paste0(format(as.Date(d), "%Y"), "-Q", (as.integer(format(as.Date(d), "%m")) - 1) %/% 3 + 1)
  out <- list(); k <- 0
  for (smp in list(R19, R22)) {
    b <- bis_standardize(long, smp, "5J")
    for (q in QS) {
      sn <- bis_snapshot(b$claims[!is.na(claims)], b$totals[!is.na(total_claims)], as.Date(q), smp)
      for (th in c(0.01, 0.02, 0.05)) {
        k <- k + 1; A <- bis_adjacency(sn, th)
        if (length(smp) == 19 && th == 0.02) {
          # baseline row: the same network as Table 4, so the same test (4 x N_DRAWS draws)
          A2 <- nets[[qlab(q)]]$A
          stopifnot(all(A[rownames(A2), colnames(A2)] == A2))
          r <- fx[fx$network == qlab(q), ]
          out[[k]] <- data.frame(sample = 19, quarter = qlab(q), theta = th,
                                 p_fixed_pairs = r$p_fixed_pairs, p_min = r$p_fixed_pairs_min,
                                 p_max = r$p_fixed_pairs_max, draws = r$draws)
          next
        }
        tc <- constrained_test(A, sn$w, "pairs", n_draws = N_DRAWS / 2, seed = 300 + k)
        out[[k]] <- data.frame(sample = length(smp), quarter = qlab(q), theta = th,
                               p_fixed_pairs = tc$p, p_min = tc$p_range[1], p_max = tc$p_range[2],
                               draws = tc$B)
      }
    }
  }
  s12 <- do.call(rbind, out)
  write.csv(s12, file.path(DIR_TAB, "tabOA_bank_null_fixed.csv"), row.names = FALSE)
  print(reshape(s12[, c("sample", "theta", "quarter", "p_fixed_pairs")], idvar = c("sample", "theta"),
                timevar = "quarter", direction = "wide"), row.names = FALSE, digits = 3)
  r2 <- s12[s12$theta == 0.02, "p_fixed_pairs"]; r15 <- s12[s12$theta != 0.02, "p_fixed_pairs"]
  cat(sprintf("fixed pairs: p at 2%% between %.3f and %.3f; at 1%% and 5%% between %.3f and %.3f\n",
              min(r2), max(r2), min(r15), max(r15)))
}
log_msg("fixed-reciprocity tests written")
