# =============================================================================
# 16_bank_null_robustness.R -- Section 5.4 / Appendix B.4: is the alignment of
# mutual exposures with size beyond connectivity robust to the threshold and to
# the sample of reporting countries?
#
# For the five quarters of Table 6, thresholds theta in {1%, 2%, 5%} and two
# samples (the 19 continuous reporters; 22 = 19 + Canada, Finland, India, which
# also report bilateral ultimate-risk claims since 2006), the reciprocal
# component psi_F^(1) is tested against the directed binary configuration model
# (Definition 1(b)), with the weights held fixed.
# Reads the raw CBS zip (data/raw/bis/WS_CBS_PUB_csv_col.zip), ~3 minutes.
# Output: tabOA_bank_null_robustness.tex/csv; bis_fixed_support.csv (links only on
#         the bilateral cells published in every quarter); tabOA_bank_missing.tex/csv
#         (bound on unpublished cells, without Switzerland / the Netherlands)
# =============================================================================
source("R/utils.R"); source("R/data_loaders.R"); source("R/core_indices.R")
source("R/null_models.R"); source("R/fixes.R")
suppressPackageStartupMessages(library(data.table))
B_NULL <- as.integer(Sys.getenv("NCI_B_BANKNULL", 499))

cbs_zip <- list.files(file.path(DIR_RAW, "bis"), "^WS_CBS_PUB.*\\.zip$", full.names = TRUE)
if (!length(cbs_zip)) {
  log_msg("no CBS zip in data/raw/bis: bank null robustness skipped")
} else {
# same filter as 00_get_data.R
BIS_FILTER <- list(L_MEASURE = "S", CBS_BANK_TYPE = "4R", CBS_BASIS = "U", L_POSITION = "C",
                   L_INSTR = "A", REM_MATURITY = "A", CURR_TYPE_BOOK = "TO1", L_CP_SECTOR = "A")
TOTAL_CP <- "5J"
long <- suppressWarnings(read_bis_wide(cbs_zip[1], BIS_FILTER, from_year = 2006))
R19 <- c("AU","AT","BE","CL","FR","DE","GR","IE","IT","JP","NL","PT","ES","SE","CH","TW","TR","GB","US")
R22 <- sort(c(R19, "CA", "FI", "IN"))
QS  <- c("2007-01-01", "2009-10-01", "2012-10-01", "2019-10-01", "2025-10-01")
qlab <- function(d) paste0(format(as.Date(d), "%Y"), "-Q", (as.integer(format(as.Date(d), "%m")) - 1) %/% 3 + 1)

set.seed(3); out <- list()
for (smp in list(R19, R22)) {
  b <- bis_standardize(long, smp, TOTAL_CP)
  for (q in QS) {
    sn <- bis_snapshot(b$claims[!is.na(claims)], b$totals[!is.na(total_claims)], as.Date(q), smp)
    for (th in c(0.01, 0.02, 0.05)) {
      A <- bis_adjacency(sn, th); w <- sn$w
      G <- factorize(w, A * t(A))$Gamma
      db <- dbcm_test(A, function(X) psi_index(w, X * t(X)), B = B_NULL)
      out[[length(out) + 1]] <- data.frame(sample = length(smp), quarter = qlab(q), theta = th,
        density = sum(A) / (length(smp) * (length(smp) - 1)), Gamma_rec = G,
        p_dbcm = db$p, z_dbcm = db$z)
    }
  }
}
res <- do.call(rbind, out)
utils::write.csv(res, file.path(DIR_TAB, "tabOA_bank_null_robustness.csv"), row.names = FALSE)

# table: rows = sample x theta, columns = quarters, cells = Gamma (p; z)
cell <- function(r) sprintf("%.2f (%.3f; %.1f)", r$Gamma_rec, r$p_dbcm, r$z_dbcm)
rows <- list()
for (n in c(19, 22)) for (th in c(0.01, 0.02, 0.05)) {
  r <- res[res$sample == n & res$theta == th, ]
  rows[[length(rows) + 1]] <- data.frame(Sample = sprintf("%d countries", n),
    theta = sprintf("%g\\%%", 100 * th), t(sapply(seq_len(nrow(r)), function(i) cell(r[i, ]))))
}
tab <- do.call(rbind, rows); names(tab)[-(1:2)] <- qlab(QS)
write_tex_table(tab, "tabOA_bank_null_robustness.tex",
  caption = "Reciprocal exposures against the configuration model: thresholds and samples",
  label = "tab:oa-bank-null", align = "llccccc",
  header = c("Sample", "$\\theta$", qlab(QS)),
  notes = paste0("Cells: alignment factor of reciprocal exposures; in parentheses, $p$-value and $z$-score of $\\psi_F^{(1)}$ against ",
    B_NULL, " draws of the directed binary configuration model, weights held fixed. 22 countries: the 19 reporters plus Canada, Finland and India."))
print(res, row.names = FALSE)

# ---- fixed observable support ----------------------------------------------------
# The 19 reporters are the same in every quarter, but some bilateral cells are
# published only in part of the period (confidentiality); unpublished cells are
# treated as no link.  Here links are allowed only on the ordered pairs published
# in every quarter, so that changes cannot come from changes in publication.
cl0 <- readRDS(file.path(DIR_PROC, "bis_claims.rds")); tt0 <- readRDS(file.path(DIR_PROC, "bis_totals.rds"))
Q0 <- sort(unique(cl0$quarter)); Q0 <- Q0[Q0 >= as.Date("2007-01-01")]
pub <- cl0[quarter %in% Q0, .(nq = uniqueN(quarter)), by = .(creditor, counterparty)][nq == length(Q0)]
R0 <- sort(unique(tt0$creditor))
mask <- matrix(0, length(R0), length(R0), dimnames = list(R0, R0))
mask[cbind(match(pub$creditor, R0), match(pub$counterparty, R0))] <- 1
fs <- do.call(rbind, lapply(QS, function(q) {
  sn <- bis_snapshot(cl0, tt0, as.Date(q), R0)
  A_b <- bis_adjacency(sn, 0.02); sn$C <- sn$C * mask; A_f <- bis_adjacency(sn, 0.02)
  data.frame(quarter = qlab(q), Gamma_baseline = factorize(sn$w, A_b * t(A_b))$Gamma,
             Gamma_fixed_support = factorize(sn$w, A_f * t(A_f))$Gamma,
             reciprocal_pairs_baseline = sum(A_b * t(A_b)) / 2,
             reciprocal_pairs_fixed = sum(A_f * t(A_f)) / 2)
}))
utils::write.csv(fs, file.path(DIR_TAB, "bis_fixed_support.csv"), row.names = FALSE)
cat(sprintf("ordered pairs published in all %d quarters: %d of %d\n", length(Q0), sum(mask),
            length(R0) * (length(R0) - 1)))
print(fs, row.names = FALSE, digits = 3)
# ---- unpublished bilateral cells: bound and sensitivity ------------------------------
# An unpublished cell r -> c cannot exceed the unallocated part of the foreign claims
# of r (total minus the sum of all published country cells).  Upper scenario: an
# unpublished cell becomes a link when the reverse link exists and the unallocated
# share of r can accommodate it (at most floor(unallocated / theta) cells, heaviest
# counterparties first): the largest number of reciprocal pairs compatible with the
# published totals.  Also: the 19 reporters without Switzerland (reporting break in
# 2013-Q4, see below) and without the Netherlands (largest unallocated share).
TH <- 0.02; R19s <- sort(R19)
set.seed(4)
miss_q <- lapply(QS, function(q) {
  d <- long[quarter == as.Date(q) & L_REP_CTY %in% R19s]
  tot <- d[L_CP_COUNTRY == TOTAL_CP][match(R19s, L_REP_CTY), value]
  alloc <- vapply(R19s, function(r) d[L_REP_CTY == r & grepl("^[A-Z]{2}$", L_CP_COUNTRY) &
                                      L_CP_COUNTRY != r, sum(value)], numeric(1))
  unalloc <- pmax(0, 1 - alloc / tot)
  S <- matrix(NA_real_, 19, 19, dimnames = list(R19s, R19s))
  b <- d[L_CP_COUNTRY %in% R19s & L_CP_COUNTRY != L_REP_CTY]
  S[cbind(match(b$L_REP_CTY, R19s), match(b$L_CP_COUNTRY, R19s))] <- b$value / tot[match(b$L_REP_CTY, R19s)]
  diag(S) <- 0; w <- tot / sum(tot)
  A0 <- (S >= TH) * 1; A0[is.na(A0)] <- 0; A1 <- A0
  for (i in seq_along(R19s)) {
    cand <- which(is.na(S[i, ]) & A0[, i] == 1)
    k <- min(length(cand), floor(unalloc[i] / TH))
    if (k > 0) A1[i, cand[order(-w[cand])][seq_len(k)]] <- 1
  }
  st <- function(A) { db <- dbcm_test(A, function(X) psi_index(w, X * t(X)), B = B_NULL)
    c(G = factorize(w, A * t(A))$Gamma, p = db$p, pairs = sum(A * t(A)) / 2) }
  list(base = st(A0), upper = st(A1), unpublished = sum(is.na(S)),
       unalloc = setNames(unalloc, R19s))
})
drop_q <- function(drop) {
  smp <- setdiff(R19, drop); b <- bis_standardize(long, smp, TOTAL_CP)
  lapply(QS, function(q) {
    sn <- bis_snapshot(b$claims[!is.na(claims)], b$totals[!is.na(total_claims)], as.Date(q), smp)
    A <- bis_adjacency(sn, TH); w <- sn$w
    c(G = factorize(w, A * t(A))$Gamma,
      p = dbcm_test(A, function(X) psi_index(w, X * t(X)), B = B_NULL)$p)
  })
}
sens <- list("Baseline (19 countries)" = lapply(miss_q, `[[`, "base"),
             "Unpublished cells: upper scenario" = lapply(miss_q, `[[`, "upper"),
             "Without Switzerland" = drop_q("CH"),
             "Without the Netherlands" = drop_q("NL"),
             "Without Switzerland and the Netherlands" = drop_q(c("CH", "NL")))
ms <- do.call(rbind, lapply(names(sens), function(n) data.frame(case = n, quarter = qlab(QS),
  Gamma_rec = sapply(sens[[n]], `[[`, "G"), p_dbcm = sapply(sens[[n]], `[[`, "p"))))
utils::write.csv(ms, file.path(DIR_TAB, "tabOA_bank_missing.csv"), row.names = FALSE)
tm <- do.call(rbind, lapply(names(sens), function(n) {
  r <- ms[ms$case == n, ]
  data.frame(Case = n, t(sprintf("%.2f (%.3f)", r$Gamma_rec, r$p_dbcm)))
}))
names(tm)[-1] <- qlab(QS)
write_tex_table(tm, "tabOA_bank_missing.tex",
  caption = "Reciprocal exposures: unpublished bilateral cells and reporting countries",
  label = "tab:oa-bank-missing", align = "lccccc", header = c("", qlab(QS)),
  notes = paste0("Cells: alignment factor of reciprocal exposures ($\\theta=2\\%$); in parentheses, $p$-value of $\\psi_F^{(1)}$ against ",
    B_NULL, " draws of the directed binary configuration model. Upper scenario: every unpublished cell becomes a link when the reverse link exists and the unallocated share of the creditor's foreign claims can accommodate it."))
ua <- sapply(miss_q, `[[`, "unalloc"); colnames(ua) <- qlab(QS)
cat("\nunallocated share of foreign claims (upper bound for any unpublished cell), %:\n")
print(round(100 * ua[order(-ua[, 4]), ], 1))
cat(sprintf("unpublished cells among the 342: %s\n", paste(sapply(miss_q, `[[`, "unpublished"), collapse = " ")))
print(tm, row.names = FALSE)
# Swiss reporting break: unallocated share before and after 2013-Q4
ch <- long[L_REP_CTY == "CH", .(total = value[L_CP_COUNTRY == TOTAL_CP][1],
           alloc = sum(value[grepl("^[A-Z]{2}$", L_CP_COUNTRY) & L_CP_COUNTRY != "CH"])), by = quarter]
ch[, unalloc := 1 - alloc / total]
cat(sprintf("Switzerland, unallocated share: %.0f%% in 2013-Q3, %.0f%% in 2013-Q4\n",
            100 * ch[quarter == as.Date("2013-07-01"), unalloc], 100 * ch[quarter == as.Date("2013-10-01"), unalloc]))

cat(sprintf("z > 0 in %d of %d cases; p-values between %.3f and %.3f\n",
            sum(res$z_dbcm > 0), nrow(res), min(res$p_dbcm), max(res$p_dbcm)))
log_msg("bank null robustness written")
}
