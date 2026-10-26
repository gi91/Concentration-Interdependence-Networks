# =============================================================================
# 00_make_synthetic_data.R
# Synthetic data with EXACTLY the structure of the processed inputs, so the
# whole pipeline can be run and checked end-to-end before the real WIOD, BIS
# and SRISK files are plugged in.  Numbers produced from these data are
# meaningless.  Do NOT run this script once real data are in data/processed.
# =============================================================================
source("R/utils.R"); source("R/data_loaders.R")
set.seed(20261008)

# ---- WIOD-like production network -------------------------------------------
wiod_cc <- c("AUS","AUT","BEL","BGR","BRA","CAN","CHE","CHN","CYP","CZE","DEU",
             "DNK","ESP","EST","FIN","FRA","GBR","GRC","HRV","HUN","IDN","IND",
             "IRL","ITA","JPN","KOR","LTU","LUX","LVA","MEX","MLT","NLD","NOR",
             "POL","PRT","ROU","RUS","SVK","SVN","SWE","TUR","TWN","USA")
J  <- 8                                     # industries (56 in WIOD)
nC <- length(wiod_cc); nCI <- nC * J
country  <- rep(wiod_cc, each = J)
industry <- rep(sprintf("I%02d", seq_len(J)), nC)
base <- setNames(rlnorm(nC, 0, 1.1), wiod_cc)
base[c("USA","JPN","DEU","CHN","GBR","FRA","ITA")] <-
  c(40, 18, 12, 4, 9, 8, 6)
noise_fixed <- matrix(rlnorm(nCI^2, 0, 1), nCI, nCI)     # persistent pair noise
for (yr in 2000:2014) {
  s <- base
  s["CHN"] <- 4 * exp(0.10 * (yr - 2000))                # China rises
  s <- s / sum(s)
  attract <- s^0.8; attract["CHN"] <- attract["CHN"] * (1 + 0.08 * (yr - 2000))
  go <- rep(s, each = J) * rlnorm(nCI, 0, 0.3) * 1e5
  A_tech <- matrix(rexp(J * J), J, J); A_tech <- A_tech / colSums(A_tech)
  Z <- matrix(0, nCI, nCI)
  for (b in seq_len(nCI)) {
    cb <- country[b]
    alloc <- attract[country] * A_tech[cbind(rep(seq_len(J), nC),
                                             (b - 1) %% J + 1)] * noise_fixed[, b]
    alloc[country == cb] <- 0
    alloc <- alloc / sum(alloc) * 0.15                   # foreign inputs
    dom <- (country == cb) * A_tech[cbind(rep(seq_len(J), nC), (b - 1) %% J + 1)]
    dom <- dom / sum(dom) * 0.35                         # domestic inputs
    Z[, b] <- (alloc * exp(rnorm(nCI, 0, 0.1)) + dom) * go[b]
  }
  va <- go - colSums(Z)
  saveRDS(list(Z = Z, go = go, va = va, country = country,
               industry = industry), file.path(DIR_PROC, sprintf("wiod_%d.rds", yr)))
}

# ---- BIS-like banking network -------------------------------------------------
rep19 <- c("AU","AT","BE","CL","FR","DE","GR","IE","IT","JP","NL","PT","ES",
           "SE","CH","TW","TR","GB","US")
quarters <- seq(as.Date("2007-01-01"), as.Date("2025-10-01"), by = "quarter")
size <- setNames(c(3, 2, 3, .5, 9, 10, 1, 3, 4, 8, 5, 1, 5, 2, 6, .7, .5, 12, 14),
                 rep19)
eps <- matrix(rnorm(19^2, 0, 1), 19, 19)
cl_list <- list(); tot_list <- list()
for (k in seq_along(quarters)) {
  eps <- 0.97 * eps + sqrt(1 - 0.97^2) * matrix(rnorm(19^2), 19, 19)
  sz <- size * exp(0.01 * k * c(rep(0, 9), 0.6, rep(0, 9)))   # Japan grows
  C <- outer(sz, sz^0.9) * exp(eps); diag(C) <- 0
  C <- C * 100
  tot <- rowSums(C) / runif(19, 0.62, 0.75)
  ij <- which(C > 0, arr.ind = TRUE)
  cl_list[[k]] <- data.table(quarter = quarters[k], creditor = rep19[ij[, 1]],
                             counterparty = rep19[ij[, 2]], claims = C[ij])
  tot_list[[k]] <- data.table(quarter = quarters[k], creditor = rep19,
                              total_claims = tot)
}
saveRDS(rbindlist(cl_list), file.path(DIR_PROC, "bis_claims.rds"))
saveRDS(rbindlist(tot_list), file.path(DIR_PROC, "bis_totals.rds"))

# ---- SRISK-like monthly panel (18 countries: Ireland missing) ---------------
months <- seq(as.Date("2006-01-01"), as.Date("2025-12-01"), by = "month")
sr <- CJ(date = months, country = setdiff(rep19, "IE"))
sr[, srisk := pmax(0, size[country] * 5 * (1 + 2 * (date >= "2007-07-01" &
        date <= "2009-06-30") * runif(.N, 0.5, 1.5)) + rnorm(.N, 0, 5))]
saveRDS(sr, file.path(DIR_PROC, "srisk.rds"))
log_msg("synthetic data written to ", DIR_PROC)
