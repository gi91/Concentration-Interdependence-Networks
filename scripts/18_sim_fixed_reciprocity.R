# =============================================================================
# 18_sim_fixed_reciprocity.R -- experiments E7 and E8 of Table S16
#
# E7  dyad-independent p1 networks (Holland and Leinhardt, 1981): the four dyad
#     states (none, i->j, j->i, reciprocal) have probabilities proportional to
#     1, exp(a_i + b_j), exp(a_j + b_i), exp(a_i + b_j + a_j + b_i + lambda_ij).
#     a_i = c + s z_i, b_i = c + s z'_i, corr(z, z') = 0.9, s = 1 (N = 19) or 1.2
#     (N = 80), c set to give densities 0.35 (N = 19) and 0.07 (N = 80);
#     log-normal weights (sigma = 1.2) correlated with z (rho = 0.8).
#     Cells: random reciprocity (lambda = 0); abundant reciprocity (lambda = 1 for
#     N = 19, 2.5 for N = 80); node-level propensity (abundant lambda plus
#     0.6 (r_i + r_j), r = standardised log weight); preference among the four
#     heaviest nodes (abundant lambda plus 2 for those pairs).  With a common lambda the network is
#     uniform given its degrees and number of reciprocal pairs, so (c) and (d)
#     are exact nulls.
# E8  designs of Table S3 (09_sim_design.R) with N = 19, rho = 0.8, delta = 0,
#     core-periphery and Chung-Lu; networks with at least one reciprocal pair
#     and no isolated node.
# Tests of psi_F^(1): fixed pairs (S1(c)), fixed partners (S1(d)), exact rewiring
# (S1(a)), each with 199 draws of one chain (500|E| burn-in, 100|E| thinning).
# Replications: NCI_REPS_E7 (default 500) per cell, NCI_REPS_E8 (400) per topology.
# Output: simE7_fixed_reciprocity.csv, simE8_fixed_reciprocity.csv, and rows
#         E7-E8 appended to tabS2_robustness_summary (run after 08).
# =============================================================================
source("R/utils.R"); source("R/core_indices.R"); source("R/fixes.R")
source("R/constrained_nulls.R")
cn_compile()
REPS7 <- as.integer(Sys.getenv("NCI_REPS_E7", 500))
REPS8 <- as.integer(Sys.getenv("NCI_REPS_E8", 400))
BD    <- 199
CORES <- as.integer(Sys.getenv("NCI_CORES", 1))
par_lapply <- function(X, FUN) if (CORES > 1 && .Platform$OS.type == "unix")
  parallel::mclapply(X, FUN, mc.cores = CORES, mc.preschedule = FALSE) else lapply(X, FUN)
ptest <- function(A, w, null, seed)
  constrained_test(A, w, null, n_chains = 1, n_draws = BD, seed = seed, cores = 1)$p

# ---- E7: p1 networks ---------------------------------------------------------------
p1_probs <- function(a, b, Lam) {           # matrices of state probabilities, i < j
  e_ij <- exp(outer(a, b, "+")); e_ji <- t(e_ij)
  e_r <- e_ij * e_ji * exp(Lam); Z <- 1 + e_ij + e_ji + e_r
  list(ij = e_ij / Z, ji = e_ji / Z, r = e_r / Z)
}
p1_density <- function(cc, z1, z2, s, Lam) {
  N <- length(z1); P <- p1_probs(cc + s * z1, cc + s * z2, Lam); up <- upper.tri(Lam)
  sum((P$ij + P$ji + 2 * P$r)[up]) / (N * (N - 1))
}
draw_p1 <- function(N, cell) {
  s <- if (N == 19) 1 else 1.2; dens <- if (N == 19) 0.35 else 0.07
  lam_ab <- if (N == 19) 1 else 2.5
  z1 <- rnorm(N); z2 <- 0.9 * z1 + sqrt(1 - 0.81) * rnorm(N)
  zw <- 0.8 * z1 + 0.6 * rnorm(N)
  wv <- sort(rlnorm(N, 0, 1.2)); w <- wv[rank(zw)]; w <- w / sum(w)
  r <- as.numeric(scale(log(w)))
  Lam <- switch(cell,
    random = matrix(0, N, N), abundant = matrix(lam_ab, N, N),
    propensity = lam_ab + 0.6 * outer(r, r, "+"),
    preference = { top <- order(w, decreasing = TRUE)[1:4]; L <- matrix(lam_ab, N, N)
                   L[top, top] <- lam_ab + 2; L })
  cc <- uniroot(function(x) p1_density(x, z1, z2, s, Lam) - dens, c(-15, 10))$root
  P <- p1_probs(cc + s * z1, cc + s * z2, Lam)
  A <- matrix(0, N, N); up <- which(upper.tri(A), arr.ind = TRUE)
  u <- runif(nrow(up))
  pij <- P$ij[up]; pji <- P$ji[up]; pr <- P$r[up]
  st <- ifelse(u < pij, 1, ifelse(u < pij + pji, 2, ifelse(u < pij + pji + pr, 3, 0)))
  A[up[st %in% c(1, 3), , drop = FALSE]] <- 1
  A[up[st %in% c(2, 3), 2:1, drop = FALSE]] <- 1
  list(A = A, w = w)
}
TABLE_ONLY <- Sys.getenv("NCI_E78_TABLE_ONLY") == "1"   # redraw Table S16 rows from the saved csv
cells7 <- expand.grid(cell = c("random", "abundant", "propensity", "preference"),
                      N = c(19, 80), stringsAsFactors = FALSE)
if (TABLE_ONLY) {
  s7 <- read.csv(file.path(DIR_TAB, "simE7_fixed_reciprocity.csv"))
  s8 <- read.csv(file.path(DIR_TAB, "simE8_fixed_reciprocity.csv"))
} else {
e7 <- do.call(rbind, lapply(seq_len(nrow(cells7)), function(g) {
  cl <- cells7$cell[g]; N <- cells7$N[g]
  r <- par_lapply(seq_len(REPS7), function(k) {
    set.seed(700000 + 1000 * g + k)
    repeat { x <- draw_p1(N, cl); if (sum(x$A * t(x$A)) > 0) break }
    data.frame(cell = cl, N = N, density = sum(x$A) / (N * (N - 1)),
               reciprocity = sum(x$A * t(x$A)) / sum(x$A),
               p_pairs = ptest(x$A, x$w, "pairs", 7 * k + g),
               p_partners = ptest(x$A, x$w, "partners", 7 * k + g + 1),
               p_rewire = if (cl %in% c("random", "abundant")) ptest(x$A, x$w, "rewire", 7 * k + g + 2) else NA)
  })
  log_msg("E7 ", cl, " N = ", N)
  do.call(rbind, r)
}))
s7 <- aggregate(cbind(density, reciprocity, rej_pairs = p_pairs <= .05,
                      rej_partners = p_partners <= .05, rej_rewire = p_rewire <= .05) ~ cell + N,
                data = e7, FUN = mean, na.action = na.pass)
write.csv(s7, file.path(DIR_TAB, "simE7_fixed_reciprocity.csv"), row.names = FALSE)
print(s7, digits = 3, row.names = FALSE)

# ---- E8: designs of Table S3 -------------------------------------------------------
KBAR <- 6
gen_network <- function(topo, N, f) {            # as in 09_sim_design.R
  p0 <- KBAR / (N - 1)
  P <- switch(topo,
    "Chung-Lu" = { P <- outer(f, f); P * (p0 * N * (N - 1) / sum(P[row(P) != col(P)])) },
    "Core-periphery" = { core <- f >= quantile(f, .8)
      P <- outer(core, core, function(a, b) ifelse(a & b, 8, ifelse(a | b, 3, 1)))
      P * (p0 * N * (N - 1) / sum(P[row(P) != col(P)])) })
  P <- pmin(P, 1); diag(P) <- 0
  A <- (matrix(runif(N * N), N, N) < P) * 1; diag(A) <- 0; A
}
draw_case <- function(topo, N = 19, rho = 0.8) {
  zf <- rnorm(N); f <- exp(zf)
  zw <- rho * zf + sqrt(1 - rho^2) * rnorm(N)
  wv <- sort(rlnorm(N, 0, 1.2)); w <- wv[rank(zw)]; w <- w / sum(w)
  list(A = gen_network(topo, N, f), w = w)
}
e8 <- do.call(rbind, lapply(c("Core-periphery", "Chung-Lu"), function(topo) {
  r <- par_lapply(seq_len(REPS8), function(k) {
    set.seed(800000 + 1000 * (topo == "Chung-Lu") + k)
    repeat { x <- draw_case(topo)
      if (sum(x$A * t(x$A)) > 0 && all(rowSums(x$A) + colSums(x$A) > 0)) break }
    data.frame(topo = topo, p_pairs = ptest(x$A, x$w, "pairs", 9 * k),
               p_partners = ptest(x$A, x$w, "partners", 9 * k + 1),
               p_rewire = ptest(x$A, x$w, "rewire", 9 * k + 2))
  })
  log_msg("E8 ", topo)
  do.call(rbind, r)
}))
s8 <- aggregate(cbind(rej_pairs = p_pairs <= .05, rej_partners = p_partners <= .05,
                      rej_rewire = p_rewire <= .05) ~ topo, data = e8, FUN = mean)
write.csv(s8, file.path(DIR_TAB, "simE8_fixed_reciprocity.csv"), row.names = FALSE)
print(s8, digits = 3, row.names = FALSE)
}

# ---- rows E7-E8 of Table S16 -------------------------------------------------------
g7 <- function(cl, N, v) s7[s7$cell == cl & s7$N == N, v]
f3 <- function(x) sprintf("%.3f", x)
rng <- function(x) paste(f3(min(x)), f3(max(x)), sep = "--")
rows <- rbind(
  c("E7", "Fixed reciprocal pairs, size at 5\\%, $N=19$ / $N=80$ (random or abundant reciprocity)",
    paste(rng(c(g7("random", 19, "rej_pairs"), g7("abundant", 19, "rej_pairs"))),
          rng(c(g7("random", 80, "rej_pairs"), g7("abundant", 80, "rej_pairs"))), sep = " / ")),
  c("E7", "Exact rewiring on $\\psi_F^{(1)}$, rejection with abundant reciprocity, $N=19$ / $N=80$",
    paste(f3(g7("abundant", 19, "rej_rewire")), f3(g7("abundant", 80, "rej_rewire")), sep = " / ")),
  c("E7", "Node-level propensity: fixed pairs / fixed partners, $N=19$; $N=80$",
    sprintf("%s / %s; %s / %s", f3(g7("propensity", 19, "rej_pairs")), f3(g7("propensity", 19, "rej_partners")),
            f3(g7("propensity", 80, "rej_pairs")), f3(g7("propensity", 80, "rej_partners")))),
  c("E7", "Preference among the 4 heaviest: fixed pairs / fixed partners, $N=19$; $N=80$",
    sprintf("%s / %s; %s / %s", f3(g7("preference", 19, "rej_pairs")), f3(g7("preference", 19, "rej_partners")),
            f3(g7("preference", 80, "rej_pairs")), f3(g7("preference", 80, "rej_partners")))),
  c("E8", "Fixed pairs / fixed partners / exact rewiring, size at 5\\%: core--periphery; Chung--Lu ($N=19$)",
    sprintf("%s / %s / %s; %s / %s / %s",
            f3(s8$rej_pairs[s8$topo == "Core-periphery"]), f3(s8$rej_partners[s8$topo == "Core-periphery"]),
            f3(s8$rej_rewire[s8$topo == "Core-periphery"]), f3(s8$rej_pairs[s8$topo == "Chung-Lu"]),
            f3(s8$rej_partners[s8$topo == "Chung-Lu"]), f3(s8$rej_rewire[s8$topo == "Chung-Lu"]))))
rows <- as.data.frame(rows); names(rows) <- c("Experiment", "Quantity", "Value")
write.csv(rows, file.path(DIR_TAB, "tabS2_rows_E7_E8.csv"), row.names = FALSE)
f08 <- file.path(DIR_TAB, "tabS2_robustness_summary.csv")
if (file.exists(f08)) {
  old <- read.csv(f08, check.names = FALSE); names(rows) <- names(old)
  tabS2 <- rbind(old, rows)
  write_tex_table(tabS2, "tabS2_robustness_summary.tex",
    caption = "Monte Carlo robustness experiments: summary", label = "tab:sim-robust",
    align = "llc", header = names(tabS2),
    notes = "E1--E6: scripts/08\\_sim\\_robustness.R; E7--E8: scripts/18\\_sim\\_fixed\\_reciprocity.R.")
}
print(rows, row.names = FALSE)
log_msg("E7-E8 written")
