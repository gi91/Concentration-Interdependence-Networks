# =============================================================================
# 07_sim_propositions.R -- numerical verification of Propositions 1-5 and of
# the auxiliary claims of the paper on random networks.
# Output: output/tables/tabS1_proposition_checks.tex (+ .csv)
# Every exact identity must hold to machine precision; stochastic statements
# (expectations) are checked with Monte Carlo t-statistics.
# =============================================================================
source("R/utils.R"); source("R/core_indices.R"); source("R/null_models.R")
source("R/fixes.R")
set.seed(42)
REPS <- as.integer(Sys.getenv("NCI_REPS", 200))

rdir <- function(N, p, weighted = FALSE) {
  A <- matrix(runif(N * N) < p, N, N) * if (weighted) matrix(runif(N * N), N, N) else 1
  diag(A) <- 0; A
}
rw <- function(N, a) { g <- rgamma(N, a); g / sum(g) }     # Dirichlet(a)
design <- expand.grid(N = c(8, 19, 43), p = c(0.1, 0.3), a = c(0.2, 1, 5),
                      weighted = c(FALSE, TRUE))
err <- function(x) max(abs(x))
checks <- list()
add <- function(id, claim, value, tol, type = "max abs error") {
  checks[[length(checks) + 1]] <<- data.frame(id = id, claim = claim, type = type,
    value = value, pass = if (type == "max abs error") value < tol else abs(value) < tol)
}

# ---- draw a common set of random cases --------------------------------------------
cases <- lapply(seq_len(REPS), function(r) {
  d <- design[(r - 1) %% nrow(design) + 1, ]
  list(A = rdir(d$N, d$p, d$weighted), w = rw(d$N, d$a), d = d)
})

# ---- Proposition 1 ----------------------------------------------------------------
e1 <- sapply(cases, function(cs) { f <- factorize(cs$w, cs$A)
  if (is.na(f$Gamma)) 0 else f$psi - f$mbar * f$Gamma })
add("P1", "Psi = mbar * Gamma (exact factorization)", err(e1), 1e-12)
e1b <- sapply(cases, function(cs) { N <- length(cs$w)
  g <- alignment(rep(1 / N, N), cs$A); if (is.na(g)) 0 else g - 1 })
add("P1", "Gamma = 1 under uniform weights", err(e1b), 1e-12)
e1c <- sapply(cases, function(cs) { mm <- mantel_moments(cs$w, cs$A)
  if (!is.finite(mm$mean)) 0 else mm$mean - 1 })
add("P1", "E_pi Gamma = 1 under random permutation (exact moment)", err(e1c), 1e-12)
# Monte Carlo check of the permutation mean and of the analytic variance
mc <- t(sapply(cases[1:min(40, REPS)], function(cs) {
  pa <- perm_alignment(cs$w, cs$A, B = 4000); mm <- mantel_moments(cs$w, cs$A)
  if (is.na(pa$Gamma)) return(c(0, 0))
  c(t_mean = (mean(pa$null) - 1) / (sd(pa$null) / sqrt(4000)),
    rel_sd = sd(pa$null) / mm$sd - 1) }))
add("P1", "E_pi Gamma = 1 (Monte Carlo, max |t|, 4000 perms)", max(abs(mc[, 1])),
    4, "max |t|")
add("P1", "Analytic permutation sd of Gamma vs Monte Carlo (max rel. error)",
    max(abs(mc[, 2])), 0.1, "max rel. error")
# Erdos-Renyi: E(Psi | w) = p for a symmetric ER network independent of w
er_t <- sapply(1:20, function(r) {
  N <- sample(c(8, 19, 43), 1); p <- runif(1, .05, .5); w <- rw(N, .5)
  ps <- replicate(2000, { U <- matrix(runif(N * N) < p, N, N) * 1
    S <- U * upper.tri(U); S <- S + t(S); psi_index(w, S) })
  (mean(ps) - p) / (sd(ps) / sqrt(2000)) })
add("P1", "E(Psi|w) = p in Erdos-Renyi networks (max |t|, 20 x 2000 draws)",
    max(abs(er_t)), 4, "max |t|")

# ---- Proposition 2 ----------------------------------------------------------------
e2 <- t(sapply(cases, function(cs) {
  tp <- topological_profile(cs$w, cs$A, "binary")
  N <- length(cs$w); tu <- topological_profile(rep(1 / N, N), cs$A, "binary")
  c(sum(tp$profile$psi) + tp$psi_inf - 1, sum(tp$profile$mbar) + tp$mbar_inf - 1,
    tu$dw - tu$dbar, if (any(tp$profile$psi < 0)) 1 else 0) }))
add("P2", "sum_h psi_G^(h) (+ disconnected share) = 1", err(e2[, 1]), 1e-12)
add("P2", "sum_h mbar_h (+ disconnected share) = 1", err(e2[, 2]), 1e-12)
add("P2", "uniform weights: economic mean distance = average path length",
    err(e2[, 3]), 1e-12)
add("P2", "psi_G^(h) >= 0", err(e2[, 4]), 1e-12)
e2b <- sapply(cases, function(cs) { tp <- topological_profile(cs$w, cs$A, "strength")
  if (any(tp$profile$psi < -1e-15 | tp$profile$psi > 1 + 1e-15)) 1 else 0 })
add("Sec. 2", "strength-adjusted psi_G^(h) in [0, 1]", err(e2b), 1e-12)
e2c <- sapply(cases, function(cs) {
  la <- c(.5, .3, .2); tp <- topological_profile(cs$w, cs$A, "binary")
  H <- min(3, nrow(tp$profile)); lam <- la[1:H] / sum(la[1:H])
  aggregate_topological(cs$w, cs$A, lam) - sum(lam * tp$profile$psi[1:H]) })
add("Sec. 2", "aggregate index = weighted average of the profile", err(e2c), 1e-12)

# ---- Proposition 3 ----------------------------------------------------------------
pcases <- lapply(cases, function(cs) { A <- cs$A
  for (j in which(colSums(A) == 0)) A[sample(setdiff(seq_len(nrow(A)), j), 1), j] <- 1
  list(A = A, w = cs$w) })
e3 <- t(sapply(pcases, function(cs) {
  out <- sapply(1:4, function(k) { pm <- propagation_measures(cs$w, cs$A, k)
    N <- length(cs$w)
    c(pm$identity_error, as.numeric(pm$psiO < -1e-15 | pm$psiO > 1 + 1e-15),
      as.numeric(any(pm$Y < 1 / N - 1e-12 | pm$Y > 1 + 1e-12))) })
  apply(abs(out), 1, max) }))
add("P3(i)", "HHI(P^k w) = own + (1 - HHI) psi_O^(k), k = 1..4", err(e3[, 1]), 1e-12)
add("P3(i)", "0 <= psi_O^(k) <= 1", err(e3[, 2]), 1e-12)
add("P3(i)", "1/N <= Y_j^(k) <= 1", err(e3[, 3]), 1e-12)
e3b <- sapply(1:30, function(r) {          # ergodic: dense positive matrices
  N <- sample(c(8, 19), 1); A <- matrix(runif(N * N), N, N); diag(A) <- 0
  w <- rw(N, .5); st <- stationary_dist(propagation_matrix(A)); k <- 60
  pm <- propagation_measures(w, A, k)
  max(abs(c(pm$psiO, pm$CP, pm$hhi_wk, pm$Y) - st$hhi_pi)) })
add("P3(ii)", "ergodic limit: Y, C_P, psi_O, HHI(w^(k)) -> sum pi^2 (k = 60)",
    err(e3b), 1e-8)

# ---- Proposition 4 ----------------------------------------------------------------
e4 <- t(sapply(cases, function(cs) {
  A <- pmin(cs$A, 1); N <- nrow(A)
  sd <- structural_decomposition(cs$w, A, H = 3)
  inrange <- all(sapply(sd$mats, function(M) all(M >= -1e-15 & M <= 1 + 1e-15)))
  cf <- if (N <= 19) max(sapply(2:3, function(h)
    max(abs(path_matrix(A, h, "closed") - path_matrix(A, h, "dfs"))))) else 0
  # kappa invariance: unnormalised path matrices give the same Gamma
  ginv <- max(sapply(2:3, function(h) {
    Tb <- path_matrix(A, h); g1 <- alignment(cs$w, sym_part(t(A) * Tb))
    g2 <- alignment(cs$w, sym_part(t(A) * Tb / kappa_Nh(N, h)))
    if (is.na(g1)) 0 else g1 - g2 }))
  c(as.numeric(!inrange), sd$residual, cf, ginv) }))
add("P4(i)", "entries of M_T, M_R, M_F in [0, 1]", err(e4[, 1]), 1e-12)
add("P4(ii)", "psi_TOT = one-way + reciprocal + open chains + circuits", err(e4[, 2]), 1e-12)
add("P4(iii)", "Gamma invariant to kappa_{N,h}", err(e4[, 4]), 1e-10)
add("Sec. 4", "closed-form T^(2), T^(3) = DFS enumeration of simple paths",
    err(e4[, 3]), 1e-10)
e4b <- sapply(c(5, 8, 10), function(N) { A <- matrix(1, N, N); diag(A) <- 0
  max(sapply(1:min(4, N - 1), function(h) {
    Tb <- if (h <= 3) path_matrix(A, h) else simple_path_matrix(A, h)
    max(abs(Tb[row(Tb) != col(Tb)] - kappa_Nh(N, h))) })) })
add("Sec. 4", "kappa_{N,h} attained in the complete network (h <= 4)", err(e4b), 1e-9)
# circuits are invariant to the size normalisation Z_ji / (X_i^a X_j^(1-a))
e4c <- sapply(1:30, function(r) {
  N <- 12; Z <- matrix(rexp(N * N), N, N) * (runif(N * N) < .5); diag(Z) <- 0
  X <- rowSums(Z) + colSums(Z) + rexp(N) * 10; w <- rw(N, 1)
  g <- sapply(c(0, .25, .5, .75, 1), function(a) {
    M <- t(Z) / outer(X^a, X^(1 - a)); M <- M / max(M)
    tb <- structural_decomposition(w, M, H = 3)$table
    tb$Gamma[tb$comp == "F"] })
  max(abs(g - g[, 1]), na.rm = TRUE) })
add("Sec. 7.1", "Gamma of reciprocal links and circuits invariant to normalisation exponent a",
    err(e4c), 1e-10)

# ---- Proposition 5 and the intertemporal index ------------------------------------
e5 <- t(sapply(cases, function(cs) {
  N <- length(cs$w); w1 <- rw(N, cs$d$a)
  A1 <- cs$A; flip <- matrix(runif(N * N) < .2, N, N); A1[flip] <- 1 - (A1[flip] > 0)
  diag(A1) <- 0
  if (sum(cs$A) == 0 || sum(A1) == 0) return(c(0, 0, 0))
  d01 <- dynamic_decomposition(cs$w, cs$A, w1, A1)
  d10 <- dynamic_decomposition(w1, A1, cs$w, cs$A)
  c(d01$residual,
    d01$availability + d10$availability, d01$network + d10$network) }))
add("P5", "Delta psi = availability + weights + network", err(e5[, 1]), 1e-12)
add("P5", "time reversal flips every component (symmetry)", err(e5[, 2:3]), 1e-12)
e6 <- t(sapply(cases, function(cs) {
  N <- length(cs$w); G <- pmin(cs$A, 1); wl <- rw(N, 1)
  it <- psi_intertemporal(cs$w, wl, G)
  same <- psi_intertemporal(cs$w, cs$w, sym_part(G))$psi - psi_index(cs$w, sym_part(G))
  e <- embeddedness(cs$w, G)
  dec <- sum(cs$w * (1 - cs$w) * e) / (1 - sum(cs$w^2)) - psi_index(cs$w, G)
  c(as.numeric(it$psi < -1e-15 | it$psi > 1 + 1e-15), same, dec) }))
add("Sec. 5", "intertemporal index in [0, 1]", err(e6[, 1]), 1e-12)
add("Sec. 5", "intertemporal index = quadratic index when w_t = w_{t-h}", err(e6[, 2]), 1e-12)
add("Sec. 7.3", "node decomposition: Psi = sum_i w_i (1-w_i) e_i / (1-HHI)", err(e6[, 3]), 1e-12)

# ---- identity used in the critique: psi_T^(1) = (psi_D + psi_F^(1)) / 2 ---------
e7 <- sapply(cases, function(cs) { A <- (cs$A > 0) * 1
  tb <- structural_decomposition(cs$w, A, H = 1, alpha = 1)$table
  tb$psi[tb$comp == "T"] - (baseline_index(cs$w, A)$psiD + tb$psi[tb$comp == "F"]) / 2 })
add("Remark", "binary: psi_T^(1) = (psi_D + psi_F^(1)) / 2, so psi_R^(1) = (psi_D - psi_F^(1)) / 2",
    err(e7), 1e-12)

# ---- tools introduced to address the review (R/fixes.R) --------------------------
e8 <- t(sapply(cases, function(cs) {
  N <- length(cs$w); w1 <- rw(N, cs$d$a); A <- (cs$A > 0) * 1
  A1 <- A; fl <- matrix(runif(N * N) < .2, N, N); A1[fl] <- 1 - A1[fl]; diag(A1) <- 0
  if (sum(A) == 0 || sum(A1) == 0) return(c(0, 0, 0, 0))
  d <- dynamic_decomposition(cs$w, A, w1, A1); g <- shapley_gap(cs$w, A, w1, A1)
  F <- A * t(A)
  pc <- if (sum(F) > 0) sum(pair_contributions(cs$w, F, top = Inf)$contribution) -
    alignment(cs$w, F) else 0
  sd <- structural_decomposition(cs$w, A, H = 3)
  wo <- within_order_shares(sd)
  # within-order share computed without kappa (raw path counts)
  raw <- sapply(2:3, function(h) { Tb <- path_matrix(A, h)
    sT <- psi_index(cs$w, sym_part(Tb)); if (sT == 0) return(NA)
    psi_index(cs$w, sym_part(t(A) * Tb)) / sT })
  c(d$weights - d$shapley_weights - g[["gap_weights"]], pc,
    { dd <- abs(wo$closed_share[2:3] - raw); if (all(is.na(dd))) 0 else max(dd, na.rm = TRUE) },
    0) }))
add("Fix P4", "nested - Shapley weight effect = -(1/4) Delta mbar x interaction", err(e8[, 1]), 1e-12)
add("Fix P8", "pair contributions sum to Gamma", err(e8[, 2]), 1e-12)
add("Fix P3", "within-order closed share independent of kappa", err(e8[, 3]), 1e-12)
# (N >= 19: in very small sparse graphs the ML solution can lie on the boundary,
#  e.g. a node linked to every node that has positive in-degree)
e9 <- sapply(Filter(function(cs) cs$d$N >= 19, cases), function(cs) {
  A <- (cs$A > 0) * 1; dbcm_fit(A)$max_degree_error })
add("Fix P5", "DBCM reproduces expected in- and out-degrees (N >= 19)", err(e9), 1e-6)

res <- do.call(rbind, checks)
print(res, right = FALSE)
# LaTeX labels for the paper, in the order of the checks above, with the
# numbering of the paper (Methodology 2.1-2.7, applications 5.x)
tex_id <- c(rep("Prop.~1", 4), "Prop.~6", "Prop.~1", rep("Prop.~2", 4),
            "Sec.~2.2", "Sec.~2.2", rep("Prop.~3(i)", 3), "Prop.~3(ii)",
            "Prop.~4(i)", "Prop.~4(ii)", "Prop.~4(iii)", "Sec.~2.4", "Sec.~2.4",
            "Sec.~5.1", "Prop.~5", "Prop.~5", "Sec.~2.5", "Sec.~2.5", "Cor.~1",
            "Rem.~1", "Cor.~3", "Cor.~1", "Cor.~2", "Def.~1")
tex_claim <- c(
  "$\\Psi=\\bar m_M\\Gamma_M$ (exact factorization)",
  "$\\Gamma_M=1$ under uniform weights",
  "$E_\\pi\\Gamma_M=1$ under random permutation (exact moment)",
  "$E_\\pi\\Gamma_M=1$ (Monte Carlo, max $|t|$, 4000 permutations)",
  "Analytic permutation s.d.\\ of $\\Gamma_M$ vs Monte Carlo (max rel.\\ error)",
  "$E(\\Psi\\mid w)=p$ in Erd\\H{o}s--R\\'enyi networks (max $|t|$, $20\\times2000$ draws)",
  "$\\sum_h\\psi^{(h)}_G=1$ (with disconnected share)",
  "$\\sum_h\\bar m_h=1$ (with disconnected share)",
  "Uniform weights: $\\bar d_w=\\bar d$",
  "$\\psi^{(h)}_G\\ge0$",
  "Strength-adjusted $\\psi^{(h)}_G\\in[0,1]$",
  "Aggregate index $=$ weighted average of the profile",
  "$\\mathrm{HHI}(P^kw)=\\sum_jw_j^2Y_j^{(k)}+(1-\\mathrm{HHI})\\psi^{(k)}_O$, $k=1,\\dots,4$",
  "$0\\le\\psi^{(k)}_O\\le1$",
  "$1/N\\le Y^{(k)}_j\\le1$",
  "Ergodic limit of $Y$, $C_P$, $\\psi_O$, $\\mathrm{HHI}(w^{(k)})$ ($k=60$)",
  "Entries of $M_T$, $M_R$, $M_F$ in $[0,1]$",
  "$\\psi_{\\mathrm{TOT}}=$ one-way $+$ reciprocal $+$ open chains $+$ circuits",
  "$\\Gamma$ invariant to $\\kappa_{N,h}$",
  "Closed forms of $\\bar T^{(2)}$, $\\bar T^{(3)}$ $=$ enumeration of simple paths",
  "$\\kappa_{N,h}$ attained in the complete network ($h\\le4$)",
  "$\\Gamma$ of reciprocal links and circuits invariant to the normalization exponent $a$",
  "$\\Delta\\psi=$ availability $+$ weights $+$ network",
  "Time reversal flips every component",
  "Intertemporal index in $[0,1]$",
  "Intertemporal index $=$ quadratic index when $w_t=w_{t-h}$",
  "$\\Psi=\\sum_iw_i(1-w_i)e_i/(1-\\mathrm{HHI})$",
  "Binary: $\\psi^{(1)}_R=(\\psi_D-\\psi^{(1)}_F)/2$",
  "Nested $-$ Shapley weight effect $=-\\frac14\\Delta\\bar m\\, I_t$",
  "Pair contributions sum to $\\Gamma$",
  "Within-order closed share independent of $\\kappa_{N,h}$",
  "DBCM reproduces expected in- and out-degrees ($N\\ge19$)")
stopifnot(length(tex_id) == nrow(res), length(tex_claim) == nrow(res))
res$id <- tex_id; res$claim <- tex_claim
tabS1 <- data.frame(Result = res$id, Claim = res$claim, Statistic = res$type,
                    Value = formatC(res$value, format = "e", digits = 1),
                    Pass = ifelse(res$pass, "yes", "\\textbf{NO}"))
write_tex_table(tabS1, "tabS1_proposition_checks.tex",
  caption = "Numerical verification of the propositions on random networks",
  label = "tab:sim-props", align = "llllc",
  header = c("Result", "Claim", "Statistic", "Value", "Pass"),
  notes = paste0(REPS, " random directed networks with $N\\in\\{8,19,43\\}$, density",
    " $\\{0.1,0.3\\}$, Dirichlet weights with concentration $\\{0.2,1,5\\}$,",
    " binary and weighted links. Exact identities: maximum absolute error.",
    " Expectations: maximum $|t|$ of Monte Carlo means."))
log_msg("Proposition checks: ", sum(res$pass), "/", nrow(res), " passed")
