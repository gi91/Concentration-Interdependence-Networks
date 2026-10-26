# =============================================================================
# fixes.R -- tools that address the problems raised in the review.
# None of them changes a definition of the paper: they add diagnostics,
# alternative nulls and κ-free summaries that can be reported alongside.
# =============================================================================

# ---- (P3) κ-free within-order shares ---------------------------------------
# At a given order h, the one-way/open and reciprocal/circuit components share
# the same normalising constant κ_{N,h}, so share_F(h) = ψ_F(h) / ψ_T(h) does
# not depend on κ (nor on N through κ).  Cross-order shares of ψ_TOT do.
within_order_shares <- function(sd) {
  tb <- sd$table
  hs <- sort(unique(tb$h))
  data.frame(h = hs,
    closed_share = sapply(hs, function(h) {
      pT <- tb$psi[tb$comp == "T" & tb$h == h]
      if (pT > 0) tb$psi[tb$comp == "F" & tb$h == h] / pT else NA }),
    Gamma_open   = sapply(hs, function(h) tb$Gamma[tb$comp == "R" & tb$h == h]),
    Gamma_closed = sapply(hs, function(h) tb$Gamma[tb$comp == "F" & tb$h == h]))
}

# ---- (P8) who drives an alignment factor? ----------------------------------
# Gamma_M(w) = sum over unordered pairs {i,j} of c_ij, with
# c_ij = 2 M_ij w_i w_j / (xbar * sum(M)) for symmetric M; c_ij sums to Gamma.
pair_contributions <- function(w, M, labels = names(w), top = 10) {
  M <- sym_part(M); diag(M) <- 0
  N <- length(w); K <- N * (N - 1)
  xbar <- (1 - sum(w^2)) / K
  ij <- which(upper.tri(M) & M > 0, arr.ind = TRUE)
  if (is.null(labels)) labels <- as.character(seq_len(N))
  c_ij <- 2 * M[ij] * w[ij[, 1]] * w[ij[, 2]] / (xbar * sum(M))
  out <- data.frame(pair = paste(labels[ij[, 1]], labels[ij[, 2]], sep = "-"),
                    i = ij[, 1], j = ij[, 2], contribution = c_ij)
  out <- out[order(-out$contribution), ]
  out$cum_share <- cumsum(out$contribution) / sum(out$contribution)
  attr(out, "Gamma") <- sum(c_ij)
  head(out, top)
}

# ---- (P5) directed binary configuration model (soft degree constraints) ----
# Convention A[to, from].  p[to, from] = x_from y_to / (1 + x_from y_to),
# fitted so that expected out- and in-degrees equal the observed ones
# (Squartini and Garlaschelli, 2011).  Fast to sample: suited to large
# Monte Carlo designs where exact rewiring is too slow.
dbcm_fit <- function(A, tol = 1e-12, max_iter = 50000) {
  A <- (A > 0) * 1; diag(A) <- 0
  kout <- colSums(A); kin <- rowSums(A); N <- nrow(A)
  x <- kout / sqrt(sum(A) + 1); y <- kin / sqrt(sum(A) + 1)
  for (it in seq_len(max_iter)) {
    D <- outer(y, x)                       # [to, from] = y_to x_from
    Q <- 1 / (1 + D); diag(Q) <- 0
    xn <- ifelse(kout > 0, kout / colSums(Q * y), 0)
    D <- outer(y, xn); Q <- 1 / (1 + D); diag(Q) <- 0
    yn <- ifelse(kin > 0, kin / rowSums(sweep(Q, 2, xn, "*")), 0)
    if (max(abs(xn - x), abs(yn - y)) < tol * (1 + max(x, y))) { x <- xn; y <- yn; break }
    x <- xn; y <- yn
  }
  P <- outer(y, x); P <- P / (1 + P)
  # nodes linked to (or from) every other node have p = 1 (x or y infinite)
  P[, kout == N - 1] <- 1; P[kin == N - 1, ] <- 1; diag(P) <- 0
  list(P = P, x = x, y = y, iterations = it,
       max_degree_error = max(abs(colSums(P) - kout), abs(rowSums(P) - kin)))
}
dbcm_sample <- function(fit) {
  N <- nrow(fit$P); A <- (matrix(runif(N * N), N, N) < fit$P) * 1; diag(A) <- 0; A
}

# Test of a statistic stat_fun(A) (weights held fixed) against DBCM draws
dbcm_test <- function(A, stat_fun, B = 499, fit = NULL) {
  if (is.null(fit)) fit <- dbcm_fit(A)
  obs <- stat_fun(A)
  null <- vapply(seq_len(B), function(b) stat_fun(dbcm_sample(fit)), numeric(1))
  list(obs = obs, null = null, p = (1 + sum(null >= obs - 1e-12)) / (B + 1),
       z = (obs - mean(null)) / sd(null))
}

# ---- (P1) node-level embeddedness in excess of its degree-conditional mean --
# e_i(F) observed minus its mean over degree-preserving networks: removes the
# part of reciprocal embeddedness that is implied by how many systems each
# system is exposed to and is exposed by (in/out-degree), i.e. by its
# importance as a counterparty in the binary network.
excess_embeddedness <- function(w, A, B = 500, null = c("rewire", "dbcm"),
                                seed = NULL) {
  null <- match.arg(null)
  if (!is.null(seed)) set.seed(seed)
  eF <- function(X) { F <- X * t(X); embeddedness(w, F) }
  obs <- eF(A)
  fit <- if (null == "dbcm") dbcm_fit(A) else NULL
  sims <- vapply(seq_len(B), function(b) eF(if (null == "dbcm") dbcm_sample(fit)
                                            else rewire_directed(A)),
                 numeric(length(w)))
  m <- rowMeans(sims); s <- apply(sims, 1, sd)
  data.frame(e_obs = obs, e_null_mean = m, e_excess = obs - m,
             e_z = ifelse(s > 0, (obs - m) / s, 0))
}

# ---- (P6) threshold profile of an alignment factor --------------------------
gamma_threshold_profile <- function(w, S, thetas, comp = "F", h = 1) {
  sapply(thetas, function(th) {
    A <- (S >= th) * 1
    if (sum(A) == 0) return(NA)
    tb <- structural_decomposition(w, A, H = h, alpha = rep(1 / h, h))$table
    tb$Gamma[tb$comp == comp & tb$h == h]
  })
}
# threshold that reproduces a target directed density (fixed-density design)
threshold_for_density <- function(S, density) {
  v <- sort(S[row(S) != col(S)], decreasing = TRUE)
  v[max(1, round(density * length(v)))]
}

# ---- (P4) exact gap between Prop. 5 and the two-player Shapley value --------
# nested weight effect - Shapley weight effect = -(1/4) (mbar_t - mbar_{t-1}) I,
# I = Gamma(w_t,M_t) - Gamma(w_t,M_{t-1}) - Gamma(w_{t-1},M_t) + Gamma(w_{t-1},M_{t-1})
shapley_gap <- function(w0, M0, w1, M1) {
  G <- function(w, M) factorize(w, M)$Gamma
  I <- G(w1, M1) - G(w1, M0) - G(w0, M1) + G(w0, M0)
  dm <- factorize(w1, M1)$mbar - factorize(w0, M0)$mbar
  c(interaction = I, delta_mbar = dm, gap_weights = -dm * I / 4)
}
