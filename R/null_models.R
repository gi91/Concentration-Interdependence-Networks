# =============================================================================
# null_models.R -- inference for the alignment factors (Proposition 1)
#   (1) permutation of the weights across nodes  -> E(Gamma) = 1
#   (2) analytic permutation moments (Mantel / QAP, Hubert 1987)
#   (3) directed degree-preserving rewiring (in- and out-degree fixed)
# =============================================================================

# ---- (1) permutation test ----------------------------------------------------
# p_share: share of permutations with Gamma >= observed (definition used in
#          Table 4 of the paper);  p_plus1: (1 + #) / (B + 1), never zero.
perm_alignment <- function(w, M, B = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  M <- sym_part(M); diag(M) <- 0
  if (sum(M) == 0) return(list(Gamma = NA, null = NA, p_share = NA,
                                p_plus1 = NA))
  W <- tcrossprod(w); diag(W) <- 0
  den <- sum(W) / length(w) / (length(w) - 1) * sum(M)
  G_obs <- sum(M * W) / den
  G_null <- vapply(seq_len(B), function(b) {
    wp <- w[sample.int(length(w))]
    Wp <- tcrossprod(wp); diag(Wp) <- 0
    sum(M * Wp) / den
  }, numeric(1))
  cnt <- sum(G_null >= G_obs - 1e-12)
  list(Gamma = G_obs, null = G_null, p_share = cnt / B,
       p_plus1 = (1 + cnt) / (B + 1))
}

# ---- (2) analytic permutation moments of Gamma -------------------------------
# Gamma_M(w) is, up to a constant, the Mantel/QAP statistic
#   sum_{i != j} a_ij b_ij,  a = sym(M), b = w w' (zero diagonals),
# under a uniformly random relabelling of b.  Exact mean and variance:
#   E  = A1 B1 / n(n-1)
#   E2 = 2 A2 B2 / n(n-1) + 4 A3 B3 / n(n-1)(n-2) + A4 B4 / n(n-1)(n-2)(n-3)
# with A1 = sum a, A2 = sum a^2, A3 = sum_i (a_i+)^2 - A2,
#      A4 = A1^2 - 2 A2 - 4 A3   (sums over ordered off-diagonal pairs).
mantel_moments <- function(w, M) {
  a <- sym_part(M); diag(a) <- 0
  b <- tcrossprod(w); diag(b) <- 0
  n <- length(w)
  mom <- function(x) {
    x1 <- sum(x); x2 <- sum(x^2); x3 <- sum(rowSums(x)^2) - x2
    c(x1, x2, x3, x1^2 - 2 * x2 - 4 * x3)
  }
  A <- mom(a); Bm <- mom(b)
  n2 <- n * (n - 1); n3 <- n2 * (n - 2); n4 <- n3 * (n - 3)
  E  <- A[1] * Bm[1] / n2
  E2 <- 2 * A[2] * Bm[2] / n2 + 4 * A[3] * Bm[3] / n3 + A[4] * Bm[4] / n4
  scale <- (Bm[1] / n2) * A[1]          # Gamma = statistic / scale
  obs <- sum(a * b) / scale
  sdG <- sqrt(max(E2 - E^2, 0)) / scale
  list(Gamma = obs, mean = E / scale, sd = sdG,
       z = (obs - E / scale) / sdG)
}

# ---- (3) directed degree-preserving rewiring --------------------------------
# A: binary, A[i, j] = 1 for a link j -> i.  Double-edge swaps
# (a->b, c->d) -> (a->d, c->b) preserve every in- and out-degree; they are
# combined with reversals of directed 3-cycles, without which the chain is
# not ergodic on directed graphs (Rao, Jana and Bandyopadhyay, 1996).
rewire_directed <- function(A, n_swap = 10 * sum(A), p_triangle = 0.05) {
  A <- (A > 0) * 1; diag(A) <- 0
  N <- nrow(A)
  E <- which(A > 0, arr.ind = TRUE)          # columns: row = to, col = from
  to <- E[, 1]; from <- E[, 2]; m <- length(to)
  if (m < 2) return(A)
  for (s in seq_len(n_swap)) {
    if (runif(1) < p_triangle) {
      e <- sample.int(m, 1); a <- from[e]; b <- to[e]
      c3 <- sample.int(N, 1)
      if (c3 != a && c3 != b && A[c3, b] == 1 && A[a, c3] == 1 &&
          A[a, b] == 0 && A[b, c3] == 0 && A[c3, a] == 0) {
        # cycle a->b->c->a  becomes  a->c->b->a
        A[b, a] <- 0; A[c3, b] <- 0; A[a, c3] <- 0
        A[a, b] <- 1; A[b, c3] <- 1; A[c3, a] <- 1
        E <- which(A > 0, arr.ind = TRUE); to <- E[, 1]; from <- E[, 2]
      }
      next
    }
    e <- sample.int(m, 2)
    a <- from[e[1]]; b <- to[e[1]]; c <- from[e[2]]; d <- to[e[2]]
    if (a == d || c == b || a == c || b == d) next
    if (A[d, a] == 1 || A[b, c] == 1) next
    A[b, a] <- 0; A[d, c] <- 0; A[d, a] <- 1; A[b, c] <- 1
    to[e[1]] <- d; to[e[2]] <- b
  }
  A
}

# stat_fun(A) must return a named numeric vector.
rewiring_test <- function(A, stat_fun, B = 1000, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  obs <- stat_fun(A)
  null <- matrix(vapply(seq_len(B), function(b) stat_fun(rewire_directed(A, ...)),
                        numeric(length(obs))), ncol = length(obs), byrow = TRUE)
  colnames(null) <- names(obs)
  cnt <- colSums(sweep(null, 2, obs - 1e-12, ">="))
  z <- (obs - colMeans(null)) / apply(null, 2, sd)
  # p_share: definition of Table 4 in the paper; p_plus1 = (1 + #)/(B + 1) is
  # exact under the null and should be used with small B
  list(obs = obs, null = null, p_share = cnt / B, p_plus1 = (1 + cnt) / (B + 1),
       z = z)
}

# Statistics used in Table 4: first-order reciprocal component, its
# availability, and the circuit component of order h = 2 (length-3 circuits)
structural_stats <- function(w) {
  function(A) {
    sd <- structural_decomposition(w, A, H = 2, alpha = c(0.5, 0.5))
    tb <- sd$table
    c(psi_F1  = tb$psi[tb$comp == "F" & tb$h == 1],
      mbar_F1 = tb$mbar[tb$comp == "F" & tb$h == 1],
      psi_R1  = tb$psi[tb$comp == "R" & tb$h == 1],
      psi_F2  = tb$psi[tb$comp == "F" & tb$h == 2])
  }
}
