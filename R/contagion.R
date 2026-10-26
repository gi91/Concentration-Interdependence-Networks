# =============================================================================
# contagion.R -- Section 7.4: linear DebtRank and reciprocity-tilted MCMC
#
# Here X[i, j] = 1 if system i is EXPOSED TO (holds claims on) j.  In the
# structural convention of Section 4 a link j -> i ("losses travel from the
# debtor j to the creditor i") is A[i, j] = 1, so X and A coincide.
# =============================================================================

# Lambda_ij = c X_ij / sum_k X_ik : equal interbank leverage c for every system,
# so that the spectral radius of Lambda equals c in every network.
debtrank_lambda <- function(X, c) {
  rs <- rowSums(X)
  if (any(rs == 0)) stop("every system needs at least one exposure")
  c * X / rs
}

# Fixed point h = min(1, h0 + Lambda h) for the failure of each system in turn.
# Returns the size-weighted contagion loss sum_j w_j L_j and the vector L.
debtrank_loss <- function(X, w, c, tol = 1e-12, max_iter = 10000) {
  Lam <- debtrank_lambda(X, c)
  N <- nrow(X); e <- w                      # equity proportional to size
  H0 <- diag(N); H <- H0
  for (it in seq_len(max_iter)) {
    Hn <- pmin(H0 + Lam %*% H, 1)
    if (max(abs(Hn - H)) < tol) { H <- Hn; break }
    H <- Hn
  }
  L <- vapply(seq_len(N), function(j) sum(e[-j] * H[-j, j]) / sum(e[-j]),
              numeric(1))
  list(total = sum(w * L), L = L, iterations = it)
}

# ---- MCMC on directed networks with fixed undirected support ----------------
# State: each undirected edge {i, j} is reciprocal, i -> j or j -> i.
# The number of reciprocal pairs is held fixed, so density, baseline index and
# the Garlaschelli-Loffredo reciprocity are constant along the chain.
# Moves (symmetric proposals):
#   (a) swap the status of a reciprocal and a one-way pair (new one-way
#       direction drawn at random);  (b) reverse a one-way link.
# Target: pi(X) proportional to exp(beta * Gamma_F(X)), Gamma_F the alignment
# of reciprocal exposures with w.  Constraint: every row of X non-empty.
reciprocity_mcmc <- function(X0, w, beta, n_draws = 150, burn = 2000,
                             thin = 200, p_swap = 0.5, seed = NULL,
                             require_rows = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  X <- (X0 > 0) * 1; diag(X) <- 0
  W <- tcrossprod(w); diag(W) <- 0
  gammaF <- function(X) { F <- X * t(X); if (sum(F) == 0) return(NA)
    (sum(F * W) / sum(F)) / (sum(W) / (nrow(X) * (nrow(X) - 1))) }
  pairs <- which(upper.tri(X) & (X + t(X)) > 0, arr.ind = TRUE)
  status <- function(X) ifelse(X[cbind(pairs[, 1], pairs[, 2])] == 1 &
                                 X[cbind(pairs[, 2], pairs[, 1])] == 1, "R", "O")
  g <- gammaF(X)
  out <- vector("list", n_draws); k <- 0; acc <- 0
  total <- burn + n_draws * thin
  for (s in seq_len(total)) {
    st <- status(X)
    Y <- X
    rec <- which(st == "R"); one <- which(st == "O")
    if (runif(1) < p_swap && length(rec) > 0 && length(one) > 0) {
      r <- rec[sample.int(length(rec), 1)]; o <- one[sample.int(length(one), 1)]
      i <- pairs[r, 1]; j <- pairs[r, 2]
      if (runif(1) < 0.5) Y[i, j] <- 0 else Y[j, i] <- 0
      Y[pairs[o, 1], pairs[o, 2]] <- 1; Y[pairs[o, 2], pairs[o, 1]] <- 1
    } else if (length(one) > 0) {
      o <- one[sample.int(length(one), 1)]
      i <- pairs[o, 1]; j <- pairs[o, 2]
      Y[i, j] <- 1 - X[i, j]; Y[j, i] <- 1 - X[j, i]
    }
    if (!require_rows || all(rowSums(Y) > 0)) {
      gy <- gammaF(Y)
      if (log(runif(1)) < beta * (gy - g)) { X <- Y; g <- gy; acc <- acc + 1 }
    }
    if (s > burn && (s - burn) %% thin == 0) { k <- k + 1; out[[k]] <- X }
  }
  attr(out, "acceptance") <- acc / total
  out
}
