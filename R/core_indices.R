# =============================================================================
# core_indices.R
# Core library for "Beyond Direct Connectivity: Higher-Order and Dynamic
# Concentration in Weighted Networks" (NCI part II).
#
# CONVENTION (as in Sections 3-4 of the paper): A[i, j] > 0 denotes a
# relationship FROM node j TO node i.  Column j = what originates from j.
# All interaction matrices M have zero diagonal and entries in [0, 1]
# (unless stated otherwise).  Weights w are non-negative and sum to one.
# =============================================================================

# ---- helpers ----------------------------------------------------------------

.check_w <- function(w) {
  stopifnot(is.numeric(w), all(w >= 0))
  if (abs(sum(w) - 1) > 1e-10) stop("weights must sum to one")
  if (sum(w > 0) < 2) stop("at least two components of w must be positive")
  invisible(TRUE)
}

.zero_diag <- function(M) {
  if (any(diag(M) != 0)) {
    warning("non-zero diagonal set to zero")
    diag(M) <- 0
  }
  M
}

sym_part <- function(M) (M + t(M)) / 2

# undirected binary support of a (possibly directed, weighted) matrix
support_undirected <- function(A) {
  S <- ((A + t(A)) > 0) * 1
  diag(S) <- 0
  S
}

hhi <- function(w) sum(w^2)

# ---- Section 1: Psi and the availability-alignment factorization -----------

# Psi(w; M, B_C) = w'Mw / (1 - sum w^2)  (eq. 1 with complete benchmark)
psi_index <- function(w, M) {
  .check_w(w)
  M <- .zero_diag(M)
  as.numeric(crossprod(w, M %*% w)) / (1 - sum(w^2))
}

# Proposition 1: Psi = mbar_M * Gamma_M(w)
factorize <- function(w, M) {
  .check_w(w)
  M <- .zero_diag(M)
  N <- length(w)
  K <- N * (N - 1)
  sM <- sum(M)
  W <- tcrossprod(w); diag(W) <- 0
  xbar <- sum(W) / K                          # = (1 - sum w^2) / K
  mbar <- sM / K
  Gamma <- if (sM > 0) (sum(M * W) / sM) / xbar else NA_real_
  psi <- sum(M * W) / sum(W)
  list(psi = psi, mbar = mbar, Gamma = Gamma)
}

alignment <- function(w, M) factorize(w, M)$Gamma

# ---- Node-level contributions (Section 7.3) --------------------------------
# e_i(M) = sum_j M_ij w_j / (1 - w_i);  Psi = sum_i w_i (1 - w_i) e_i / (1 - HHI)
embeddedness <- function(w, M) {
  M <- .zero_diag(M)
  as.numeric(M %*% w) / (1 - w)
}

# ---- Section 2: higher-order topological concentration ---------------------

# BFS from every node on the undirected support of A.  Returns the geodesic
# distance matrix D, the number of shortest paths S and the strength-adjusted
# intensity Omega (average product of A along the shortest paths).
geodesic_structure <- function(A) {
  Au <- pmax(A, t(A))             # symmetric intensities on the support
  diag(Au) <- 0
  N <- nrow(Au)
  nb <- lapply(seq_len(N), function(i) which(Au[i, ] > 0))
  D <- matrix(Inf, N, N); S <- matrix(0, N, N); Th <- matrix(0, N, N)
  for (s in seq_len(N)) {
    dist <- rep(Inf, N); sig <- rep(0, N); th <- rep(0, N)
    dist[s] <- 0; sig[s] <- 1; th[s] <- 1
    frontier <- s; d <- 0
    while (length(frontier) > 0) {
      nxt <- integer(0)
      for (u in frontier) for (v in nb[[u]]) {
        if (dist[v] == Inf) { dist[v] <- d + 1; nxt <- c(nxt, v) }
        if (dist[v] == d + 1) {
          sig[v] <- sig[v] + sig[u]
          th[v]  <- th[v] + th[u] * Au[u, v]
        }
      }
      frontier <- unique(nxt); d <- d + 1
    }
    D[s, ] <- dist; S[s, ] <- sig; Th[s, ] <- th
  }
  Omega <- ifelse(S > 0, Th / pmax(S, 1), 0)
  diag(D) <- 0; diag(Omega) <- 0
  list(D = D, S = S, Omega = Omega)
}

# Order-h matrices M^(h); spec = "binary" or "strength"
order_matrices <- function(A, spec = c("binary", "strength"), gs = NULL) {
  spec <- match.arg(spec)
  if (is.null(gs)) gs <- geodesic_structure(A)
  D <- gs$D
  hmax <- max(D[is.finite(D)])
  Ms <- lapply(seq_len(hmax), function(h) {
    M <- (D == h) * if (spec == "binary") 1 else gs$Omega
    diag(M) <- 0; M
  })
  names(Ms) <- paste0("h", seq_len(hmax))
  Minf <- (D == Inf) * 1; diag(Minf) <- 0
  list(M = Ms, Minf = Minf, D = D, hmax = hmax)
}

# Topological concentration profile, Gamma_h, economic mean distance
topological_profile <- function(w, A, spec = "binary") {
  om <- order_matrices(A, spec)
  tab <- do.call(rbind, lapply(seq_along(om$M), function(h) {
    f <- factorize(w, om$M[[h]])
    data.frame(h = h, psi = f$psi, mbar = f$mbar, Gamma = f$Gamma)
  }))
  connected <- all(om$Minf == 0)
  fin <- factorize(w, om$Minf)
  # with disconnected networks, means are computed over connected pairs
  dw   <- sum(tab$h * tab$psi)  / sum(tab$psi)
  dbar <- sum(tab$h * tab$mbar) / sum(tab$mbar)
  list(profile = tab, dw = dw, dbar = dbar, connected = connected,
       psi_inf = if (connected) 0 else fin$psi,
       mbar_inf = if (connected) 0 else fin$mbar)
}

aggregate_topological <- function(w, A, lambda, spec = "binary") {
  om <- order_matrices(A, spec)
  H <- length(lambda); stopifnot(abs(sum(lambda) - 1) < 1e-12)
  M <- Reduce(`+`, Map(`*`, om$M[seq_len(min(H, om$hmax))],
                       lambda[seq_len(min(H, om$hmax))]))
  psi_index(w, M)
}

# ---- Section 3: propagation concentration ----------------------------------

propagation_matrix <- function(A) {
  s_out <- colSums(A)
  if (any(s_out <= 0)) stop("every node needs positive outgoing strength")
  sweep(A, 2, s_out, "/")
}

mat_pow <- function(P, k) {
  R <- diag(nrow(P)); for (i in seq_len(k)) R <- R %*% P; R
}

propagation_measures <- function(w, A, k) {
  P  <- propagation_matrix(A)
  Pk <- mat_pow(P, k)
  Y  <- colSums(Pk^2)
  q  <- colSums(A) / sum(A)
  CP <- sum(q * Y)
  N  <- length(w)
  O  <- crossprod(Pk); diag(O) <- 0
  fO <- factorize(w, O)
  wk <- as.numeric(Pk %*% w)
  own  <- sum(w^2 * Y)
  conv <- (1 - sum(w^2)) * fO$psi
  list(k = k, Y = Y, CP = CP, CP_tilde = (N * CP - 1) / (N - 1),
       psiO = fO$psi, mbarO = fO$mbar, GammaO = fO$Gamma,
       wk = wk, hhi_wk = sum(wk^2), own = own, conv = conv,
       share_conv = conv / sum(wk^2),
       identity_error = sum(wk^2) - own - conv)    # Prop. 3(i): must be 0
}

stationary_dist <- function(P) {
  e <- eigen(P)
  v <- Re(e$vectors[, which.max(Re(e$values))])
  pi <- v / sum(v)
  slem <- sort(Mod(e$values), decreasing = TRUE)[2]
  list(pi = pi, slem = slem, hhi_pi = sum(pi^2))
}

# ---- Section 4: structural concentration -----------------------------------

kappa_Nh <- function(N, h) exp(lfactorial(N - 2) - lfactorial(N - h - 1))

# Total order-h simple-path interaction Tbar^(h) by depth-first enumeration.
# Path j = v0 -> v1 -> ... -> vh = i with intensity prod A[v_{l+1}, v_l].
simple_path_matrix <- function(A, h) {
  N <- nrow(A)
  out_nb <- lapply(seq_len(N), function(j) which(A[, j] > 0))
  Tb <- matrix(0, N, N)
  rec <- function(path, val, depth, origin) {
    u <- path[length(path)]
    for (v in out_nb[[u]]) {
      if (v %in% path) next
      nv <- val * A[v, u]
      if (depth + 1 == h) Tb[v, origin] <<- Tb[v, origin] + nv
      else rec(c(path, v), nv, depth + 1, origin)
    }
  }
  for (j in seq_len(N)) rec(j, 1, 0, j)
  diag(Tb) <- 0
  Tb
}

# Closed forms of Section 4 (h = 2, 3); fall back to DFS for h >= 4
path_matrix <- function(A, h, method = c("closed", "dfs")) {
  method <- match.arg(method)
  diag(A) <- 0
  if (h == 1) return(A)
  if (method == "dfs" || h >= 4) return(simple_path_matrix(A, h))
  A2 <- A %*% A
  if (h == 2) { T2 <- A2; diag(T2) <- 0; return(T2) }
  A3 <- A2 %*% A
  d2 <- diag(A2)
  T3 <- A3 - A * outer(d2, d2, "+") + A^2 * t(A)
  diag(T3) <- 0
  T3
}

# Structural components for orders 1..H.  Returns per-order indices,
# alignment factors, availabilities and the aggregate decomposition.
structural_decomposition <- function(w, A, H = 3, alpha = rep(1 / H, H),
                                     method = "closed") {
  stopifnot(length(alpha) == H, abs(sum(alpha) - 1) < 1e-12,
            all(A >= 0), all(A <= 1))
  N <- nrow(A); diag(A) <- 0
  At <- t(A)
  rows <- list(); mats <- list()
  for (h in seq_len(H)) {
    Th <- path_matrix(A, h, method) / kappa_Nh(N, h)
    Rh <- (1 - At) * Th
    Fh <- At * Th
    for (X in c("T", "R", "F")) {
      M <- sym_part(switch(X, T = Th, R = Rh, F = Fh)); diag(M) <- 0
      f <- factorize(w, M)
      rows[[length(rows) + 1]] <- data.frame(h = h, comp = X, psi = f$psi,
                                             mbar = f$mbar, Gamma = f$Gamma)
      mats[[paste0(X, h)]] <- M
    }
  }
  tab <- do.call(rbind, rows)
  g <- function(X, hs) sum(alpha[hs] * tab$psi[tab$comp == X & tab$h %in% hs])
  psiTOT <- g("T", seq_len(H))
  parts <- c(one_way    = g("R", 1),
             reciprocal = g("F", 1),
             open_chains = if (H >= 2) g("R", 2:H) else 0,
             circuits    = if (H >= 2) g("F", 2:H) else 0)
  list(table = tab, psiTOT = psiTOT, parts = parts,
       shares = parts / psiTOT,
       residual = psiTOT - sum(parts),     # Prop. 4(ii): must be 0
       mats = mats)
}

# ---- Section 5: dynamic concentration --------------------------------------

# Proposition 5 (nested decomposition, as in the paper) and, for comparison,
# the two-player Shapley decomposition of Psi over the primitive changes
# (w, M), where the network effect includes availability.
dynamic_decomposition <- function(w0, M0, w1, M1) {
  G <- function(w, M) factorize(w, M)$Gamma
  m0 <- factorize(w0, M0)$mbar; m1 <- factorize(w1, M1)$mbar
  G00 <- G(w0, M0); G11 <- G(w1, M1); G10 <- G(w1, M0); G01 <- G(w0, M1)
  ms <- (m0 + m1) / 2
  d_a <- 0.5 * (m1 - m0) * (G11 + G00)
  d_w <- ms / 2 * (G10 - G00 + G11 - G01)
  d_M <- ms / 2 * (G01 - G00 + G11 - G10)
  psi0 <- m0 * G00; psi1 <- m1 * G11
  # full two-player Shapley on Psi(w, M)
  P <- function(w, M) psi_index(w, M)
  sh_w <- 0.5 * ((P(w1, M0) - P(w0, M0)) + (P(w1, M1) - P(w0, M1)))
  sh_M <- 0.5 * ((P(w0, M1) - P(w0, M0)) + (P(w1, M1) - P(w1, M0)))
  data.frame(psi0 = psi0, psi1 = psi1, delta = psi1 - psi0,
             availability = d_a, weights = d_w, network = d_M,
             residual = psi1 - psi0 - d_a - d_w - d_M,
             shapley_weights = sh_w, shapley_network = sh_M)
}

# Lag-h intertemporal transmission concentration (bilinear form)
psi_intertemporal <- function(w_t, w_lag, G) {
  diag(G) <- 0
  num <- as.numeric(crossprod(w_t, G %*% w_lag))
  den <- 1 - sum(w_t * w_lag)
  W <- outer(w_t, w_lag); diag(W) <- 0
  K <- length(w_t) * (length(w_t) - 1)
  sG <- sum(G)
  list(psi = num / den, mbar = sG / K,
       Gamma = if (sG > 0) (sum(G * W) / sG) / (sum(W) / K) else NA_real_)
}

# ---- convenience: baseline first-order quantities --------------------------
baseline_index <- function(w, A) {
  S <- support_undirected(A)
  f <- factorize(w, S)
  list(psiD = f$psi, density = f$mbar, Gamma1 = f$Gamma)
}
