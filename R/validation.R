# =============================================================================
# validation.R -- Section 7.3: embeddedness, SRISK regressions, centralities
# =============================================================================
suppressPackageStartupMessages({ library(sandwich); library(igraph) })

# OLS with heteroskedasticity-robust t-statistics.
# NB: the paper does not state which HC estimator is used.  With n = 18 and
# 3-5 regressors the size of HC1 depends on the leverage of the design: check
# it on the real design with scripts/08_sim_robustness.R, block B5.  HC3 is
# reported alongside as a check.
robust_ols <- function(formula, data, type = "HC1") {
  fit <- lm(formula, data = data)
  V <- vcovHC(fit, type = type)
  b <- coef(fit); se <- sqrt(diag(V))
  list(fit = fit, coef = b, t = b / se, r2 = summary(fit)$r.squared)
}

# Freedman and Lane (1983) permutation p-value for the t-statistic of `var`
freedman_lane <- function(formula_full, var, data, B = 5000, type = "HC1",
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  full <- robust_ols(formula_full, data, type)
  t_obs <- full$t[var]
  f_red <- update(formula_full, paste(". ~ . -", var))
  red <- lm(f_red, data = data)
  yhat <- fitted(red); res <- resid(red)
  yname <- all.vars(formula_full)[1]
  t_null <- vapply(seq_len(B), function(b) {
    d <- data; d[[yname]] <- yhat + res[sample.int(length(res))]
    robust_ols(formula_full, d, type)$t[var]
  }, numeric(1))
  list(t = t_obs, p = mean(abs(t_null) >= abs(t_obs)))
}

leave_one_out_t <- function(formula, var, data, type = "HC1") {
  vapply(seq_len(nrow(data)), function(i)
    robust_ols(formula, data[-i, ], type)$t[var], numeric(1))
}

# ---- centrality measures of Table 6 -----------------------------------------
# X: binary exposure matrix, X[i, j] = 1 if i is exposed to j (= structural A)
# C: claims matrix, C[i, j] = claims of system i on country j
# Verified on the real data (Table 6 of the paper reproduced, correlations and
# t-statistics identical): in-degree of j = number of reporting systems exposed
# to j (column sums of X); degrees are divided by N - 1, PageRank and Katz by
# their maximum (eigenvector centralities are max-normalised by igraph).
centralities <- function(X, C, total_claims = sum(C)) {
  N <- nrow(X)
  Su <- ((X + t(X)) > 0) * 1; diag(Su) <- 0
  g_u <- graph_from_adjacency_matrix(Su, mode = "undirected")
  # creditor -> debtor orientation for PageRank: edge i -> j weight C[i, j]
  g_c <- graph_from_adjacency_matrix(C * (1 - diag(N)), mode = "directed",
                                     weighted = TRUE)
  Cs <- (C + t(C)) / 2; diag(Cs) <- 0
  g_s <- graph_from_adjacency_matrix(Cs, mode = "undirected", weighted = TRUE)
  indeg <- colSums(X)                      # systems exposed to j
  # Katz on the network of incoming exposures: x = (I - a X')^{-1} 1
  Xin <- t(X); lam <- max(Mod(eigen(Xin, only.values = TRUE)$values))
  katz <- as.numeric(solve(diag(N) - (0.5 / lam) * Xin, rep(1, N)))
  data.frame(
    eig_binary   = eigen_centrality(g_u)$vector,
    out_degree   = rowSums(X) / (N - 1),
    pagerank_w   = { pr <- page_rank(g_c, damping = 0.85)$vector; pr / max(pr) },
    eig_weighted = eigen_centrality(g_s)$vector,
    in_degree    = indeg / (N - 1),
    katz         = katz / max(katz),
    size_cpty    = log1p(1000 * colSums(C * (1 - diag(N))) / total_claims))
}
