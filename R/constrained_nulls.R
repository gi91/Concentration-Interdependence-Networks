# =============================================================================
# constrained_nulls.R -- degree-conditional nulls with fixed reciprocity
# (Supplementary Definition S1(c) and (d)); exact rewiring (S1(a)) as a chain.
#
# The sampler is written in C++ (src/constrained_mcmc.cpp, compiled once per
# session with Rcpp; on Windows this needs Rtools).  Each chain starts from the
# observed network, uses 500|E| proposals of burn-in and 100|E| proposals between
# draws (|E| = number of links), and draws its random numbers from R's
# generator, so set.seed() makes the results reproducible, also in parallel.
# =============================================================================
.cn_env <- new.env()
cn_compile <- function() {
  if (!isTRUE(.cn_env$ok)) {
    if (!requireNamespace("Rcpp", quietly = TRUE))
      stop("package Rcpp is required (install.packages('Rcpp'); on Windows also Rtools)")
    Rcpp::sourceCpp(file.path("src", "constrained_mcmc.cpp"), env = .cn_env)
    .cn_env$ok <- TRUE
  }
  invisible(TRUE)
}
# TRUE if the C++ sampler can be compiled (Rcpp and a compiler are available)
cn_available <- function() isTRUE(tryCatch(cn_compile(), error = function(e) {
  message("constrained nulls unavailable: ", conditionMessage(e)); FALSE }))
NULL_MODE <- c(rewire = 0L, pairs = 1L, partners = 2L)

# one chain; A binary (any direction convention), w weights summing to one
constrained_chain <- function(A, w, null = c("pairs", "partners", "rewire"),
                              n_draws = 1000, burn_mult = 500, thin_mult = 100,
                              keep_mean = FALSE, keep_codes = FALSE) {
  cn_compile(); null <- match.arg(null)
  A <- (A > 0) * 1L; diag(A) <- 0L; storage.mode(A) <- "integer"
  .cn_env$constrained_chain_cpp(A, as.numeric(w), NULL_MODE[[null]], as.integer(n_draws),
                                burn_mult, thin_mult, keep_mean, keep_codes)
}

# test of psi_F^(1) = w' (A * t(A)) w / (1 - HHI) against n_chains independent
# chains; p = (1 + #{T_b >= T_0}) / (B + 1) over all draws, and per chain
constrained_test <- function(A, w, null = c("pairs", "partners", "rewire"),
                             n_chains = 4, n_draws = 1000, seed = 1,
                             burn_mult = 500, thin_mult = 100, keep_mean = FALSE,
                             cores = as.integer(Sys.getenv("NCI_CORES", 1))) {
  null <- match.arg(null); cn_compile()
  A <- (A > 0) * 1; diag(A) <- 0
  T0 <- psi_index(w, A * t(A))
  run1 <- function(k) { set.seed(seed * 1000 + k)
    constrained_chain(A, w, null, n_draws, burn_mult, thin_mult, keep_mean) }
  ch <- if (cores > 1 && .Platform$OS.type == "unix")
    parallel::mclapply(seq_len(n_chains), run1, mc.cores = min(cores, n_chains))
  else lapply(seq_len(n_chains), run1)
  tol <- 1e-12 * max(1, abs(T0))
  p_chain <- vapply(ch, function(x) (1 + sum(x$psi >= T0 - tol)) / (length(x$psi) + 1), 0)
  allp <- unlist(lapply(ch, `[[`, "psi"))
  # checks: degrees and the constrained quantity are preserved
  kin <- colSums(A); kout <- rowSums(A)
  stopifnot(all(vapply(ch, function(x) all(x$k_out == kout) && all(x$k_in == kin), TRUE)))
  if (null == "pairs") stopifnot(all(unlist(lapply(ch, `[[`, "n_pairs")) == sum(A * t(A)) / 2))
  list(T0 = T0, p = (1 + sum(allp >= T0 - tol)) / (length(allp) + 1),
       p_chain = p_chain, p_range = range(p_chain), null_mean = mean(allp),
       z = (T0 - mean(allp)) / sd(allp), B = length(allp),
       accept = vapply(ch, `[[`, 0, "accept"),
       mean_rec = if (keep_mean) Reduce(`+`, lapply(ch, `[[`, "mean_rec")) / n_chains)
}

# node-level reciprocal embeddedness minus its mean under the constrained null
excess_embeddedness_constrained <- function(w, A, null = c("pairs", "partners"),
                                            n_chains = 4, n_draws = 250, seed = 7, ...) {
  null <- match.arg(null)
  A <- (A > 0) * 1; diag(A) <- 0
  ct <- constrained_test(A, w, null, n_chains = n_chains, n_draws = n_draws,
                         seed = seed, keep_mean = TRUE, ...)
  obs <- embeddedness(w, A * t(A)); m <- embeddedness(w, ct$mean_rec)
  data.frame(e_obs = obs, e_null_mean = m, e_excess = obs - m)
}
