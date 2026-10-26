# =============================================================================
# 05_srisk_validation.R -- Section 7.3: Tables 5 and 6
# =============================================================================
source("R/utils.R"); source("R/core_indices.R"); source("R/validation.R")
source("R/null_models.R"); source("R/fixes.R"); source("R/data_loaders.R")
source("R/constrained_nulls.R")   # rows (12c)-(12d): Definition S1(c)-(d), needs Rcpp
CN_OK <- cn_available()

HC <- Sys.getenv("NCI_HC", "HC1")         # robust covariance used in the tables
B_FL <- as.integer(Sys.getenv("NCI_B_FL", 5000))
B_EX <- as.integer(Sys.getenv("NCI_B_EX", 1000))   # rewirings for excess embeddedness
bis <- readRDS(file.path(DIR_PROC, "results_bis.rds"))
sr  <- readRDS(file.path(DIR_PROC, "srisk.rds"))
REP <- bis$reporters
# Sample of the regressions.  The paper uses 18 systems (Ireland had no SRISK
# series at the time); V-Lab now publishes Ireland.  Default = sample of the paper;
# set NCI_SRISK_EXCLUDE="" to include Ireland.  Row (14) reports the other choice.
EXCL <- strsplit(Sys.getenv("NCI_SRISK_EXCLUDE", "IE"), ",")[[1]]
ext_ok <- file.exists(file.path(DIR_PROC, "bis_claims_ext.rds"))
if (ext_ok) { cl_ext <- readRDS(file.path(DIR_PROC, "bis_claims_ext.rds"))
              tt_ext <- readRDS(file.path(DIR_PROC, "bis_totals_ext.rds")) }
qdate <- function(q) as.Date(paste0(substr(q, 1, 4), "-",
                                    (as.integer(substr(q, 7, 7)) - 1) * 3 + 1, "-01"))

WINDOWS <- list(
  GFC   = list(q = "2007-Q2", peak = c("2007-07-01", "2009-06-01")),
  EA    = list(q = "2010-Q1", peak = c("2010-04-01", "2012-12-01")),
  COVID = list(q = "2019-Q4", peak = c("2020-01-01", "2020-12-01")),
  Calm  = list(q = "2014-Q2", peak = c("2014-07-01", "2016-06-01")))
WIN_LAB <- c(GFC = "Global financial crisis 2007--2009",
             EA = "Euro-area crisis 2010--2012", COVID = "COVID-19 2020",
             Calm = "Calm period 2014--2016")

# Missing months are treated as zero shortfall, as in the paper.
# NOTE: this conflates "not reported" with "no shortfall"; see the
# robustness block B5 in 08_sim_robustness.R and consider a robustness check
# that drops country-windows with missing peak months.
sr_full <- CJ(date = sort(unique(sr$date)), country = unique(sr$country))
sr_full <- merge(sr_full, sr, all.x = TRUE)[is.na(srisk), srisk := 0]

# circuit embeddedness: sum of the circuit matrices of orders h = 2, 3 (circuits
# of length three and four); verified on the real data (row 5 of the paper).
build_data <- function(win, excl = EXCL, with_india = FALSE) {
  m <- bis$mats[[win$q]]; REPS <- REP
  if (with_india) {                      # 20-node network including India
    REPS <- sort(c(REP, "IN"))
    sn <- bis_snapshot(cl_ext, tt_ext, qdate(win$q), REPS)
    m <- list(w = sn$w, A = bis_adjacency(sn, 0.02), C = sn$C, total = sn$total,
              S = bis_adjacency(sn, weighted = TRUE))
  }
  sdc <- structural_decomposition(m$w, m$A, H = 3)
  M_dir <- support_undirected(m$A)
  M_rec <- sdc$mats$F1
  M_cir <- sdc$mats$F2 + sdc$mats$F3
  p0 <- as.Date(win$peak[1]); p1 <- as.Date(win$peak[2])
  pre <- seq(p0, by = "-1 month", length.out = 13)[-1]
  peak <- sr_full[date >= p0 & date <= p1, .(peak = max(srisk)), by = country]
  prev <- sr_full[date %in% pre, .(pre = mean(srisk)), by = country]
  d <- data.frame(country = REPS, fc = m$total,
                  e_dir = embeddedness(m$w, M_dir),
                  e_rec = embeddedness(m$w, M_rec),
                  e_cir = embeddedness(m$w, M_cir),
                  # Online Appendix OA.3: weighted reciprocal embeddedness, on the
                  # mutual exposure shares S_ij S_ji (no threshold)
                  e_rec_w = embeddedness(m$w, m$S * t(m$S)),
                  # reciprocal embeddedness net of its degree-conditional mean
                  e_rec_ex = excess_embeddedness(m$w, m$A, B = B_EX, seed = 7)$e_excess,
                  centralities(m$A, m$C, total_claims = sum(m$total)))
  d <- merge(d, merge(peak, prev, by = "country"), by = "country")
  d <- d[!d$country %in% excl, ]
  d$y <- log1p(d$peak); d$lfc <- log(d$fc); d$lpre <- log1p(d$pre)
  # row 9: SRISK (US$ bn) relative to foreign claims (BIS: US$ mn -> bn)
  d$y_rel <- log1p(d$peak / (d$fc / 1000)); d$lpre_rel <- log1p(d$pre / (d$fc / 1000))
  d
}

ctrl <- "lfc + lpre"
spec <- function(rhs) as.formula(paste("y ~", ctrl, if (nzchar(rhs)) "+", rhs))
res5 <- list(); res6 <- NULL
for (nm in names(WINDOWS)) {
  d <- build_data(WINDOWS[[nm]])
  r1 <- robust_ols(spec(""), d, HC)
  r2 <- robust_ols(spec("e_dir"), d, HC)
  r3 <- robust_ols(spec("e_rec"), d, HC)
  r4 <- robust_ols(spec("e_dir + e_rec"), d, HC)
  r5 <- robust_ols(spec("e_cir"), d, HC)
  r6 <- robust_ols(spec("e_rec + I(lfc^2)"), d, HC)
  fl <- freedman_lane(spec("e_rec"), "e_rec", d, B = B_FL, type = HC, seed = 1)
  lo <- range(leave_one_out_t(spec("e_rec"), "e_rec", d, HC))
  r9 <- robust_ols(y_rel ~ lfc + lpre_rel + e_rec, d, HC)
  r11 <- robust_ols(y ~ lfc + e_rec, d, HC)
  r3_hc3 <- robust_ols(spec("e_rec"), d, "HC3")
  # rows (12c)-(12d): excess over 1,000 draws (4 chains x 250) with fixed reciprocal
  # pairs or partners (Definition S1(c)-(d))
  # (skipped, with "--" in the table, if Rcpp cannot compile the sampler)
  if (CN_OK) {
    mq <- bis$mats[[WINDOWS[[nm]]$q]]
    for (nl in c("pairs", "partners")) {
      ex <- excess_embeddedness_constrained(mq$w, mq$A, nl, n_chains = 4, n_draws = 250,
                                            seed = if (nl == "pairs") 71 else 72)$e_excess
      d[[paste0("e_rec_ex_", nl)]] <- ex[match(d$country, REP)]
    }
    r3c <- robust_ols(spec("e_rec_ex_pairs"), d, HC)
    flc <- freedman_lane(spec("e_rec_ex_pairs"), "e_rec_ex_pairs", d, B = B_FL, type = HC, seed = 4)
    r3d <- robust_ols(spec("e_rec_ex_partners"), d, HC)
    fld <- freedman_lane(spec("e_rec_ex_partners"), "e_rec_ex_partners", d, B = B_FL, type = HC, seed = 5)
  }
  r3x <- robust_ols(spec("e_rec_ex"), d, HC)
  r3xc <- robust_ols(spec("e_rec_ex + size_cpty"), d, HC)
  flx <- freedman_lane(spec("e_rec_ex"), "e_rec_ex", d, B = B_FL, type = HC, seed = 2)
  r10 <- if (ext_ok) robust_ols(spec("e_rec"), build_data(WINDOWS[[nm]], with_india = TRUE), HC)
  d14 <- build_data(WINDOWS[[nm]], excl = if (length(EXCL) && nzchar(EXCL[1])) character(0) else "IE")
  r14 <- robust_ols(spec("e_rec"), d14, HC)
  ct <- function(r, v) fmt_coef_t(r$coef[v], r$t[v])
  res5[[nm]] <- c(
    r2_1 = fmt_num(r1$r2, 2),
    dir_2 = ct(r2, "e_dir"), r2_2 = fmt_num(r2$r2, 2),
    rec_3 = ct(r3, "e_rec"), r2_3 = fmt_num(r3$r2, 2),
    dir_4 = ct(r4, "e_dir"), rec_4 = ct(r4, "e_rec"), r2_4 = fmt_num(r4$r2, 2),
    cir_5 = ct(r5, "e_cir"), r2_5 = fmt_num(r5$r2, 2),
    sq_6 = ct(r6, "e_rec"), fl_7 = fmt_num(fl$p, 3),
    loo_8 = sprintf("%s to %s", fmt_num(lo[1], 2), fmt_num(lo[2], 2)),
    rel_9 = ct(r9, "e_rec"),
    india_10 = if (is.null(r10)) "--" else ct(r10, "e_rec"),
    nopre_11 = ct(r11, "e_rec"),
    hc3_3 = ct(r3_hc3, "e_rec"),
    ex_3 = ct(r3x, "e_rec_ex"), ex_fl = fmt_num(flx$p, 3),
    exc_pairs = if (CN_OK) ct(r3c, "e_rec_ex_pairs") else "--",
    exc_pairs_fl = if (CN_OK) fmt_num(flc$p, 3) else "--",
    exc_partners = if (CN_OK) ct(r3d, "e_rec_ex_partners") else "--",
    exc_partners_fl = if (CN_OK) fmt_num(fld$p, 3) else "--",
    exc_3 = ct(r3xc, "e_rec_ex"), exc_c = ct(r3xc, "size_cpty"),
    ex_cor = fmt_num(cor(d$e_rec_ex, d$size_cpty), 2),
    ie_14 = sprintf("%s [n = %d]", ct(r14, "e_rec"), nrow(d14)),
    n = as.character(nrow(d)))
  if (nm == "GFC") {                                  # Table 6
    cents <- c(eig_binary = "Eigenvector, binary", out_degree = "Out-degree",
               pagerank_w = "PageRank, weighted", eig_weighted = "Eigenvector, weighted",
               in_degree = "In-degree", katz = "Katz", size_cpty = "Size as a counterparty")
    res6 <- do.call(rbind, lapply(names(cents), function(cn) {
      a <- robust_ols(spec(cn), d, HC)
      j <- robust_ols(spec(paste(cn, "+ e_rec")), d, HC)
      data.frame(Measure = cents[[cn]], corr = fmt_num(cor(d[[cn]], d$e_rec), 2),
                 alone = ct(a, cn), joint_c = ct(j, cn), joint_r = ct(j, "e_rec"))
    }))
    saveRDS(d, file.path(DIR_PROC, "validation_data_gfc.rds"))
  }
}
# ---- Online Appendix OA.3: weighted reciprocal embeddedness ------------------
oa3 <- do.call(rbind, lapply(names(WINDOWS), function(nm) {
  d <- build_data(WINDOWS[[nm]])
  r <- robust_ols(spec("e_rec_w"), d, HC)
  rj <- robust_ols(spec("e_rec_w + e_rec"), d, HC)
  rc <- robust_ols(spec("e_rec_w + size_cpty"), d, HC)
  fl <- freedman_lane(spec("e_rec_w"), "e_rec_w", d, B = B_FL, type = HC, seed = 3)
  ct <- function(r, v) fmt_coef_t(r$coef[v], r$t[v])
  data.frame(window = WIN_LAB[[nm]], alone = ct(r, "e_rec_w"), r2 = fmt_num(r$r2, 2),
             fl = fmt_num(fl$p, 3), joint_w = ct(rj, "e_rec_w"), joint_b = ct(rj, "e_rec"),
             with_cpty = ct(rc, "e_rec_w"), corr_bin = fmt_num(cor(d$e_rec_w, d$e_rec), 2))
}))
write_tex_table(oa3, "tabOA3_weighted_embeddedness.tex",
  caption = "Weighted reciprocal embeddedness and peak systemic risk",
  label = "tab:oa3", align = "lccccccc",
  header = c("Window", "Weighted emb.", "$R^2$", "Perm.\\ $p$", "Joint: weighted",
             "Joint: binary", "With size as cpty", "Corr.\\ with binary"),
  notes = paste0("Specification (3) of Table 5 with reciprocal embeddedness computed on the",
    " matrix of mutual exposure shares $S_{ij}S_{ji}$, without threshold. ", HC,
    " robust $t$-statistics in parentheses; Freedman--Lane permutation $p$-value."))

R5 <- do.call(cbind, res5)
lab5 <- c("(1) Size and pre-crisis SRISK: $R^2$", "(2) Direct embeddedness", "$R^2$",
          "(3) Reciprocal embeddedness", "$R^2$", "(4) Direct embeddedness",
          "\\phantom{(4)} Reciprocal embeddedness", "$R^2$",
          "(5) Circuit embeddedness", "$R^2$", "(6) With squared log size",
          "(7) Permutation $p$-value", "(8) Leave-one-out $t$, range",
          "(9) SRISK relative to foreign claims", "(10) Including India",
          "(11) Without pre-crisis SRISK", "(3) with HC3 standard errors",
          "(12) Excess reciprocal embeddedness", "\\phantom{(12)} Permutation $p$-value",
          "(12c) Excess, fixed reciprocal pairs", "\\phantom{(12c)} Permutation $p$-value",
          "(12d) Excess, fixed reciprocal partners", "\\phantom{(12d)} Permutation $p$-value",
          "(13) Excess reciprocal embeddedness", "\\phantom{(13)} Size as a counterparty",
          "\\phantom{(13)} Corr.\\ excess emb., size as counterparty",
          if (length(EXCL) && nzchar(EXCL[1])) "(14) Including Ireland" else "(14) Excluding Ireland",
          "Number of banking systems")
tab5 <- data.frame(Specification = lab5, R5, check.names = FALSE)
panel5 <- c(rep("Network embeddedness", 10),
            rep("Robustness of reciprocal embeddedness, specification (3)", 7),
            rep("Reciprocal embeddedness net of its degree-conditional expectation", 9),
            rep("Sample", 2))
write_tex_table(tab5, "tab5_srisk.tex", panel = panel5,
  caption = "Network embeddedness before a crisis and peak systemic risk during the crisis",
  label = "tab:srisk", align = "lcccc",
  header = c("", WIN_LAB[names(WINDOWS)]),
  notes = paste0("Cross-sectional OLS of $\\log(1+\\text{peak SRISK})$ on log foreign",
    " claims, $\\log(1+\\text{pre-crisis SRISK})$ and embeddedness in the quarter",
    " preceding each window. Coefficients with ", HC, " robust $t$-statistics in",
    " parentheses (last row: HC3). Row (7): two-sided Freedman--Lane permutation",
    " $p$-value, ", B_FL, " permutations. Row (10): 20-node network including India.",
    " SRISK blanks in V-Lab (no shortfall) are set to zero."))
write_tex_table(res6, "tab6_centrality.tex",
  caption = "Reciprocal embeddedness and standard centrality measures, global financial crisis",
  label = "tab:centrality", align = "lcccc",
  header = c("Centrality measure", "Corr.\\ with rec.\\ emb.", "Centrality alone",
             "Joint: centrality", "Joint: rec.\\ emb."),
  notes = "Specification (3) of Table 5, GFC window. In-degree and Katz orientation: see R/validation.R.")
log_msg("Tables 5 and 6 written")
