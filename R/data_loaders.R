# =============================================================================
# data_loaders.R -- standardized inputs used by every analysis script
#
# The analysis scripts read ONLY the processed files below.  Raw-data readers
# are thin adapters: if the layout of your raw files differs, change the
# adapter, not the analysis code.
#
#  data/processed/wiod_<year>.rds  list(Z, go, va, country, industry)
#     Z[a, b]  intermediate inputs supplied by country-industry a to b
#              (WIOD convention: rows = supplier, columns = user), ROW excluded
#     go, va   gross output and value added of each country-industry
#  data/processed/bis_claims.rds    data.table(quarter, creditor, counterparty,
#                                              claims)   -- reporters x reporters
#  data/processed/bis_totals.rds    data.table(quarter, creditor, total_claims)
#                                   -- total foreign claims on ALL counterparties
#  data/processed/srisk.rds         data.table(date, country, srisk)  US$ bn
# =============================================================================

# ---- WIOD 2016 release (R format: WIOT<year>_October16_ROW.RData) -----------
# ASSUMED layout (verify on your files): object `wiot` with identifier columns
# IndustryCode, IndustryDescription, Country, RNr, Year; then use columns named
# <ISO3><1..56> (intermediate) and <ISO3><57..61> (final), plus TOT; rows with
# RNr 1..2464 are country-industries, and rows whose IndustryCode is "VA" and
# "GO" contain value added and gross output.
read_wiod_rdata <- function(file, drop = "ROW", va_code = "VA",
                            go_code = "GO") {
  env <- new.env(); load(file, envir = env)
  wiot <- get(ls(env)[1], envir = env)
  wiot <- as.data.frame(wiot)
  ci <- wiot[wiot$RNr <= 2464, ]
  keep <- ci$Country != drop
  cols <- paste0(ci$Country, ((ci$RNr - 1) %% 56) + 1)
  stopifnot(all(cols %in% names(wiot)))
  Z <- as.matrix(ci[keep, cols[keep]])
  va <- as.numeric(wiot[wiot$IndustryCode == va_code, cols[keep]])
  go <- as.numeric(wiot[wiot$IndustryCode == go_code, cols[keep]])
  stopifnot(length(va) == nrow(Z), length(go) == nrow(Z))
  list(Z = unname(Z), go = go, va = va, country = ci$Country[keep],
       industry = ci$IndustryCode[keep])
}

# Production-network source: NCI_PROD = "wiod" (2016 release, 2000-2014) or
# "icio" (OECD ICIO, 1995-2022).  Both are stored in the same processed format.
# default: OECD ICIO when its processed files exist, WIOD otherwise
PROD_SOURCE <- Sys.getenv("NCI_PROD", if (length(list.files(file.path("data", "processed"),
  "^icio_[0-9]{4}\\.rds$"))) "icio" else "wiod")
load_wiod_year <- function(year, source = PROD_SOURCE)
  readRDS(file.path(DIR_PROC, sprintf("%s_%d.rds", source, year)))

# ---- OECD ICIO (one CSV per year, e.g. inside the per-period zip files) -------
# ASSUMED layout (verify on the first file): first column = row labels; rows
# and columns of the intermediate block labelled <ISO3>_<industry code>
# (e.g. "AUS_01T02"); final-demand columns <ISO3>_<HFCE|NPISH|GGFC|GFCF|INVNT|
# DPABR>; rows "VA" (value added) and "OUT" (output).  China and Mexico are
# split into CN1/CN2 and MX1/MX2 (processing trade): they are summed back
# into CHN and MEX.  The rest of the world (ROW) is dropped.
read_icio_csv <- function(file, drop = "ROW", va_row = "VA", out_row = "OUT") {
  d <- fread(file, header = TRUE)
  lab <- as.character(d[[1]]); d <- as.matrix(d[, -1, with = FALSE])
  rownames(d) <- lab
  is_ci <- function(x) grepl("^[A-Z0-9]{3}_", x) &
    !grepl("_(HFCE|NPISH|GGFC|GFCF|INVNT|DPABR|TAXSUB|TLS|VA|OUT)$", x)
  rci <- lab[is_ci(lab)]; cci <- colnames(d)[is_ci(colnames(d))]
  stopifnot(length(rci) > 0, setequal(rci, cci))
  if (!all(c(va_row, out_row) %in% lab))
    stop("rows '", va_row, "' / '", out_row, "' not found: check the labels")
  # map split economies back to ISO3
  iso <- function(x) { cty <- sub("_.*$", "", x)
    cty[cty %in% c("CN1", "CN2")] <- "CHN"; cty[cty %in% c("MX1", "MX2")] <- "MEX"; cty }
  ind <- sub("^[^_]*_", "", rci)
  key <- paste(iso(rci), ind, sep = "_")
  keep <- iso(rci) != drop
  rk <- rci[keep]; kk <- key[keep]
  ukey <- unique(kk)
  Zr <- d[rk, rk]
  if (anyDuplicated(kk)) {               # merge CN1/CN2, MX1/MX2 (O(N^2) sums)
    Z <- t(rowsum(t(rowsum(Zr, kk, reorder = FALSE)), kk, reorder = FALSE))
    va <- as.numeric(rowsum(d[va_row, rk], kk, reorder = FALSE))
    go <- as.numeric(rowsum(d[out_row, rk], kk, reorder = FALSE))
  } else {                               # "SML" files: already merged
    Z <- Zr; va <- as.numeric(d[va_row, rk]); go <- as.numeric(d[out_row, rk])
  }
  if (any(va > go + 1e-6)) warning("value added exceeds output for some nodes")
  list(Z = unname(Z), go = go, va = va, country = sub("_.*$", "", ukey),
       industry = sub("^[^_]*_", "", ukey))
}

# Aggregation to the country level: Zc[j, i] supplied by country j to i,
# X = gross output by country, VA by country.
wiod_country <- function(d) {
  cc <- sort(unique(d$country))
  G <- sapply(cc, function(x) as.numeric(d$country == x))   # CI x C
  Zc <- t(G) %*% d$Z %*% G
  dimnames(Zc) <- list(cc, cc)
  list(Z = Zc, X = as.numeric(t(G) %*% d$go), VA = as.numeric(t(G) %*% d$va),
       country = cc)
}

# Structural adjacency of Section 7.1 with convention A[i, j]: link j -> i.
# intensity: "share" = Z_ji / X_i ; "neutral" = Z_ji / sqrt(X_i X_j)
wiod_intensity <- function(dc, type = c("share", "neutral")) {
  type <- match.arg(type)
  Zji <- t(dc$Z)                                 # [i, j] = Z_ji
  M <- if (type == "share") Zji / dc$X else Zji / sqrt(outer(dc$X, dc$X))
  diag(M) <- 0
  M
}

# ---- BIS consolidated banking statistics ------------------------------------
# Template for the flat CSV downloaded from data.bis.org (CBS, foreign claims,
# ultimate-risk basis, all sectors).  Column names differ across downloads:
# map them in `colmap`.
read_bis_cbs_csv <- function(file, colmap = c(quarter = "TIME_PERIOD",
                              creditor = "L_REP_CTY", counterparty = "L_CP_COUNTRY",
                              claims = "OBS_VALUE"), total_code = "5J") {
  d <- fread(file)
  setnames(d, colmap, names(colmap))
  d <- d[, names(colmap), with = FALSE]
  d[, quarter := as.Date(paste0(substr(quarter, 1, 4), "-",
                         (as.integer(substr(quarter, 7, 7)) - 1) * 3 + 1, "-01"))]
  totals <- d[counterparty == total_code, .(quarter, creditor,
                                             total_claims = claims)]
  list(claims = d[counterparty != total_code], totals = totals)
}

# Claims matrix C[i, j] (claims of system i on country j) and weights at date q
bis_snapshot <- function(claims, totals, q, reporters) {
  cl <- claims[quarter == q & creditor %in% reporters &
                 counterparty %in% reporters]
  C <- matrix(0, length(reporters), length(reporters),
              dimnames = list(reporters, reporters))
  C[cbind(match(cl$creditor, reporters), match(cl$counterparty, reporters))] <-
    cl$claims
  diag(C) <- 0
  tot <- totals[quarter == q][match(reporters, creditor), total_claims]
  list(C = C, total = tot, w = tot / sum(tot), quarter = q)
}

# Exposure intensity: share of foreign claims of r on c.  Structural matrix
# A[r, c] = 1(C[r, c] / total_r >= theta)  (link c -> r, r exposed to c).
bis_adjacency <- function(snap, theta = 0.02, weighted = FALSE) {
  S <- snap$C / snap$total
  if (weighted) S else (S >= theta) * 1
}


# ---- BIS bulk files in "csv_col" (wide) format: WS_CBS_PUB / WS_LBS_D_PUB ----
# One row per series, dimension codes in UPPERCASE columns (each followed by
# its label column), one column per quarter.  `filter` is a named list of
# dimension = code; quarters from `from_year`.  Returns a long data.table.
read_bis_wide <- function(zipfile, filter, from_year = 2007) {
  # portable (Windows too): extract to a temporary folder with R's own unzip
  tmp <- tempfile(); dir.create(tmp)
  csv <- file.path(tmp, utils::unzip(zipfile, exdir = tmp)[1] |> basename())
  on.exit(unlink(tmp, recursive = TRUE))
  hdr <- names(fread(csv, nrows = 0))
  dims <- hdr[grepl("^[A-Z_]+$", hdr) &
              !hdr %in% c("FREQ", "TIME_FORMAT", "COLLECTION", "ORG_VISIBILITY")]
  qs <- hdr[grepl("^[0-9]{4}-Q[1-4]$", hdr) & as.integer(substr(hdr, 1, 4)) >= from_year]
  # na.strings = "": the ISO code of Namibia ("NA") must not be read as missing
  d <- fread(csv, select = c(dims, qs), na.strings = "",
             showProgress = FALSE, colClasses = list(character = dims))
  miss <- setdiff(names(filter), dims)
  if (length(miss)) warning("filter dimensions not in file: ", paste(miss, collapse = ", "))
  for (v in intersect(names(filter), dims)) d <- d[get(v) == filter[[v]]]
  long <- melt(d, id.vars = dims, measure.vars = qs, variable.name = "period",
               value.name = "value", na.rm = TRUE)
  long[, quarter := as.Date(paste0(substr(period, 1, 4), "-",
                      (as.integer(substr(period, 7, 7)) - 1) * 3 + 1, "-01"))]
  long[]
}

# claims/totals in the standard processed format from a long BIS table
bis_standardize <- function(long, reporters, total_code = "5J") {
  cl <- long[L_REP_CTY %in% reporters & L_CP_COUNTRY %in% reporters &
             L_CP_COUNTRY != L_REP_CTY,
             .(quarter, creditor = L_REP_CTY, counterparty = L_CP_COUNTRY, claims = value)]
  tt <- long[L_REP_CTY %in% reporters & L_CP_COUNTRY == total_code,
             .(quarter, creditor = L_REP_CTY, total_claims = value)]
  list(claims = cl, totals = tt)
}
