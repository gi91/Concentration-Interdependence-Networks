# =============================================================================
# 13_gdp_comovement.R -- Online Appendix OA.4: do circuits carry information on
# the co-movement of real GDP growth?  Pair-level regression (Frankel and Rose,
# 1998) of the correlation of real GDP growth of countries i and j over a window
# on the intensity of the circuits linking them at the start of the window,
# controlling for bilateral trade intensity and country fixed effects.
#
# Input: data/raw/gdp/  World Bank WDI "GDP growth (annual %)", indicator
#        NY.GDP.MKTP.KD.ZG, CSV as downloaded (API_NY.GDP.MKTP.KD.ZG_*.csv:
#        4 header lines, then Country Name, Country Code, ..., 1960, 1961, ...)
#        or a long CSV with columns iso3, year, growth.
# Output: tabOA4_gdp_comovement.tex/csv.  Skipped (with a message) if no GDP file.
# =============================================================================
source("R/utils.R"); source("R/core_indices.R"); source("R/data_loaders.R")
suppressPackageStartupMessages(library(sandwich))

gdir <- file.path(DIR_RAW, "gdp")
files <- list.files(gdir, "\\.csv$", full.names = TRUE)
# the World Bank download is a zip (API_NY.GDP.MKTP.KD.ZG_DS2_..._csv_v2_....zip):
# read the data csv inside it without extracting by hand
for (z in list.files(gdir, "\\.zip$", full.names = TRUE)) {
  inside <- utils::unzip(z, list = TRUE)$Name
  inside <- inside[grepl("^API_.*\\.csv$", basename(inside))]
  if (length(inside)) files <- c(files, utils::unzip(z, files = inside, exdir = tempdir()))
}
files <- files[!grepl("^Metadata", basename(files))]
# NB: no quit() here -- run() sources the script inside the R session, and
# quit() would close the whole session (RStudio included)
if (!length(files)) {
  log_msg("no GDP file in data/raw/gdp: robustness section on GDP skipped")
} else {

read_gdp <- function(f) {
  first <- readLines(f, n = 5, warn = FALSE)
  if (any(grepl("Country Code", first))) {                # WDI wide format
    x <- fread(f, skip = grep("Country Code", first)[1] - 1, header = TRUE)
    yrs <- grep("^(19|20)[0-9]{2}$", names(x), value = TRUE)
    x[, (yrs) := lapply(.SD, as.numeric), .SDcols = yrs]    # empty years are logical
    melt(x[, c("Country Code", yrs), with = FALSE], id.vars = "Country Code",
         variable.name = "year", value.name = "growth")[,
      .(iso3 = `Country Code`, year = as.integer(as.character(year)), growth = as.numeric(growth))]
  } else fread(f)[, .(iso3, year = as.integer(year), growth = as.numeric(growth))]
}
gdp <- rbindlist(lapply(files, read_gdp))[!is.na(growth)]

# ICIO country codes are ISO3; Taiwan is not in WDI (no World Bank member)
WINDOWS <- list(
  "Full 2000--2014"   = list(base = 2000, years = 2000:2014),
  "Before 2008"       = list(base = 2000, years = 2000:2007),
  "2008--2014"        = list(base = 2007, years = 2008:2014),
  "Extended 1995--2022" = list(base = 1995, years = 1995:2022))

pair_data <- function(win) {
  d  <- load_wiod_year(win$base, "icio"); dc <- wiod_country(d)
  w  <- dc$VA / sum(dc$VA); cc <- dc$country
  S  <- wiod_intensity(dc, "share")                     # weighted: Z_ji / X_i
  sw <- structural_decomposition(w, pmin(S / max(S), 1), H = 3)
  circ <- sw$mats$F2 + sw$mats$F3                       # circuits of length 3 and 4
  A  <- (S >= 0.005) * 1; rec <- A * t(A); one <- ((A + t(A)) == 1) * 1
  trade <- (dc$Z + t(dc$Z)) / outer(dc$VA, dc$VA, "+")   # Frankel-Rose intensity
  g  <- dcast(gdp[year %in% win$years & iso3 %in% cc], year ~ iso3, value.var = "growth")
  ok <- intersect(cc, names(g)[colSums(!is.na(g)) >= length(win$years) - 1])
  G  <- as.matrix(g[, ok, with = FALSE])
  R  <- cor(G, use = "pairwise.complete.obs")
  idx <- which(upper.tri(R), arr.ind = TRUE)
  i <- match(ok[idx[, 1]], cc); j <- match(ok[idx[, 2]], cc)
  data.frame(ci = ok[idx[, 1]], cj = ok[idx[, 2]], corr = R[idx],
             circ = circ[cbind(i, j)] / sd(circ[upper.tri(circ)]),   # per s.d.
             ltrade = log(trade[cbind(i, j)]),
             reciprocal = rec[cbind(i, j)], oneway = one[cbind(i, j)])
}

fit_pairs <- function(f, d) {
  allc <- union(d$ci, d$cj)
  FE <- sapply(allc, function(k) as.numeric(d$ci == k | d$cj == k))  # country FE
  dd <- cbind(d, FE[, -1, drop = FALSE])
  ff <- update(f, paste(". ~ . +", paste0("`", colnames(FE)[-1], "`", collapse = " + ")))
  m <- lm(ff, data = dd)
  V <- vcovCL(m, cluster = dd[, c("ci", "cj")], multi0 = FALSE)   # two-way clustered
  # two-way clustered V need not be positive definite in the fixed-effect block:
  # standard errors are computed only for the coefficients reported
  keep <- setdiff(names(coef(m)), c("(Intercept)", paste0("`", colnames(FE)[-1], "`"), colnames(FE)[-1]))
  b <- coef(m); se <- setNames(rep(NA_real_, length(b)), names(b))
  se[keep] <- sqrt(diag(V)[keep])
  list(b = b, t = b / se, V = V, n = nrow(d), r2 = summary(m)$r.squared)
}

rows <- lapply(names(WINDOWS), function(nm) {
  d <- pair_data(WINDOWS[[nm]])
  m1 <- fit_pairs(corr ~ circ + ltrade, d)
  m2 <- fit_pairs(corr ~ reciprocal + oneway + ltrade, d)
  ct <- function(m, v) sprintf("%.3f (%.2f)", m$b[v], m$t[v])
  data.frame(Window = nm, Pairs = m1$n, Circuits = ct(m1, "circ"),
             Trade = ct(m1, "ltrade"), Reciprocal = ct(m2, "reciprocal"),
             Oneway = ct(m2, "oneway"),
             Diff = {L <- c(reciprocal = 1, oneway = -1)
                     dif <- sum(L * m2$b[names(L)])
                     sprintf("%.3f (%.2f)", dif, dif / sqrt(drop(L %*% m2$V[names(L), names(L)] %*% L)))})
})
tab <- do.call(rbind, rows)
write_tex_table(tab, "tabOA4_gdp_comovement.tex",
  caption = "Circuits and the co-movement of real GDP growth across country pairs",
  label = "tab:oa4", align = "lcccccc",
  header = c("Window", "Pairs", "Circuits (per s.d.)", "Log trade intensity",
             "Reciprocal link", "One-way link", "Rec.\\ $-$ one-way"),
  notes = paste("Dependent variable: correlation of annual real GDP growth (WDI) of the two",
    "countries over the window. Circuit intensity of the pair (weighted circuits of",
    "length 3 and 4, OECD ICIO) and links measured in the first year of the window;",
    "trade intensity $(Z_{ij}+Z_{ji})/(VA_i+VA_j)$. Country fixed effects; $t$-statistics",
    "clustered by both countries of the pair."))
print(tab)
log_msg("Appendix B.7 (GDP comovement) written: tabOA4_gdp_comovement")
}  # end if GDP file found
