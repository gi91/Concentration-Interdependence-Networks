# =============================================================================
# 00_get_data.R -- download / import the real data and build data/processed/
# Run on your own computer from the project root:  Rscript scripts/00_get_data.R
#
#   BIS CBS   bulk file WS_CBS_PUB_csv_col.zip into data/raw/bis/ (or downloaded
#             automatically with the CRAN package "BIS" if the zip is absent)
#   BIS LBS   optional, robustness: WS_LBS_D_PUB_csv_col.zip into data/raw/bis/
#   OECD ICIO zip files from the OECD ICIO page into data/raw/icio/, as they
#             are (e.g. 2016-2022_SML.zip): no need to unzip
#   WIOD 2016 optional (robustness 2000-2014): WIOT<year>_October16_ROW.RData
#             into data/raw/wiod/
#   SRISK     V-Lab (free login): one export PER COUNTRY ("<country> Financials -
#             Total SRISK"), saved as data/raw/srisk/<ISO2>.csv (e.g. US.csv).
#             The world aggregate is not usable.
#
# Each block is skipped if its raw files are missing, and prints diagnostics
# to check the layout before the analysis scripts are run.
# =============================================================================
source("R/utils.R"); source("R/data_loaders.R")
for (d in c("bis", "icio", "wiod", "srisk")) dir.create(file.path(DIR_RAW, d),
                                                        recursive = TRUE, showWarnings = FALSE)
REPORTERS <- c("AU","AT","BE","CL","FR","DE","GR","IE","IT","JP","NL","PT","ES",
               "SE","CH","TW","TR","GB","US")

# ---- 1. BIS consolidated banking statistics ----------------------------------------
# Foreign claims, ultimate-risk basis, domestic banks, all sectors, all
# instruments, amounts outstanding.  Verified on WS_CBS_PUB_csv_col.zip (BIS
# release of 14 Sept 2026): this filter reproduces Table 3 of the paper exactly.
BIS_FILTER <- list(
  L_MEASURE      = "S",    # amounts outstanding / stocks
  CBS_BANK_TYPE  = "4R",   # domestic banks, excluding domestic positions
  CBS_BASIS      = "U",    # guarantor = ultimate-risk basis ("F" = immediate counterparty)
  L_POSITION     = "C",    # total claims
  L_INSTR        = "A",    # all instruments
  REM_MATURITY   = "A",    # all maturities
  CURR_TYPE_BOOK = "TO1",  # all currencies
  L_CP_SECTOR    = "A")    # all counterparty sectors
TOTAL_CP <- "5J"           # counterparty "all countries"

cbs_zip <- list.files(file.path(DIR_RAW, "bis"), "^WS_CBS_PUB.*\\.zip$", full.names = TRUE)
if (length(cbs_zip)) {
  # (a) bulk file downloaded by hand (data.bis.org > Bulk downloads, CSV "csv_col")
  long <- read_bis_wide(cbs_zip[1], BIS_FILTER, from_year = 2006)
  log_msg("CBS read from ", basename(cbs_zip[1]))
} else {
  # (b) otherwise download with the BIS package (long format)
  if (!requireNamespace("BIS", quietly = TRUE)) install.packages("BIS")
  ds <- BIS::get_datasets()
  url <- ds$url[grepl("consolidated", ds$name, ignore.case = TRUE)][1]
  log_msg("downloading BIS CBS from ", url)
  long <- as.data.table(BIS::get_bis(url)); setnames(long, toupper(names(long)))
  for (v in intersect(names(BIS_FILTER), names(long))) long <- long[get(v) == BIS_FILTER[[v]]]
  long <- long[grepl("^[0-9]{4}-Q[1-4]$", TIME_PERIOD)]
  long[, `:=`(value = as.numeric(OBS_VALUE), quarter = as.Date(paste0(
    substr(TIME_PERIOD, 1, 4), "-", (as.integer(substr(TIME_PERIOD, 7, 7)) - 1) * 3 + 1, "-01")))]
}
cbs <- bis_standardize(long, REPORTERS, TOTAL_CP)
claims <- cbs$claims; totals <- cbs$totals
cat("CBS reporters in the file:", paste(sort(unique(long$L_REP_CTY)), collapse = " "), "\n")
# keep the quarters in which all 19 reporters report a total
ok_q <- totals[!is.na(total_claims), .N, by = quarter][N == length(REPORTERS), quarter]
claims <- claims[quarter %in% ok_q & !is.na(claims)]; totals <- totals[quarter %in% ok_q]
sh <- claims[, .(b = sum(claims)), by = .(quarter, creditor)][totals, on = .(quarter, creditor)]
cat(sprintf("CBS: %d quarters with all %d reporters (%s to %s); claims on the reporters as a share of total foreign claims: %.2f-%.2f (paper: 0.62-0.75)\n",
            length(ok_q), length(REPORTERS), min(ok_q), max(ok_q),
            min(sh[, sum(b, na.rm = TRUE) / sum(total_claims), by = quarter]$V1),
            max(sh[, sum(b, na.rm = TRUE) / sum(total_claims), by = quarter]$V1)))
# missing bilateral cells (confidential): treated as no link, as in the paper
miss <- claims[quarter == max(ok_q), .(n_cp = .N), by = creditor][n_cp < length(REPORTERS) - 1]
if (nrow(miss)) cat("  reporters with unpublished bilateral cells (last quarter):",
                    paste(miss$creditor, miss$n_cp, sep = ":", collapse = " "), "\n")
saveRDS(claims, file.path(DIR_PROC, "bis_claims.rds"))
saveRDS(totals, file.path(DIR_PROC, "bis_totals.rds"))
# extended sample with India (row 10 of Table 5): 20 reporters
ext <- bis_standardize(long, c(REPORTERS, "IN"), TOTAL_CP)
saveRDS(ext$claims[!is.na(claims)], file.path(DIR_PROC, "bis_claims_ext.rds"))
saveRDS(ext$totals[!is.na(total_claims)], file.path(DIR_PROC, "bis_totals_ext.rds"))

# ---- 1b. BIS locational banking statistics (robustness only) -------------------
# Cross-border total claims of all reporting banks resident in each country,
# all instruments, currencies and sectors.  Residence basis, immediate
# counterparty: NOT the data of the main analysis.
LBS_FILTER <- list(L_MEASURE = "S", L_POSITION = "C", L_INSTR = "A", L_DENOM = "TO1",
                   L_CURR_TYPE = "A", L_PARENT_CTY = "5J", L_REP_BANK_TYPE = "A",
                   L_CP_SECTOR = "A", L_POS_TYPE = "N")
lbs_zip <- list.files(file.path(DIR_RAW, "bis"), "^WS_LBS_D_PUB.*\\.zip$", full.names = TRUE)
if (length(lbs_zip)) {
  lbs <- bis_standardize(read_bis_wide(lbs_zip[1], LBS_FILTER), REPORTERS)
  cov <- lbs$claims[, .(n_cp = uniqueN(counterparty)), by = creditor]
  cat("LBS: reporters with bilateral claims on the other reporters:",
      paste(sort(cov$creditor), collapse = " "), "\n  missing:",
      paste(setdiff(REPORTERS, cov$creditor), collapse = " "), "\n")
  saveRDS(lbs$claims, file.path(DIR_PROC, "lbs_claims.rds"))
  saveRDS(lbs$totals, file.path(DIR_PROC, "lbs_totals.rds"))
}

# ---- 2. OECD ICIO -------------------------------------------------------------------
# zip files are read one CSV at a time (no need to unzip by hand)
zips <- list.files(file.path(DIR_RAW, "icio"), "\\.zip$", full.names = TRUE)
items <- do.call(rbind, lapply(zips, function(z)
  data.frame(zip = z, csv = unzip(z, list = TRUE)$Name)))
csvs <- list.files(file.path(DIR_RAW, "icio"), "\\.csv$", full.names = TRUE, recursive = TRUE)
items <- rbind(items, data.frame(zip = rep(NA, length(csvs)), csv = csvs))
items <- items[grepl("\\.csv$", items$csv), ]
for (k in seq_len(nrow(items))) {
  f <- items$csv[k]
  yr <- as.integer(regmatches(basename(f), regexpr("(19|20)[0-9]{2}", basename(f))))
  if (length(yr) == 0) next
  out <- file.path(DIR_PROC, sprintf("icio_%d.rds", yr))
  if (file.exists(out)) next
  if (!is.na(items$zip[k])) {
    tmp <- tempdir(); unzip(items$zip[k], files = f, exdir = tmp)
    path <- file.path(tmp, f)
  } else path <- f
  d <- read_icio_csv(path)
  if (!is.na(items$zip[k])) unlink(path)
  saveRDS(d, out)
  log_msg(sprintf("ICIO %d: %d countries x %d industries, VA share of output %.2f",
                  yr, length(unique(d$country)), length(unique(d$industry)),
                  sum(d$va) / sum(d$go)))
}

# ---- 3. WIOD 2016 release (optional) ---------------------------------------------------
for (f in list.files(file.path(DIR_RAW, "wiod"), "\\.RData$", full.names = TRUE)) {
  yr <- as.integer(regmatches(basename(f), regexpr("(19|20)[0-9]{2}", basename(f))))
  out <- file.path(DIR_PROC, sprintf("wiod_%d.rds", yr))
  if (!file.exists(out)) { saveRDS(read_wiod_rdata(f), out); log_msg("WIOD ", yr) }
}

# ---- 4. SRISK (V-Lab exports, one CSV per country) --------------------------------------
# V-Lab "custom" export: line 1 = ',"<Country> ... - Total SRISK (US$ billion)"',
# then date,value.  The country is read from the header, so files can keep the
# names given by V-Lab.  Aggregates ("World", "Selected Countries") are rejected.
SRISK_COUNTRIES <- c("Australia" = "AU", "Austria" = "AT", "Belgium" = "BE", "Chile" = "CL",
  "France" = "FR", "Germany" = "DE", "Greece" = "GR", "Ireland" = "IE", "Italy" = "IT",
  "Japan" = "JP", "Netherlands" = "NL", "Portugal" = "PT", "Spain" = "ES", "Sweden" = "SE",
  "Switzerland" = "CH", "Taiwan" = "TW", "Turkey" = "TR", "Turkiye" = "TR",
  "United Kingdom" = "GB", "United States" = "US", "India" = "IN")
# CSV files, also inside zip archives (e.g. V-Lab.zip as downloaded), read directly
sdir <- file.path(DIR_RAW, "srisk"); stmp <- file.path(tempdir(), "srisk_unzipped")
for (z in list.files(sdir, "\\.zip$", full.names = TRUE, recursive = TRUE))
  utils::unzip(z, exdir = stmp)
sf <- c(list.files(sdir, "\\.csv$", full.names = TRUE, recursive = TRUE),
        list.files(stmp, "\\.csv$", full.names = TRUE, recursive = TRUE))
if (length(sf)) {
  sr <- rbindlist(lapply(sf, function(f) {
    head1 <- readLines(f, n = 1, warn = FALSE)
    if (grepl("World|Selected Countries|Region", head1, ignore.case = TRUE)) {
      warning(basename(f), ": aggregate series (", head1, ") - skipped"); return(NULL) }
    hit <- names(SRISK_COUNTRIES)[vapply(names(SRISK_COUNTRIES),
             function(n) grepl(n, head1, ignore.case = TRUE), logical(1))]
    iso <- if (length(hit)) SRISK_COUNTRIES[[hit[which.max(nchar(hit))]]]
           else sub("\\.csv$", "", basename(f))          # fall back to <ISO2>.csv
    unit <- if (grepl("million", head1, ignore.case = TRUE)) 1000 else 1
    # some exports add "SRISK Capacity" / "Probability of Crisis": keep column 2
    x <- fread(f, header = TRUE, select = 1:2, col.names = c("date", "srisk"))
    # V-Lab leaves the cell blank when SRISK <= 0 (no exact zeros occur in the
    # exports): blanks are kept as NA here and read as zero shortfall downstream
    data.table(date = as.Date(x$date), country = iso,
               srisk = pmax(as.numeric(x$srisk) / unit, 0))
  }), fill = TRUE)
  sr[, date := as.Date(format(date, "%Y-%m-01"))]
  sr <- sr[, .(srisk = srisk[.N]), by = .(date, country)]
  cov <- sr[!is.na(srisk), .(from = min(date), to = max(date)), by = country]
  print(cov[order(country)])
  cat("missing SRISK countries:", paste(setdiff(unname(unique(SRISK_COUNTRIES)), cov$country),
                                         collapse = " "), "\n")
  saveRDS(sr, file.path(DIR_PROC, "srisk.rds"))
}
log_msg("done: now run  NCI_PROD=icio Rscript main.R  (or NCI_PROD=wiod)")
