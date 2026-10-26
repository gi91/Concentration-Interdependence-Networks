# =============================================================================
# 14_descriptives.R -- Section 4.4: descriptive statistics of the networks
# Output: tab_data_descriptives.tex/csv
# Uses the results saved by 02, 03, 12 and the processed SRISK file.
# =============================================================================
source("R/utils.R")
rng <- function(x, d = 3) sprintf("%s--%s", fmt_num(min(x, na.rm = TRUE), d),
                                    fmt_num(max(x, na.rm = TRUE), d))
qlab <- function(d) paste0(format(d, "%Y"), "-Q", (as.integer(format(d, "%m")) - 1) %/% 3 + 1)

rows <- list()
pc <- readRDS(file.path(DIR_PROC, "results_icio.rds"))$ts
rows[[1]] <- data.frame(Network = "Production, countries", Source = "OECD ICIO 2025",
  Nodes = "80 economies", Period = sprintf("%d--%d, annual", min(pc$year), max(pc$year)),
  Weights = "Value added", Link = "$Z_{ji}/X_i \\ge 0.5\\%$",
  Density = rng(pc$density), HHI = rng(pc$hhi))
f12 <- file.path(DIR_PROC, "results_icio_sectors.rds")
if (file.exists(f12)) {
  ps <- readRDS(f12)$ts
  for (nt in unique(ps$net)) { r <- ps[ps$net == nt, ]
    rows[[length(rows) + 1]] <- data.frame(Network = paste0("Production, industries (", tolower(nt), ")"),
      Source = "OECD ICIO 2025", Nodes = "50 industries",
      Period = sprintf("%d--%d, annual", min(r$year), max(r$year)), Weights = "Value added",
      Link = if (nt == "All flows") "$Z_{ji}/X_i \\ge 2\\%$" else "same density",
      Density = rng(r$density), HHI = rng(r$hhi)) }
}
bb <- readRDS(file.path(DIR_PROC, "results_bis.rds"))$ts
rows[[length(rows) + 1]] <- data.frame(Network = "Banking, consolidated", Source = "BIS CBS",
  Nodes = "19 banking systems", Period = sprintf("%s--%s, quarterly", qlab(min(bb$quarter)), qlab(max(bb$quarter))),
  Weights = "Foreign claims", Link = "claims share $\\ge 2\\%$",
  Density = rng(bb$density), HHI = rng(bb$hhi))
f11 <- file.path(DIR_PROC, "results_lbs.rds")
if (file.exists(f11)) { l <- readRDS(f11); l <- l[l$source == "LBS, residence", ]
  rows[[length(rows) + 1]] <- data.frame(Network = "Banking, locational (robustness)", Source = "BIS LBS",
    Nodes = "17 reporting countries", Period = sprintf("%s--%s, quarterly",
      qlab(min(as.Date(l$quarter, origin = "1970-01-01"))), qlab(max(as.Date(l$quarter, origin = "1970-01-01")))),
    Weights = "Cross-border claims", Link = "claims share $\\ge 2\\%$",
    Density = rng(l$density), HHI = rng(l$hhi)) }
tab <- do.call(rbind, rows)
sr <- readRDS(file.path(DIR_PROC, "srisk.rds"))
write_tex_table(tab, "tab_data_descriptives.tex",
  caption = "Networks used in the applications", label = "tab:data",
  align = "llllllcc",
  header = c("Network", "Source", "Nodes", "Period", "Weights", "Link", "Density", "HHI of $w$"),
  notes = paste0("Ranges over the sample period. Systemic risk: SRISK (V-Lab, NYU Stern), monthly, ",
    length(unique(sr$country)), " countries, ", format(min(sr$date), "%Y"), "--",
    format(max(sr$date), "%Y"), ", US\\$ billion."))
print(tab[, c("Network", "Nodes", "Period", "Density", "HHI")])
log_msg("descriptive table written")
