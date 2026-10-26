# =============================================================================
# main.R -- reproduces all tables, figures and simulations of
#   Riso, L. and Zoia, M. G., "Concentration and changing interdependence in
#   production and banking networks"
#
# Run from the root of the repository:
#   RStudio:  open nci2.Rproj, then source("main.R")
#   Terminal: Rscript main.R
#
# Options (environment variables, set before running):
#   NCI_CORES=4      parallel chains in scripts 17-18 (Linux/macOS only)
#   NCI_QUICK=1      smoke test with few draws (minutes; numbers will differ)
#   NCI_REIMPORT=1   re-import data/raw even if data/processed is complete
# Without raw data in data/raw/, synthetic data are generated so that the code
# can be tested; the numbers are then not those of the paper.
# =============================================================================
if (Sys.getenv("NCI_QUICK") == "1")
  Sys.setenv(NCI_B = 100, NCI_B_FL = 300, NCI_DRAWS = 20, NCI_REPS = 60,
             NCI_DRAWS_SIM = 15, NCI_B_EX = 200, NCI_REPS_DESIGN = 15,
             NCI_B_DESIGN = 49, NCI_B_SECT = 50, NCI_B_FIXED = 100,
             NCI_REPS_E7 = 10, NCI_REPS_E8 = 10, NCI_B_BANKNULL = 99)
proc <- file.path("data", "processed")
have_raw  <- length(list.files(file.path("data", "raw"), "\\.(zip|csv|RData)$",
                               recursive = TRUE, ignore.case = TRUE)) > 0
have_proc <- file.exists(file.path(proc, "bis_claims.rds")) &&
             length(list.files(proc, "^icio_[0-9]{4}\\.rds$")) > 0
real <- have_raw || have_proc

steps <- c(
  if (have_raw && (!have_proc || Sys.getenv("NCI_REIMPORT") == "1")) "00_get_data.R",
  if (!real) "00_make_synthetic_data.R",
  "01_illustrative_example.R",   # Tables S1-S2, Figure S1
  "02_wiod.R",                   # Table 2, Figure 2, Tables S5-S6, Figures S4a, S5
  "03_bis.R",                    # Table 3, Figure 3, Tables S8-S9, Figures S4b, S7
  "04_null_tests.R",             # Table 4, Table S11, Figure S9
  "05_srisk_validation.R",       # Tables 5-6, Tables S4, S14
  "06_contagion.R",              # Figure 4
  "07_sim_propositions.R",       # Table S15
  "08_sim_robustness.R",         # Table S16 (E1-E6), Figures S10, S12
  "09_sim_design.R",             # Figure 1, Figures S2, S11, Table S17
  if (real) c(
    "11_lbs_robustness.R",       # Table S10, Figure S8
    "12_icio_sectors.R",         # Table S7, Figures S3, S6
    "13_gdp_comovement.R",       # Table S18 (needs data/raw/gdp)
    "14_descriptives.R",         # Table 1
    "16_bank_null_robustness.R", # Tables S12 (upper panel), S13
    "17_fixed_reciprocity.R"),   # Table 4 and S11 (fixed reciprocity), S12 (lower panel)
  "18_sim_fixed_reciprocity.R",  # Table S16 (E7-E8)
  "15_paper_assets.R")           # copies figures and tables to paper/

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
dir.create(file.path("output", "logs"), recursive = TRUE, showWarnings = FALSE)
t_all <- Sys.time()
for (s in steps) {
  cat("\n================", s, "================\n")
  t0 <- Sys.time()
  status <- system2(rscript, file.path("scripts", s))
  if (status != 0) stop("step failed: ", s)
  cat(sprintf("  -> %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
writeLines(capture.output(sessionInfo()), file.path("output", "logs", "sessionInfo.txt"))
cat(sprintf("\nAll steps completed in %.0f min. Figures: output/figures, tables: output/tables\n",
            as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
