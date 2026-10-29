# Concentration and changing interdependence in production and banking networks

Replication code for

>  Quatto P., Riso L., and Zoia, M. G. *Concentration and changing interdependence in production and banking networks.*

The code reproduces every table, figure and simulation of the paper and of its Supplementary Material, from the raw public data to the final numbers. It is written in R. One sampler is written in C++ and called from R through Rcpp.

## Requirements

- R ≥ 4.3 (tested with R 4.3.3)
- R packages: `data.table`, `ggplot2`, `patchwork`, `ggrepel`, `scales`, `igraph`, `sandwich`, `Rcpp`

  ```r
  install.packages(c("data.table", "ggplot2", "patchwork", "ggrepel", "scales",
                     "igraph", "sandwich", "Rcpp"))
  ```
- A C++ compiler for Rcpp, used by the samplers with fixed reciprocity in `src/constrained_mcmc.cpp`:
  - Windows: [Rtools](https://cran.r-project.org/bin/windows/Rtools/), the version that matches your R;
  - macOS: `xcode-select --install`;
  - Linux: `g++`, usually already installed.

## Data

The raw data are public, but they are not redistributed here. Download them and place the files, unzipped or as downloaded, in the folders below.

| Folder | Source | Files |
|---|---|---|
| `data/raw/icio/` | [OECD Inter-Country Input-Output tables, 2025 edition](https://www.oecd.org/en/data/datasets/inter-country-input-output-tables.html), SML format | `1995-2000_SML.zip`, `2001-2005_SML.zip`, `2006-2010_SML.zip`, `2011-2015_SML.zip`, `2016-2022_SML.zip` |
| `data/raw/bis/` | [BIS data portal, bulk download](https://data.bis.org/bulkdownload), CSV (`csv_col`) | `WS_CBS_PUB_csv_col.zip` (consolidated banking statistics); `WS_LBS_D_PUB_csv_col.zip` (locational banking statistics, Table S10) |
| `data/raw/srisk/` | [V-Lab, NYU Stern, country SRISK](https://vlab.stern.nyu.edu/srisk) (free registration) | one CSV export per country, e.g. `vlab-srisk-usa-YYYYMMDD.csv` |
| `data/raw/gdp/` | [World Bank WDI, GDP growth (annual %)](https://data.worldbank.org/indicator/NY.GDP.MKTP.KD.ZG), CSV | `API_NY.GDP.MKTP.KD.ZG_DS2_en_csv_v2_*.zip` (Table S18 only) |

The paper uses the files downloaded in October 2026: BIS release of 14 September 2026 and V-Lab exports of 8 October 2026. Later releases revise past observations, so the numbers can differ slightly. Script `00_get_data.R` reads the zip files directly, builds `data/processed/`, and prints diagnostics on the layout of the files.

## How to run

From the root of the repository:

```r
source("main.R")          # in RStudio, after opening nci2.Rproj
```

```bash
Rscript main.R            # from a terminal
NCI_CORES=4 Rscript main.R   # parallel chains in scripts 17-18 (Linux/macOS)
NCI_QUICK=1 Rscript main.R   # smoke test with few draws (a few minutes)
```

`main.R` runs the scripts in `scripts/` in order, each in a fresh R session, and stops at the first error. A complete run takes two to three hours on a laptop. The longest steps are `09_sim_design.R` (about 35 minutes), `18_sim_fixed_reciprocity.R` (about 25 minutes on two cores, about 50 minutes on Windows, where the chains run one after the other) and `02_wiod.R` (about 15 minutes). Every script fixes its random seeds, so repeated runs give identical numbers.

Without files in `data/raw/`, `main.R` generates synthetic data. This lets you test the code, but the numbers are then not those of the paper.

Outputs: tables in `output/tables/` (`.tex` and `.csv`), figures in `output/figures/` (`.pdf` and `.png`), and the R session information in `output/logs/sessionInfo.txt`.

## Where each result comes from

### Main text

| Paper | Script | Output file |
|---|---|---|
| Table 1 | `14_descriptives.R` | `tab_data_descriptives` |
| Table 2, Figure 2 | `02_wiod.R` | `tab2_icio`, `fig2_icio` |
| Table 3, Figure 3 | `03_bis.R` | `tab3_bis`, `fig3_bis` |
| Table 4 | `04_null_tests.R`, `17_fixed_reciprocity.R` | `tab4_null_models` |
| Tables 5–6 | `05_srisk_validation.R` | `tab5_srisk`, `tab6_centrality` |
| Figure 1 | `09_sim_design.R` | `figS4_size_topology` |
| Figure 4 | `06_contagion.R` | `fig4_contagion` |

### Supplementary Material

| Paper | Script | Output file |
|---|---|---|
| Tables S1–S2, Figure S1 | `01_illustrative_example.R` | `tab1_example`, `tab1b_example_dynamic`, `fig1_example` |
| Table S3, Table S17, Figures S2, S11 | `09_sim_design.R` | `simS3_design`, `tabS3_design`, `figS6_comparability`, `figS5_power_topology` |
| Figure S3, Table S7, Figure S6 | `12_icio_sectors.R` | `fig6_levels`, `tabOA_sectors`, `fig6b_pairs_sectors` |
| Figure S4 | `02_wiod.R`, `03_bis.R` | `fig5a_pairs_icio`, `fig5b_pairs_bis` |
| Table S4, Table S14 | `05_srisk_validation.R` | `tab5_srisk`, `tabOA3_weighted_embeddedness` |
| Tables S5–S6, Figure S5 | `02_wiod.R` | `tabOA1_icio_robustness`, `tabOA_icio_within_order`, `figOA_threshold_icio` |
| Tables S8–S9, Figure S7 | `03_bis.R` | `tabOA_bis_threshold`, `tabOA_bis_within_order`, `figOA_threshold_bis` |
| Table S10, Figure S8 | `11_lbs_robustness.R` | `tabOA_lbs`, `figOA_lbs` |
| Table S11, Figure S9 | `04_null_tests.R`, `17_fixed_reciprocity.R` | `tabOA2_null_models_full.csv`, `figOA_power_icio2022` |
| Table S12 | `16_bank_null_robustness.R` (upper panel), `17_fixed_reciprocity.R` (lower panel) | `tabOA_bank_null_robustness`, `tabOA_bank_null_fixed.csv` |
| Table S13 | `16_bank_null_robustness.R` | `tabOA_bank_missing` |
| Table S15 | `07_sim_propositions.R` | `tabS1_proposition_checks` |
| Table S16, Figures S10, S12 | `08_sim_robustness.R` (E1–E6), `18_sim_fixed_reciprocity.R` (E7–E8) | `tabS2_robustness_summary`, `figS1_power_permutation`, `figS3_shares_vs_N` |
| Table S18 | `13_gdp_comovement.R` | `tabOA4_gdp_comovement` |

`15_paper_assets.R` copies the figures and the bodies of the tables to `paper/`.

## Code

| File | Content |
|---|---|
| `R/core_indices.R` | quadratic-form concentration index, availability–alignment factorisation, distance and propagation decompositions, structural decomposition (one-way, reciprocal, open chains, circuits), decomposition of changes, embeddedness |
| `R/null_models.R` | permutation test and analytic permutation moments; exact degree-preserving rewiring (Definition S1(a)) |
| `R/fixes.R` | directed binary configuration model (Definition S1(b)), pair contributions, within-order shares, excess embeddedness |
| `R/constrained_nulls.R`, `src/constrained_mcmc.cpp` | samplers with fixed degrees and fixed reciprocal pairs or partners (Definition S1(c)–(d)); four chains, 500\|E\| proposals of burn-in and 100\|E\| between draws |
| `R/contagion.R` | linear DebtRank and the reciprocity-tilted Markov chain of the contagion experiment |
| `R/validation.R` | robust OLS, Freedman–Lane permutation, leave-one-out, centralities |
| `R/data_loaders.R` | readers for the OECD, BIS, V-Lab and World Bank files |
| `R/utils.R`, `R/plot_style.R`, `R/figures.R` | paths, table writer, figure style |

Script `17_fixed_reciprocity.R` first checks the samplers by complete enumeration of all directed networks with five nodes. Every draw must lie in the constrained set, and the networks it reaches must be visited with equal frequency. The results are in `fixed_recip_check_enumeration.csv`.

## License

Code released under the MIT License (see `LICENSE`). The data remain subject to the terms of their providers.

## Citation

If you use this code, please cite the paper (see `CITATION.cff`).
