# =============================================================================
# utils.R  -- paths, plotting theme, figure saving, LaTeX (booktabs) tables
# =============================================================================

# UTF-8 is needed for the Greek letters and dashes used in figure labels
invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "C.UTF-8")))
if (!l10n_info()$`UTF-8`)
  invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "en_US.UTF-8")))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

DIR_FIG  <- file.path("output", "figures")
DIR_TAB  <- file.path("output", "tables")
DIR_PROC <- file.path("data", "processed")
DIR_RAW  <- file.path("data", "raw")
for (d in c(DIR_FIG, DIR_TAB, DIR_PROC)) dir.create(d, recursive = TRUE,
                                                     showWarnings = FALSE)

# ---- plotting ----------------------------------------------------------------
# Colour system: see R/plot_style.R (palette validated for colour-vision
# deficiency; every colour has a fixed meaning across all figures).
source(file.path("R", "plot_style.R"))

save_fig <- function(p, name, width = 6.5, height = 4.5) {
  dev <- if (capabilities("cairo")) cairo_pdf else pdf
  ggsave(file.path(DIR_FIG, paste0(name, ".pdf")), p, width = width,
         height = height, device = dev)
  ggsave(file.path(DIR_FIG, paste0(name, ".png")), p, width = width,
         height = height, dpi = 200)
  invisible(file.path(DIR_FIG, name))
}

# ---- LaTeX tables ------------------------------------------------------------
fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "--", formatC(x, format = "f", digits = digits))
}
fmt_coef_t <- function(b, t, digits = 2) {
  ifelse(is.na(b), "--",
         sprintf("%s (%s)", formatC(b, format = "f", digits = digits),
                 formatC(t, format = "f", digits = digits)))
}

# df: data.frame of character (already formatted) cells.
# panel: optional vector (same length as nrow(df)); a new panel title row is
# written each time the value changes.
write_tex_table <- function(df, file, caption = NULL, label = NULL,
                            notes = NULL, align = NULL, panel = NULL,
                            header = names(df), size = "\\footnotesize") {
  nc <- ncol(df)
  if (is.null(align)) align <- paste0("l", strrep("c", nc - 1))
  # escape bare % only (cells may already contain \\%), and typeset numeric
  # minus signs outside math as $-$
  esc <- function(x) {
    x <- gsub("(?<!\\\\)%", "\\\\%", x, perl = TRUE)
    ifelse(grepl("$", x, fixed = TRUE), x,
           gsub("(^|[\\s(\\[])-(?=[0-9.])", "\\1$-$", x, perl = TRUE))
  }
  L <- c("\\begin{table}[!htbp]", "\\centering", size)
  if (!is.null(caption)) L <- c(L, sprintf("\\caption{%s}", caption))
  if (!is.null(label))   L <- c(L, sprintf("\\label{%s}", label))
  L <- c(L, sprintf("\\begin{tabular}{%s}", align), "\\toprule",
         paste0(paste(header, collapse = " & "), " \\\\"), "\\midrule")
  prev <- NULL
  for (i in seq_len(nrow(df))) {
    if (!is.null(panel) && !identical(panel[i], prev)) {
      if (!is.null(prev)) L <- c(L, "\\midrule")
      L <- c(L, sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\", nc,
                        panel[i]))
      prev <- panel[i]
    }
    L <- c(L, paste0(paste(esc(as.character(unlist(df[i, ]))),
                           collapse = " & "), " \\\\"))
  }
  L <- c(L, "\\bottomrule", "\\end{tabular}")
  if (!is.null(notes))
    L <- c(L, sprintf("\\par\\smallskip\\parbox{\\linewidth}{\\scriptsize \\textit{Notes.} %s}", notes))
  L <- c(L, "\\end{table}")
  writeLines(L, file.path(DIR_TAB, file))
  utils::write.csv(df, file.path(DIR_TAB, sub("\\.tex$", ".csv", file)),
                   row.names = FALSE)
  invisible(L)
}

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n",
                             sep = "")
