# =============================================================================
# plot_style.R -- colour system and theme shared by every figure
#
# Each colour has ONE meaning in the whole paper:
#   mechanisms   one-way = blue, reciprocal = orange, open chains = aqua,
#                circuits = violet  (also used on the network drawings)
#   effects      availability = yellow, weights = magenta, network = green
#   direct / baseline concentration = near-black; marginal (Herfindahl) = grey
# Palettes validated with the dataviz validator (CVD-safe on adjacent pairs,
# light surface).  Colours below 3:1 contrast are always direct-labelled.
# =============================================================================
suppressPackageStartupMessages(library(ggrepel))

INK   <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#8a8984"
GRID  <- "#e8e7e3"; SHADE <- "#f0efec"; SURF <- "#fcfcfb"

COL_MECH <- c("One-way" = "#2a78d6", "Reciprocal" = "#eb6834",
              "Open chains" = "#1baf7a", "Circuits" = "#4a3aa7")
COL_EFF  <- c("Availability" = "#eda100", "Weights" = "#e87ba4",
              "Network" = "#008300")
COL_NET  <- c("Network A" = "#0b0b0b", "Network B" = "#e34948")   # not mechanism colours
COL_BASE <- INK; COL_HHI <- MUTED
SEQ3     <- c("#86b6ef", "#2a78d6", "#104281")      # ordinal blue ramp
CAT      <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
              "#008300", "#4a3aa7", "#e34948")       # fixed categorical order

theme_paper <- function(base_size = 9) {
  theme_minimal(base_size = base_size) +
    theme(plot.title = element_text(face = "bold", size = base_size + 0.5,
                                    hjust = 0, colour = INK),
          plot.subtitle = element_text(size = base_size - 1, colour = INK2,
                                       hjust = 0, margin = margin(b = 4)),
          plot.caption = element_text(size = base_size - 2, colour = MUTED,
                                      hjust = 0),
          plot.title.position = "plot", plot.caption.position = "plot",
          axis.title = element_text(size = base_size - 1, colour = INK2),
          axis.text = element_text(size = base_size - 1.5, colour = INK2),
          panel.grid.major = element_line(colour = GRID, linewidth = 0.25),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(colour = MUTED, linewidth = 0.3),
          legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = base_size - 1.5, colour = INK2),
          legend.key.width = unit(1.3, "lines"),
          legend.margin = margin(0, 0, 0, 0),
          strip.text = element_text(face = "bold", colour = INK, hjust = 0),
          plot.background = element_rect(fill = "white", colour = NA))
}

# crisis shading for the banking figures, with labels
# time is plotted as decimal years (2007-Q1 = 2007.00, 2007-Q3 = 2007.50)
dec_year <- function(d) as.numeric(format(d, "%Y")) +
  (as.numeric(format(d, "%m")) - 1) / 12
crisis_rects <- data.frame(xmin = c(2007.5, 2010.25), xmax = c(2009.5, 2013.0),
                           lab = c("GFC", "EA crisis"))

geom_crisis <- function(label_y = Inf, ymin = -Inf) list(
  geom_rect(data = transform(crisis_rects, ymin = ymin),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = Inf),
            inherit.aes = FALSE, fill = SHADE),
  geom_text(data = crisis_rects, aes(x = xmin + (xmax - xmin) / 2, y = label_y,
            label = lab), inherit.aes = FALSE, vjust = 1.4, size = 2.1,
            colour = MUTED))

# Direct labels at the right end of each line (identity never colour-only)
end_labels <- function(df, x, y, group, colours, nudge = 0.02, size = 2.5,
                       parse = FALSE, box.padding = 0.3) {
  last <- do.call(rbind, lapply(split(df, df[[group]]), function(d)
    d[which.max(d[[x]]), ]))
  geom_text_repel(data = last, aes(x = .data[[x]], y = .data[[y]],
                  label = .data[[group]], colour = .data[[group]]),
                  inherit.aes = FALSE, direction = "y", hjust = 0, nudge_x = nudge,
                  segment.colour = GRID, size = size, show.legend = FALSE,
                  min.segment.length = 0, box.padding = box.padding, force = 2, parse = parse,
                  xlim = c(NA, Inf), seed = 1)
}

# reference line Gamma = 1 (expected alignment under random weights)
ref_gamma1 <- function(x = -Inf) list(
  geom_hline(yintercept = 1, colour = MUTED, linewidth = 0.3, linetype = "22"),
  annotate("text", x = x, y = 1, label = "E[Gamma] == 1~'(random weights)'",
           parse = TRUE, hjust = -0.05, vjust = -0.5, size = 2.2, colour = MUTED))

# axis breaks computed on the data range (not on the space left for labels)
data_breaks <- function(x, n = 5) {
  b <- pretty(range(x, na.rm = TRUE), n = n)
  b[b >= min(x, na.rm = TRUE) - 1e-9 & b <= max(x, na.rm = TRUE) + 1e-9]
}
