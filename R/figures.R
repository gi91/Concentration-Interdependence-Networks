# =============================================================================
# figures.R -- figure builders shared by the two applications
# =============================================================================

# long format helper: vars -> labels
to_long <- function(df, x, vars, labs) do.call(rbind, Map(function(v, l)
  data.frame(x = df[[x]], v = df[[v]], s = l), vars, labs))

# (a) levels: Herfindahl, direct (psi_T^(1) or psi_D), one-way, reciprocal
fig_levels <- function(df, x, vars, labs, cols, ltys, title, subtitle = NULL,
                       nudge, crisis = FALSE, ylim = NULL) {
  L <- to_long(df, x, vars, labs); L$s <- factor(L$s, labs)
  ggplot(L, aes(x, v, colour = s, linetype = s)) +
    (if (crisis) geom_crisis() else NULL) +
    geom_line(linewidth = .65) +
    end_labels(L, "x", "v", "s", cols, nudge = nudge, box.padding = 0.55) +
    scale_colour_manual(values = setNames(cols, labs)) +
    scale_linetype_manual(values = setNames(ltys, labs)) +
    scale_x_continuous(breaks = data_breaks(L$x), expand = expansion(mult = c(.02, .32))) +
    coord_cartesian(ylim = ylim, clip = "off") +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_paper() + theme(legend.position = "none")
}

# (b) alignment factors with the 5-95% permutation band of reciprocal links
fig_alignment <- function(df, x, band, title, subtitle, nudge, crisis = FALSE,
                          log_y = FALSE) {
  labs <- c("Direct", "Reciprocal", "Circuits (3)", "Circuits (4)")
  cols <- c(INK, COL_MECH[["Reciprocal"]], COL_MECH[["Circuits"]], COL_MECH[["Circuits"]])
  L <- to_long(df, x, c("Gamma1", "G_rec", "G_c3", "G_c4"), labs)
  L$s <- factor(L$s, labs)
  ggplot(L, aes(x, v, colour = s, linetype = s)) +
    (if (crisis) geom_crisis() else NULL) +
    geom_ribbon(data = band, aes(x = x, ymin = lo, ymax = hi), inherit.aes = FALSE,
                fill = COL_MECH[["Reciprocal"]], alpha = .14) +
    geom_hline(yintercept = 1, colour = MUTED, linewidth = .3, linetype = "22") +
    geom_line(linewidth = .65) +
    end_labels(L, "x", "v", "s", cols, nudge = nudge, box.padding = 0.55) +
    scale_colour_manual(values = setNames(cols, labs)) +
    scale_linetype_manual(values = setNames(c("solid", "solid", "solid", "42"), labs)) +
    scale_x_continuous(breaks = data_breaks(L$x), expand = expansion(mult = c(.02, .32))) +
    (if (log_y) scale_y_log10(breaks = c(.25, .5, 1, 2, 4, 8, 16)) else NULL) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = if (log_y) "log scale" else NULL, title = title, subtitle = subtitle) +
    theme_paper() + theme(legend.position = "none")
}

# (c) Proposition 5: stacked effects per component and period, total marked
fig_decomp <- function(dec, title, subtitle = NULL) {
  dec$comp <- factor(ifelse(dec$comp == "OW", "One-way", "Reciprocal"),
                     c("One-way", "Reciprocal"))
  if (length(unique(dec$period)) > 3)   # long labels: "1995-07" -> "1995\n-07"
    dec$period <- sub("-", "\n\u2013", dec$period)
  dec$period <- factor(dec$period, unique(dec$period))
  L <- rbind(data.frame(comp = dec$comp, period = dec$period, v = dec$availability, s = "Availability"),
             data.frame(comp = dec$comp, period = dec$period, v = dec$weights, s = "Weights"),
             data.frame(comp = dec$comp, period = dec$period, v = dec$network, s = "Network"))
  L$s <- factor(L$s, names(COL_EFF))
  top <- pmax(0, dec$availability, 0) + pmax(dec$weights, 0) + pmax(dec$network, 0)
  rng <- diff(range(c(L$v, 0, top)))
  tot <- data.frame(comp = dec$comp, period = dec$period, v = dec$delta,
                    ylab = pmax(top, dec$delta) + 0.06 * rng,
                    lab = sprintf("%+.3f", dec$delta))
  ggplot(L, aes(period, v, fill = s)) +
    geom_hline(yintercept = 0, colour = MUTED, linewidth = .3) +
    geom_col(width = .62, colour = "white", linewidth = .5) +
    geom_point(data = tot, aes(period, v), inherit.aes = FALSE, shape = 23,
               size = 2.6, fill = "white", colour = INK, stroke = .8) +
    geom_text(data = tot, aes(period, ylab, label = lab), inherit.aes = FALSE,
              size = if (nlevels(dec$period) > 3) 1.9 else 2.3, colour = INK,
              fontface = "bold") +
    facet_wrap(~ comp, nrow = 1) +
    scale_fill_manual(values = COL_EFF) +
    labs(x = NULL, y = expression(Delta*psi), title = title, subtitle = subtitle,
         caption = "Diamond: total change; bars: availability, weight and network effects") +
    theme_paper() + theme(strip.text = element_text(colour = INK))
}

# pair contributions to Gamma of reciprocal links, by date
fig_pairs <- function(pc, title, subtitle, nrow = 1) {
  pc$date <- factor(pc$date, unique(pc$date))
  pc$pair <- tidytext_reorder(pc$pair, pc$contribution, pc$date)
  ggplot(pc, aes(contribution, pair)) +
    geom_col(fill = COL_MECH[["Reciprocal"]], width = .7) +
    geom_text(aes(label = sprintf("%.0f%%", 100 * share)), hjust = -.15, size = 2.2,
              colour = INK2) +
    facet_wrap(~ date, scales = "free_y", nrow = nrow) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_x_continuous(expand = expansion(mult = c(0, .25))) +
    labs(x = expression("Contribution to "*Gamma*" of reciprocal links"), y = NULL,
         title = title, subtitle = subtitle) + theme_paper()
}
# reorder within facets (avoids a dependency on tidytext)
tidytext_reorder <- function(x, by, within) {
  key <- paste(x, within, sep = "___")
  factor(key, levels = unique(key[order(within, by)]))
}

# threshold profile of Gamma (binary at several theta, weighted as reference)
fig_threshold <- function(tp, x, title, subtitle, crisis = FALSE, nudge) {
  tp$s <- factor(tp$s, unique(tp$s))
  n_bin <- length(unique(tp$s)) - 1
  cols <- c(colorRampPalette(c("#86b6ef", "#104281"))(n_bin), COL_MECH[["Reciprocal"]])
  names(cols) <- levels(tp$s)
  ggplot(tp, aes(x, v, colour = s)) +
    (if (crisis) geom_crisis(ymin = 0) else NULL) +
    geom_hline(yintercept = 1, colour = MUTED, linewidth = .3, linetype = "22") +
    geom_line(linewidth = .6) +
    end_labels(tp, "x", "v", "s", cols, nudge = nudge) +
    scale_colour_manual(values = cols) +
    scale_x_continuous(breaks = data_breaks(tp$x), expand = expansion(mult = c(.02, .25))) +
    scale_y_log10() + coord_cartesian(clip = "off") +
    labs(x = NULL, y = expression(Gamma*" of reciprocal links (log scale)"), title = title,
         subtitle = subtitle) + theme_paper() + theme(legend.position = "none")
}
