# =============================================================================
# 01_illustrative_example.R -- Section 6: Figure 1 and Table 1
# =============================================================================
source("R/utils.R"); source("R/core_indices.R")

N <- 8
mk <- function(edges) { A <- matrix(0, N, N)
  for (e in edges) A[e[2], e[1]] <- 1; A }            # A[i, j]: j -> i
nets <- list(
  A = mk(list(c(1,3), c(3,2), c(2,4), c(4,1), c(4,5), c(5,6), c(6,7),
              c(7,8), c(8,7))),
  B = mk(list(c(4,1), c(1,3), c(3,5), c(5,6), c(6,5), c(7,6), c(2,7),
              c(8,7), c(8,2))))
w  <- c(0.3, 0.3, rep(0.4 / 6, 6))
w1 <- c(0.36, 0.24, rep(0.4 / 6, 6))                  # dynamic illustration

# ---- all measures ------------------------------------------------------------
res <- lapply(nets, function(A) {
  b  <- baseline_index(w, A)
  tp <- topological_profile(w, A)
  p1 <- propagation_measures(w, A, 1); p3 <- propagation_measures(w, A, 3)
  sd <- structural_decomposition(w, A, H = 3)
  tb <- sd$table
  list(b = b, tp = tp, p1 = p1, p3 = p3, sd = sd,
       G_open2 = tb$Gamma[tb$comp == "R" & tb$h == 2],
       G_circ3 = tb$Gamma[tb$comp == "F" & tb$h == 3],
       psiO = sapply(1:6, function(k) propagation_measures(w, A, k)$psiO))
})

# ---- Table 1 -------------------------------------------------------------------
f3 <- function(x) fmt_num(x, 3); f2 <- function(x) fmt_num(x, 2)
row <- function(dim, lab, fa, fb) data.frame(Dimension = dim, Measure = lab,
                                             A = fa, B = fb)
r <- res
tab1 <- rbind(
  row("Marginal", "Herfindahl index of $w$", f3(hhi(w)), f3(hhi(w))),
  row("Direct", "$\\psi_D(w)$; $\\bar m_1$; $\\Gamma_1(w)$",
      paste(f3(r$A$b$psiD), f3(r$A$b$density), f3(r$A$b$Gamma1), sep = "; "),
      paste(f3(r$B$b$psiD), f3(r$B$b$density), f3(r$B$b$Gamma1), sep = "; ")),
  row("Topological", "$\\psi^{(2)}_G(w)$; $\\Gamma_2(w)$",
      paste(f3(r$A$tp$profile$psi[2]), f3(r$A$tp$profile$Gamma[2]), sep = "; "),
      paste(f3(r$B$tp$profile$psi[2]), f3(r$B$tp$profile$Gamma[2]), sep = "; ")),
  row("", "Economic mean distance $\\bar d_w$; average path length $\\bar d$",
      paste(f2(r$A$tp$dw), f2(r$A$tp$dbar), sep = "; "),
      paste(f2(r$B$tp$dw), f2(r$B$tp$dbar), sep = "; ")),
  row("Propagation", "$\\psi^{(1)}_O(w)$; $\\psi^{(3)}_O(w)$",
      paste(f3(r$A$p1$psiO), f3(r$A$p3$psiO), sep = "; "),
      paste(f3(r$B$p1$psiO), f3(r$B$p3$psiO), sep = "; ")),
  row("", "Herfindahl index of $w^{(3)} = P^3 w$", f3(r$A$p3$hhi_wk),
      f3(r$B$p3$hhi_wk)),
  row("Structural", "$\\psi_{TOT}(w)$", f3(r$A$sd$psiTOT), f3(r$B$sd$psiTOT)),
  row("", "Shares: one-way; reciprocal",
      paste(f2(r$A$sd$shares[1]), f2(r$A$sd$shares[2]), sep = "; "),
      paste(f2(r$B$sd$shares[1]), f2(r$B$sd$shares[2]), sep = "; ")),
  row("", "Shares: open chains; circuits",
      paste(f2(r$A$sd$shares[3]), f2(r$A$sd$shares[4]), sep = "; "),
      paste(f2(r$B$sd$shares[3]), f2(r$B$sd$shares[4]), sep = "; ")),
  row("", "$\\Gamma$ of open chains ($h=2$); of circuits ($h=3$)",
      paste(f2(r$A$G_open2), f2(r$A$G_circ3), sep = "; "),
      paste(f2(r$B$G_open2), f2(r$B$G_circ3), sep = "; ")))
write_tex_table(tab1, "tab1_example.tex",
  caption = "Concentration measures for the two networks of the illustrative example",
  label = "tab:example", align = "llcc",
  header = c("Dimension", "Measure", "Network A", "Network B"),
  notes = paste("$N=8$, $w_1=w_2=0.3$ and $w_i=1/15$ for $i=3,\\dots,8$.",
    "Topological indices use the binary specification. Structural indices",
    "use $H=3$ and $\\alpha_h=1/3$; shares are fractions of $\\psi_{TOT}$",
    "(Proposition 4(ii)). The circuit component of network B is empty."))

# ---- dynamic illustration (end of Section 6) -----------------------------------
SA <- support_undirected(nets$A); SB <- support_undirected(nets$B)
om <- function(A) order_matrices(A, "binary")$M[[2]]
dyn <- rbind(
  cbind(order = "first",  dynamic_decomposition(w, SB, w1, SA)),
  cbind(order = "second", dynamic_decomposition(w, om(nets$B), w1, om(nets$A))))
tabd <- data.frame(Order = c("First ($\\psi_D$)", "Second ($\\psi^{(2)}_G$)"),
  t0 = f3(dyn$psi0), t1 = f3(dyn$psi1), d = f3(dyn$delta),
  a = f3(dyn$availability), wgt = f3(dyn$weights), net = f3(dyn$network),
  shw = f3(dyn$shapley_weights), shn = f3(dyn$shapley_network))
write_tex_table(tabd, "tab1b_example_dynamic.tex",
  caption = "Illustrative example: decomposition of the change from network B ($t-1$) to network A ($t$)",
  label = "tab:example-dyn", align = "lcccccccc",
  header = c("Index", "$\\psi_{t-1}$", "$\\psi_t$", "$\\Delta\\psi_t$",
             "$\\Delta^{(a)}$", "$\\Delta^{(w)}$", "$\\Delta^{(M)}$",
             "Shapley $w$", "Shapley $M$"),
  notes = paste("Columns 5--7: Proposition 5. Columns 8--9: two-player Shapley",
    "decomposition of $\\psi$ over the primitive changes $(w, M)$, in which the",
    "network effect includes the availability effect."))

# ---- Figure 1 --------------------------------------------------------------------
pos <- list(
  A = data.frame(x = c(0, 1, .5, .5, 1, .67, .33, 0),
                 y = c(.75, .75, 1, .5, .2, 0, 0, .2)),
  B = data.frame(x = c(.15, 1, .3, 0, .45, .62, .8, 1),
                 y = c(.45, 1, .9, 1, .45, .9, .45, 0)))
# edge mechanism: reciprocal, on a directed circuit of length >= 3, or one-way
edge_mechanism <- function(A) {
  E <- which(A > 0, arr.ind = TRUE)                    # row = to, col = from
  onc <- Reduce(`+`, lapply(2:(nrow(A) - 1), function(h) simple_path_matrix(A, h)))
  sapply(seq_len(nrow(E)), function(k) {
    to <- E[k, 1]; fr <- E[k, 2]
    if (A[fr, to] > 0) "Reciprocal" else if (onc[fr, to] > 0) "Circuits" else "One-way" })
}
net_plot <- function(A, xy, title, subtitle) {
  E <- which(A > 0, arr.ind = TRUE)
  e <- data.frame(x = xy$x[E[, 2]], y = xy$y[E[, 2]], xend = xy$x[E[, 1]],
                  yend = xy$y[E[, 1]], mech = edge_mechanism(A))
  shrink <- function(e, r = 0.055) { dx <- e$xend - e$x; dy <- e$yend - e$y
    L <- sqrt(dx^2 + dy^2); e$x <- e$x + r * dx / L; e$y <- e$y + r * dy / L
    e$xend <- e$xend - r * dx / L; e$yend <- e$yend - r * dy / L; e }
  e <- shrink(e)
  v <- cbind(xy, id = seq_len(nrow(xy)), heavy = w > 0.1)
  arr <- arrow(length = unit(1.8, "mm"), type = "closed")
  ggplot() +
    geom_segment(data = e[e$mech != "Reciprocal", ],
                 aes(x, y, xend = xend, yend = yend, colour = mech),
                 arrow = arr, linewidth = 0.6) +
    geom_curve(data = e[e$mech == "Reciprocal", ],
               aes(x, y, xend = xend, yend = yend, colour = mech),
               curvature = 0.3, arrow = arr, linewidth = 0.6) +
    geom_point(data = v, aes(x, y, size = heavy, fill = heavy), shape = 21,
               colour = "white", stroke = 0.8) +
    geom_text(data = v[v$heavy, ], aes(x, y, label = id), colour = "white",
              size = 2.6, fontface = "bold", show.legend = FALSE) +
    geom_text(data = v[!v$heavy, ], aes(x, y, label = id), colour = INK,
              size = 2.6, fontface = "bold", show.legend = FALSE) +
    scale_colour_manual(values = COL_MECH, drop = FALSE,
                        limits = c("One-way", "Reciprocal", "Circuits")) +
    scale_fill_manual(values = c(`TRUE` = "#104281", `FALSE` = "#cde2fb"),
                      guide = "none") +
    scale_size_manual(values = c(`TRUE` = 9, `FALSE` = 6), guide = "none") +
    coord_equal(clip = "off", xlim = c(-.05, 1.05), ylim = c(-.05, 1.05)) +
    labs(title = title, subtitle = subtitle) + theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold", colour = INK),
          plot.subtitle = element_text(colour = INK2, size = 7.5),
          legend.position = "none")
}
dfmt <- function(x) formatC(x, format = "f", digits = 2)
sub_net <- function(nm, d12) sprintf("Heavy nodes 1, 2 (w = 0.3) at distance %d", d12)
pa <- net_plot(nets$A, pos$A, "(a) Network A", sub_net("A", 2))
pb <- net_plot(nets$B, pos$B, "(b) Network B", sub_net("B", 5))
# Shared legend of the edge mechanisms, drawn as its own strip and centred
# under panels (a)-(b).  Built by hand so that every mechanism is shown
# whatever the ggplot2 version (from 3.5 legend keys are drawn only for
# values present in a panel, and network B has no circuits).
LEG_X <- c("One-way" = 0, "Reciprocal" = 1.55, "Circuits" = 3.1)
arr <- arrow(length = unit(1.8, "mm"), type = "closed")
leg <- ggplot() +
  geom_segment(aes(x = LEG_X[c(1, 3)], xend = LEG_X[c(1, 3)] + .32, y = 0, yend = 0),
               colour = COL_MECH[c("One-way", "Circuits")], arrow = arr, linewidth = .6) +
  geom_curve(aes(x = LEG_X[2], xend = LEG_X[2] + .32, y = -.08, yend = -.08),
             colour = COL_MECH[["Reciprocal"]], curvature = -.5, arrow = arr, linewidth = .6) +
  geom_curve(aes(x = LEG_X[2] + .32, xend = LEG_X[2], y = .08, yend = .08),
             colour = COL_MECH[["Reciprocal"]], curvature = -.5, arrow = arr, linewidth = .6) +
  annotate("text", x = LEG_X + .42, y = 0, label = names(LEG_X), hjust = 0,
           size = 2.7, colour = INK2) +
  annotate("text", x = -.15, y = 0, label = "Links:", hjust = 1, size = 2.7,
           colour = INK, fontface = "bold") +
  coord_cartesian(xlim = c(-1.05, 4.35), ylim = c(-.4, .4), clip = "off") +
  theme_void()
prof <- rbind(cbind(net = "Network A", r$A$tp$profile),
              cbind(net = "Network B", r$B$tp$profile))
bars <- prof; bars$hx <- bars$h + ifelse(bars$net == "Network A", -0.17, 0.17)
ann <- data.frame(net = c("Network A", "Network B"), x = c(6.5, 6.5), y = c(.445, .415),
  lab = c(sprintf("A: economic mean distance %s < average %s", dfmt(r$A$tp$dw), dfmt(r$A$tp$dbar)),
          sprintf("B: economic mean distance %s > average %s", dfmt(r$B$tp$dw), dfmt(r$B$tp$dbar))))
pc <- ggplot(prof, aes(h, psi, colour = net)) +
  geom_col(data = bars, aes(hx, mbar, fill = net), width = .3, alpha = .22,
           colour = NA, inherit.aes = FALSE, show.legend = FALSE) +
  geom_line(linewidth = .7) + geom_point(size = 2.2, aes(shape = net)) +
  geom_text(data = ann, aes(x, y, label = lab, colour = net), hjust = 1, size = 2.3,
            inherit.aes = FALSE, show.legend = FALSE) +
  scale_colour_manual(values = COL_NET) + scale_fill_manual(values = COL_NET) +
  scale_x_continuous(breaks = 1:6) +
  labs(x = "Geodesic distance h", y = NULL,
       title = expression(bold("(c) Topological profile ")*psi[G]^{(h)}*(w)),
       subtitle = "Lines: weighted interaction at distance h. Bars: share of pairs") +
  coord_cartesian(ylim = c(0, .45)) +
  guides(fill = "none") + theme_paper()
pk <- data.frame(k = rep(1:6, 2), psiO = c(r$A$psiO, r$B$psiO),
                 net = rep(c("Network A", "Network B"), each = 6))
pd <- ggplot(pk, aes(k, psiO, colour = net, shape = net)) +
  geom_line(linewidth = .7) + geom_point(size = 2.2) +
  annotate("text", x = 4.5, y = r$B$psiO[3] - .035, hjust = .5, size = 2.3,
           colour = COL_NET[2], label = "B: trajectories of 1 and 2 converge on 5-6") +
  annotate("text", x = 6, y = .255, hjust = 1, size = 2.3,
           colour = COL_NET[1], label = "A: heavy nodes propagate along different routes") +
  scale_colour_manual(values = COL_NET) + scale_x_continuous(breaks = 1:6) +
  labs(x = "Propagation horizon k", y = NULL,
       title = expression(bold("(d) Convergence index ")*psi[O]^{(k)}*(w)),
       subtitle = "Overlap of the k-step destinations of heavy origins") +
  theme_paper()
# panels (c)-(d) share one legend (networks A and B), collected and centred
bottom <- (pc | pd) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.justification = "center")
fig1 <- wrap_elements(full = (pa | pb)) / leg / wrap_elements(full = bottom) +
  plot_layout(heights = c(1, .07, 1.05))
save_fig(fig1, "fig1_example", width = 7.2, height = 7)

# ---- check against the published values (Table S1 of the paper) ---------------
published <- c(psiD = .247, m1 = .286, G1 = .863, psiG2A = .373, G2A = 1.490,
               psiG2B = .146, G2B = .680, dwA = 2.52, dbarA = 2.64, dwB = 3.25,
               dbarB = 2.79, psiO1A = .011, psiO3A = .034, psiO1B = .048,
               psiO3B = .375, hhi3A = .184, hhi3B = .502, totA = .060,
               totB = .048, d2 = .221, d2M = .185, d2a = .038, d2w = -.002)
computed <- c(r$A$b$psiD, r$A$b$density, r$A$b$Gamma1, r$A$tp$profile$psi[2],
              r$A$tp$profile$Gamma[2], r$B$tp$profile$psi[2],
              r$B$tp$profile$Gamma[2], r$A$tp$dw, r$A$tp$dbar, r$B$tp$dw,
              r$B$tp$dbar, r$A$p1$psiO, r$A$p3$psiO, r$B$p1$psiO, r$B$p3$psiO,
              r$A$p3$hhi_wk, r$B$p3$hhi_wk, r$A$sd$psiTOT, r$B$sd$psiTOT,
              dyn$delta[2], dyn$network[2], dyn$availability[2], dyn$weights[2])
chk <- data.frame(quantity = names(published), published = published,
                  computed = round(computed, 4),
                  match = abs(published - computed) < 0.0051)
write.csv(chk, file.path(DIR_TAB, "check_example_vs_draft.csv"), row.names = FALSE)
log_msg("Example: ", sum(chk$match), "/", nrow(chk), " values match the paper")
if (!all(chk$match)) print(chk[!chk$match, ])
