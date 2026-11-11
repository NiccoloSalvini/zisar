# Explanatory figures of the README and the vignette, after the workshop
# slides: the two kinds of zeros of the zero-inflated negative binomial, and
# how far the spatial multiplier carries a change in one unit.
# Writes man/figures/{zeros,multiplier}.png (README, pkgdown home page) and
# copies them to vignettes/figures/ (vignette).
# Needs zisar, ggplot2 and ragg. Run from the package root:
#   Rscript data-raw/figures.R

library(zisar)
library(ggplot2)

col_struct <- "#D55E00"
col_sampl <- "#56B4E9"
col_count <- "grey75"
navy <- "#1a3461"
theme_zisar <- function(base_size = 13) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
          plot.title = element_text(face = "bold", colour = navy),
          plot.subtitle = element_text(colour = "grey35"),
          axis.title = element_text(colour = "grey25"),
          plot.background = element_rect(fill = "white", colour = NA))
}
save_png <- function(p, name, w, h) {
  f <- file.path("man/figures", name)
  ragg::agg_png(f, width = w, height = h, units = "in", res = 200)
  print(p)
  grDevices::dev.off()
  dir.create("vignettes/figures", showWarnings = FALSE)
  file.copy(f, file.path("vignettes/figures", name), overwrite = TRUE)
}

# 1. The zero bar of a ZINB is the sum of two kinds of zero
pi0 <- 0.25
mu <- 3
kappa <- 1.5
k <- 0:10
p_nb <- dnbinom(k, size = kappa, mu = mu)
d <- rbind(
  data.frame(k = k[-1], part = "count", p = (1 - pi0) * p_nb[-1]),
  data.frame(k = 0, part = "sampling", p = (1 - pi0) * p_nb[1]),
  data.frame(k = 0, part = "structural", p = pi0))
d$part <- factor(d$part, levels = c("structural", "sampling", "count"))
# Labels in the empty area right of the zero bar, linked to its two parts
lab <- data.frame(
  y = c(pi0 / 2 + (1 - pi0) * p_nb[1], (1 - pi0) * p_nb[1] / 2),
  y_text = c(0.31, 0.205),
  text = c("structural zero: the unit is out of the process,\nprobability \u03c0",
           "sampling zero: the negative binomial draws a zero,\nprobability (1 \u2212 \u03c0) NB(0; \u03bc, \u03ba)"),
  col = c(col_struct, "#1f78b4"))
p1 <- ggplot(d, aes(factor(k), p, fill = part)) +
  geom_col(width = 0.75) +
  # Elbow connectors up the gap between the bars of 0 and 1
  annotate("path", x = c(1.40, 1.50, 1.50, 2.18), colour = lab$col[1],
           y = c(lab$y[1], lab$y[1], lab$y_text[1], lab$y_text[1]),
           linewidth = 0.5) +
  annotate("path", x = c(1.40, 1.50, 1.50, 2.18), colour = lab$col[2],
           y = c(lab$y[2], lab$y[2], lab$y_text[2], lab$y_text[2]),
           linewidth = 0.5) +
  annotate("text", x = 2.25, y = lab$y_text, label = lab$text, hjust = 0,
           size = 3.8, colour = lab$col, family = "Helvetica",
           lineheight = 0.95) +
  scale_fill_manual(values = c(structural = col_struct, sampling = col_sampl,
                               count = col_count), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Two kinds of zeros",
       subtitle = sprintf(paste0("Zero-inflated negative binomial with ",
                                 "π = %.2f, μ = %g, κ = %g"),
                          pi0, mu, kappa),
       x = "count y", y = "probability") +
  theme_zisar()
save_png(p1, "zeros.png", 7, 4.3)

# 2. A change in one unit spreads through the spatial multiplier
n_side <- 15
n <- n_side^2
W <- lattice_weights(n_side)
grid <- expand.grid(row = seq_len(n_side), col = seq_len(n_side))
centre <- which(grid$row == 8 & grid$col == 8)
e <- replace(rep(0, n), centre, 1)
d2 <- do.call(rbind, lapply(c(0.2, 0.5, 0.8), function(rho) {
  s <- as.numeric(Matrix::solve(Matrix::Diagonal(n) - rho * W, e))
  data.frame(grid, v = log10(pmax(s, 1e-6)),
             panel = sprintf("rho == %.1f", rho))
}))
p2 <- ggplot(d2, aes(col, row, fill = v)) +
  geom_raster() +
  facet_wrap(~panel, labeller = label_parsed) +
  coord_equal() +
  scale_fill_gradientn(
    colours = c("#ffffff", "#fdd49e", "#fc8d59", "#d7301f", "#7f0000"),
    limits = c(-4, 0.2), breaks = -4:0,
    labels = c("0.0001", "0.001", "0.01", "0.1", "1"),
    oob = scales::squish, name = NULL) +
  labs(title = "How far a change spreads",
       subtitle = paste0("Effect on η of a unit change in the centre ",
                         "cell of a lattice, (I − ρW)⁻¹")) +
  theme_void(base_size = 13, base_family = "Helvetica") +
  theme(plot.title = element_text(face = "bold", colour = navy),
        plot.subtitle = element_text(colour = "grey35",
                                     margin = margin(b = 8)),
        strip.text = element_text(face = "bold", colour = navy, size = 13),
        legend.key.height = unit(0.9, "cm"),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(8, 8, 8, 8))
save_png(p2, "multiplier.png", 8, 3.4)
