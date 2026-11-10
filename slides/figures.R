# =============================================================================
#  Figures for the ZISAR workshop slides (GSSI, 30 October 2026)
#  ggplot2 + deldir (Voronoi polygons) + patchwork. Writes PDFs into figures/.
#  Monte Carlo numbers are copied from tex/montecarlo.tex (R = 100 reps).
#  Section 9 fits the real ZISAR estimator (notebooks/zisar_tmb.cpp, TMB) on the
#  synthetic map; the fit is cached in figures/est_tmb_cache.rds.
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(ggplot2); library(deldir); library(patchwork)
})
out_dir <- "figures"; dir.create(out_dir, showWarnings = FALSE)

# --- Palette and theme (Okabe-Ito, colour-blind safe) -----------------------
col_zinb  <- "#E69F00"; col_sarnb <- "#56B4E9"; col_zisar <- "#D55E00"
col_grey  <- "#666666"; navy <- "#1a3461"
model_cols <- c(ZINB = col_zinb, `SAR-NB` = col_sarnb, ZISAR = col_zisar)
model_shapes <- c(ZINB = 17, `SAR-NB` = 15, ZISAR = 16)

theme_deck <- function(base_size = 14) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey88", linewidth = 0.3),
          strip.text = element_text(face = "bold", size = base_size, colour = navy),
          plot.title = element_text(face = "bold", colour = navy, size = base_size + 1),
          plot.subtitle = element_text(colour = "grey35"),
          legend.position = "bottom", legend.title = element_blank(),
          axis.title = element_text(colour = "grey25"),
          plot.background = element_rect(fill = "white", colour = NA))
}
theme_map <- function(base_size = 13) {
  theme_void(base_family = "Helvetica", base_size = base_size) +
    theme(plot.title = element_text(face = "bold", colour = navy, size = base_size + 1, hjust = 0),
          plot.subtitle = element_text(colour = "grey35", size = base_size - 2),
          strip.text = element_text(face = "bold", size = base_size, colour = navy),
          legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = base_size - 2), legend.key.size = unit(0.45, "cm"),
          plot.background = element_rect(fill = "white", colour = NA))
}
save <- function(p, file, w, h) ggsave(file.path(out_dir, file), p, width = w, height = h, device = "pdf")

# --- Helpers -----------------------------------------------------------------
queen_W <- function(n_side) {
  n <- n_side^2; coords <- expand.grid(row = 1:n_side, col = 1:n_side)
  W <- matrix(0, n, n)
  for (i in 1:n) {
    dr <- abs(coords$row - coords$row[i]); dc <- abs(coords$col - coords$col[i])
    W[i, dr <= 1 & dc <= 1] <- 1
  }
  diag(W) <- 0; W / rowSums(W)
}
sim_zisar <- function(n_side, beta, rho, gamma, kappa, sigma, seed) {
  set.seed(seed); n <- n_side^2; W <- queen_W(n_side)
  X <- cbind(1, matrix(rnorm(n * (length(beta) - 1)), n))
  Z <- cbind(1, matrix(rnorm(n * (length(gamma) - 1)), n))
  u <- rnorm(n, 0, sigma)
  eta <- as.numeric(solve(diag(n) - rho * W, X %*% beta + u))
  s <- rbinom(n, 1, plogis(as.numeric(Z %*% gamma)))
  y <- ifelse(s == 1, 0L, rnbinom(n, size = kappa, mu = exp(eta)))
  list(y = y, eta = eta, W = W, X = X, u = u, n_side = n_side)
}
lattice_df <- function(v, n_side) data.frame(x = rep(1:n_side, n_side), y = rep(1:n_side, each = n_side), v = v)

# =============================================================================
# 1. Excess zeros: observed frequencies vs Poisson and NB with the same moments
# =============================================================================
# (fig_zeros is drawn in section 8, on the same synthetic map as the map slides)

# How a ZINB works: the zero bar is the sum of two kinds of zero
pi_z <- 0.25; mu_z <- 3; k_z <- 1.5; Kz <- 10
pz <- dnbinom(0:Kz, size = k_z, mu = mu_z)
dz <- rbind(data.frame(k = 0:Kz, part = "count from the NB", p = (1 - pi_z) * pz * (0:Kz > 0)),
            data.frame(k = 0, part = "sampling zero: NB gives 0", p = (1 - pi_z) * pz[1]),
            data.frame(k = 0, part = "structural zero: switch off", p = pi_z))
dz$part <- factor(dz$part, levels = c("structural zero: switch off", "sampling zero: NB gives 0", "count from the NB"))
pZ <- ggplot(dz, aes(factor(k), p, fill = part)) + geom_col(width = 0.75) +
  scale_fill_manual(values = c("#D55E00", "#56B4E9", "grey75")) +
  labs(x = "count y", y = "probability") + theme_deck(15) + guides(fill = guide_legend(nrow = 3))
save(pZ, "fig_zinb.pdf", 6.5, 5)

# =============================================================================
# 2. Latent field on a lattice: rho = 0 vs rho = 0.6, same X and u
# =============================================================================
n_side <- 20; set.seed(7); W <- queen_W(n_side); n <- n_side^2
Xb <- 0.5 + 0.8 * rnorm(n) - 0.5 * rnorm(n); u <- rnorm(n, 0, 0.5)
df2 <- rbind(cbind(lattice_df(Xb + u, n_side), panel = "rho = 0 (aspatial)"),
             cbind(lattice_df(as.numeric(solve(diag(n) - 0.6 * W, Xb + u)), n_side), panel = "rho = 0.6, same X*beta + u"))
p2 <- ggplot(df2, aes(x, y, fill = v)) + geom_raster() + facet_wrap(~panel) + coord_equal() +
  scale_fill_gradient2(low = "#2c7bb6", mid = "#f7f7f7", high = "#d7191c", midpoint = median(df2$v)) +
  theme_map() + theme(legend.position = "none")
save(p2, "fig_lattice.pdf", 8.5, 4.4)

# =============================================================================
# 3. Spatial multiplier: response to a unit shock in the centre of a lattice
# =============================================================================
n_side <- 15; n <- n_side^2; W <- queen_W(n_side)
centre <- (n_side %/% 2) * n_side + n_side %/% 2 + 1
e <- rep(0, n); e[centre] <- 1
df3 <- do.call(rbind, lapply(c(0.2, 0.5, 0.8), function(rho) {
  cbind(lattice_df(log10(pmax(as.numeric(solve(diag(n) - rho * W, e)), 1e-6)), n_side), panel = sprintf("rho == %.1f", rho))
}))
p3 <- ggplot(df3, aes(x, y, fill = v)) + geom_raster() + facet_wrap(~panel, labeller = label_parsed) + coord_equal() +
  scale_fill_gradientn(colours = c("#ffffff", "#fdd49e", "#fc8d59", "#d7301f", "#7f0000"), limits = c(-4, 0.2),
                       breaks = c(-4, -3, -2, -1, 0), labels = c("1e-4", "1e-3", "0.01", "0.1", "1"), name = NULL, oob = scales::squish) +
  theme_map() + theme(legend.position = "right", legend.key.height = unit(0.9, "cm"))
save(p3, "fig_multiplier.pdf", 10, 3.9)

# =============================================================================
# 4. Average direct vs total effect as a function of rho (lattice W)
# =============================================================================
rhos <- seq(0, 0.8, by = 0.05)
df4 <- do.call(rbind, lapply(rhos, function(r) {
  S <- solve(diag(n) - r * W)
  data.frame(rho = r, effect = c("total effect (own unit + spillover to all others)", "direct effect (own unit)"),
             value = c(sum(S) / n, sum(diag(S)) / n))
}))
p4 <- ggplot(df4, aes(rho, value, colour = effect)) + geom_hline(yintercept = 1, linetype = 3, colour = "grey60") +
  geom_line(linewidth = 1.3) + scale_colour_manual(values = c(col_grey, col_zisar)) +
  labs(x = expression(rho), y = "average effect multiplier",
       title = "Direct vs total effect of a unit change in a covariate") +
  theme_deck() + guides(colour = guide_legend(nrow = 2))
save(p4, "fig_effects.pdf", 7.5, 4.8)

# =============================================================================
# 5-7. Monte Carlo results (from tex/montecarlo.tex)
# =============================================================================
ns <- c(100, 225, 400); rl <- c("rho[0] == 0", "rho[0] == 0.3", "rho[0] == 0.6")
mk <- function(model, metric, vals) data.frame(model = model, metric = metric, rho = rep(rl, each = 3), n = rep(ns, 3), value = vals)
mc_rho <- rbind(
  mk("SAR-NB", "bias", c(0.066, -0.006, 0.018, -0.170, -0.180, -0.184, -0.464, -0.474, -0.479)),
  mk("ZISAR",  "bias", c(-0.019, -0.037, -0.015, -0.105, -0.032, 0.002, -0.062, -0.000, -0.007)),
  mk("SAR-NB", "rmse", c(0.241, 0.166, 0.113, 0.264, 0.204, 0.198, 0.474, 0.478, 0.481)),
  mk("ZISAR",  "rmse", c(0.390, 0.342, 0.226, 0.423, 0.199, 0.148, 0.284, 0.122, 0.089)))
mc_rho$rho <- factor(mc_rho$rho, levels = rl)
mc_plot <- function(df, ylab, title, hline = NULL) {
  p <- ggplot(df, aes(n, value, colour = model, shape = model))
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = 3, colour = "grey55")
  p + geom_line(linewidth = 1.2) + geom_point(size = 3) +
    facet_wrap(~rho, labeller = label_parsed) + scale_x_continuous(breaks = ns) +
    scale_colour_manual(values = model_cols) + scale_shape_manual(values = model_shapes) +
    labs(x = "n", y = ylab, title = title) + theme_deck() + theme(panel.spacing = unit(1.2, "lines"))
}
save(mc_plot(subset(mc_rho, metric == "bias"), expression(Bias(hat(rho))), "Bias of the spatial parameter, R = 100", 0),
     "fig_mc_rho.pdf", 10, 4.3)
save(mc_plot(subset(mc_rho, metric == "rmse"), expression(RMSE(hat(rho))), "RMSE of the spatial parameter, R = 100"),
     "fig_mc_rho_rmse.pdf", 10, 4.3)

b_bias <- rbind(mk("ZINB", "b1", c(0.018, -0.009, 0.006, -0.003, -0.000, 0.009, 0.053, 0.059, 0.055)),
                mk("SAR-NB", "b1", c(0.009, -0.017, 0.004, -0.018, -0.011, -0.009, 0.016, 0.002, 0.000)),
                mk("ZISAR", "b1", c(-0.012, -0.014, 0.008, -0.020, -0.008, -0.003, 0.011, -0.002, -0.002)))
b_bias$rho <- factor(b_bias$rho, levels = rl)
b_rmse <- rbind(
  data.frame(model = "ZINB",   par = "hat(beta)[1]", n = ns, value = c(0.167, 0.121, 0.099)),
  data.frame(model = "SAR-NB", par = "hat(beta)[1]", n = ns, value = c(0.151, 0.111, 0.090)),
  data.frame(model = "ZISAR",  par = "hat(beta)[1]", n = ns, value = c(0.154, 0.102, 0.076)),
  data.frame(model = "ZINB",   par = "hat(beta)[2]", n = ns, value = c(0.175, 0.117, 0.083)),
  data.frame(model = "SAR-NB", par = "hat(beta)[2]", n = ns, value = c(0.172, 0.116, 0.086)),
  data.frame(model = "ZISAR",  par = "hat(beta)[2]", n = ns, value = c(0.154, 0.092, 0.066)))
pA <- mc_plot(b_bias, expression(Bias(hat(beta)[1])), expression("Bias of " * hat(beta)[1] * " (true 0.8)"), 0)
pB <- ggplot(b_rmse, aes(n, value, colour = model, shape = model)) + geom_line(linewidth = 1.2) + geom_point(size = 3) +
  facet_wrap(~par, labeller = label_parsed) + scale_x_continuous(breaks = ns) +
  scale_colour_manual(values = model_cols) + scale_shape_manual(values = model_shapes) +
  labs(x = "n", y = "RMSE", title = expression("RMSE at " * rho[0] == 0.6)) + theme_deck() + theme(panel.spacing = unit(1.2, "lines"))
save((pA | pB) + plot_layout(widths = c(3, 2), guides = "collect") & theme(legend.position = "bottom"), "fig_mc_beta.pdf", 12.5, 4.4)

kap <- rbind(mk("ZINB", "k", c(-0.294, -0.648, -0.763, -0.384, -0.800, -0.873, -0.874, -1.062, -1.074)),
             mk("SAR-NB", "k", c(-1.354, -1.421, -1.443, -1.407, -1.433, -1.458, -1.489, -1.509, -1.511)),
             mk("ZISAR", "k", c(-0.033, -0.311, -0.406, -0.025, -0.421, -0.387, 0.373, 0.322, -0.202)))
kap_mean <- aggregate(value ~ model + n, kap, mean)
pK <- ggplot(kap_mean, aes(n, value, colour = model, shape = model)) + geom_hline(yintercept = 0, linetype = 3, colour = "grey55") +
  geom_line(linewidth = 1.2) + geom_point(size = 3) + scale_x_continuous(breaks = ns) +
  scale_colour_manual(values = model_cols) + scale_shape_manual(values = model_shapes) +
  labs(x = "n", y = expression(Bias(hat(kappa))), title = expression("NB shape, " * kappa[0] == 2 * " (average over " * rho[0] * ")")) + theme_deck()
tm <- data.frame(n = ns, s = c(0.3, 0.7, 1.4))
pT <- ggplot(tm, aes(n, s)) + geom_line(colour = col_zisar, linewidth = 1.2) + geom_point(colour = col_zisar, size = 3) +
  geom_text(aes(label = paste0(s, " s")), vjust = -1, colour = col_zisar, size = 4.2) +
  scale_x_continuous(breaks = ns, limits = c(70, 430)) + scale_y_continuous(limits = c(0, 1.8)) +
  labs(x = "n", y = "seconds per fit", title = "ZISAR computation time (TMB / Laplace)") + theme_deck()
save((pK | pT) + plot_layout(widths = c(1.15, 1)), "fig_mc_kappa_time.pdf", 10, 4.4)


# --- R = 1000 Monte Carlo (notebooks/mc_run_1000.r): replaces the R = 100 rho and beta figures
mc_file <- file.path("..", "notebooks", "mc_results_1000.rds")
if (file.exists(mc_file)) {
  mc <- readRDS(mc_file)
  mc$panel <- factor(paste0("n = ", mc$n), levels = paste0("n = ", c(100, 225, 400)))
  ok_z <- mc$zisar_conv %in% TRUE & is.finite(mc$zisar_rho) & mc$zisar_kappa < 100
  ok_s <- mc$sarnb_conv %in% TRUE & is.finite(mc$sarnb_rho)
  ok_i <- mc$zinb_conv %in% TRUE & is.finite(mc$zinb_b1)
  long <- rbind(data.frame(model = "ZISAR",  panel = mc$panel[ok_z], rho_true = mc$rho_true[ok_z], rho = mc$zisar_rho[ok_z], b1 = mc$zisar_b1[ok_z]),
                data.frame(model = "SAR-NB", panel = mc$panel[ok_s], rho_true = mc$rho_true[ok_s], rho = mc$sarnb_rho[ok_s], b1 = mc$sarnb_b1[ok_s]),
                data.frame(model = "ZINB",   panel = mc$panel[ok_i], rho_true = mc$rho_true[ok_i], rho = NA, b1 = mc$zinb_b1[ok_i]))
  sm <- do.call(rbind, lapply(split(long, list(long$model, long$panel, long$rho_true), drop = TRUE), function(d) data.frame(
    model = d$model[1], panel = d$panel[1], rho_true = d$rho_true[1], reps = nrow(d),
    rho_mean = mean(d$rho), rho_lo = quantile(d$rho, 0.1, na.rm = TRUE), rho_hi = quantile(d$rho, 0.9, na.rm = TRUE),
    rho_rmse = sqrt(mean((d$rho - d$rho_true)^2)), b1_mean = mean(d$b1), b1_rmse = sqrt(mean((d$b1 - 0.8)^2)))))
  write.csv(sm, file.path(out_dir, "mc_1000_summary.csv"), row.names = FALSE)
  smr <- subset(sm, model != "ZINB")
  zi <- subset(sm, model == "ZINB"); zi$rho_mean <- 0; zi$rho_lo <- 0; zi$rho_hi <- 0; zi$rho_rmse <- zi$rho_true
  smr <- rbind(smr, zi)
  pRho <- ggplot(smr, aes(rho_true, rho_mean, colour = model, fill = model, shape = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45") +
    geom_ribbon(aes(ymin = rho_lo, ymax = rho_hi), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 1.2) + geom_point(size = 2.8) + facet_wrap(~panel) +
    scale_colour_manual(values = model_cols) + scale_fill_manual(values = model_cols) + scale_shape_manual(values = model_shapes) +
    scale_x_continuous(breaks = c(0, 0.15, 0.3, 0.45, 0.6)) +
    labs(x = expression("true " * rho), y = expression("estimated " * rho)) +
    theme_deck() + theme(panel.spacing = unit(1.2, "lines"))
  save(pRho, "fig_mc_rho.pdf", 10, 4.4)
  pRm <- ggplot(smr, aes(rho_true, rho_rmse, colour = model, shape = model)) +
    geom_line(linewidth = 1.2) + geom_point(size = 2.8) + facet_wrap(~panel) +
    scale_colour_manual(values = model_cols) + scale_shape_manual(values = model_shapes) +
    scale_x_continuous(breaks = c(0, 0.15, 0.3, 0.45, 0.6)) +
    labs(x = expression("true " * rho), y = expression(RMSE(hat(rho)))) + theme_deck() + theme(panel.spacing = unit(1.2, "lines"))
  save(pRm, "fig_mc_rho_rmse.pdf", 10, 4.3)
  pB1 <- ggplot(sm, aes(rho_true, b1_mean, colour = model, shape = model)) +
    geom_hline(yintercept = 0.8, linetype = 2, colour = "grey45") +
    geom_line(linewidth = 1.2) + geom_point(size = 2.8) + facet_wrap(~panel) +
    scale_colour_manual(values = model_cols) + scale_shape_manual(values = model_shapes) +
    scale_x_continuous(breaks = c(0, 0.15, 0.3, 0.45, 0.6)) +
    labs(x = expression("true " * rho), y = expression("average " * hat(beta)[1] * "  (true 0.8)")) +
    theme_deck() + theme(panel.spacing = unit(1.2, "lines"))
  save(pB1, "fig_mc_beta.pdf", 10, 4.4)
}

# =============================================================================
# 8. The synthetic region: 300 municipalities as Voronoi polygons
# =============================================================================
set.seed(31)
n_m <- 300
centres <- cbind(runif(5, 0.15, 0.85), runif(5, 0.15, 0.85))
urban <- do.call(rbind, lapply(1:5, function(k) cbind(rnorm(25, centres[k,1], 0.05), rnorm(25, centres[k,2], 0.05))))
rural <- cbind(runif(n_m - nrow(urban)), runif(n_m - nrow(urban)))
pts <- rbind(urban, rural); pts <- pmin(pmax(pts, 0.005), 0.995)
hosp <- centres[1:4, ]
d_h  <- apply(pts, 1, function(p) min(sqrt((hosp[,1]-p[1])^2 + (hosp[,2]-p[2])^2)))
d_c  <- apply(pts, 1, function(p) min(sqrt((centres[,1]-p[1])^2 + (centres[,2]-p[2])^2)))
logpop <- as.numeric(scale(-d_c + rnorm(n_m, 0, 0.25)))
dist_h <- as.numeric(scale(d_h))
D <- as.matrix(dist(pts)); diag(D) <- Inf
Wm <- t(apply(D, 1, function(r) { w <- numeric(n_m); w[order(r)[1:6]] <- 1/6; w }))
rho_m <- 0.5
u_m <- rnorm(n_m, 0, 0.5)
eta_m <- as.numeric(solve(diag(n_m) - rho_m * Wm, 0.8 + 0.7 * logpop - 0.2 * dist_h + u_m))
pi_m <- plogis(-1.6 + 1.4 * dist_h)
s_m <- rbinom(n_m, 1, pi_m)
y_m <- ifelse(s_m == 1, 0L, rnbinom(n_m, size = 2, mu = exp(eta_m)))

vor <- deldir(pts[,1], pts[,2], rw = c(0, 1, 0, 1))
tiles <- tile.list(vor)
poly <- do.call(rbind, lapply(seq_along(tiles), function(i) data.frame(id = i, x = tiles[[i]]$x, y = tiles[[i]]$y)))
hosp_df <- data.frame(x = hosp[,1], y = hosp[,2])

map_base <- function(fillvar, title = NULL) {
  pd <- poly; pd$fill <- fillvar[pd$id]
  ggplot(pd, aes(x, y, group = id, fill = fill)) +
    geom_polygon(colour = "grey35", linewidth = 0.25) + coord_equal(expand = FALSE) +
    labs(title = title) + theme_map()
}
hosp_layer <- function() geom_point(data = hosp_df, aes(x, y), inherit.aes = FALSE, shape = 3, size = 4, stroke = 1.6, colour = "black")

count_breaks <- c(-0.5, 0.5, 2.5, 5.5, 10.5, 20.5, Inf)
count_labs <- c("0", "1-2", "3-5", "6-10", "11-20", "> 20")
seq_cols <- c("#ffffff", "#fee8c8", "#fdbb84", "#fc8d59", "#e34a33", "#b30000")
y_cls <- cut(y_m, count_breaks, labels = count_labs)

m1 <- map_base(y_cls) + hosp_layer() +
  scale_fill_manual(values = seq_cols, drop = FALSE) + guides(fill = guide_legend(nrow = 1))
save(m1, "map_1_counts.pdf", 6.4, 6.4)

m2 <- map_base(factor(ifelse(y_m == 0, "y = 0", "y > 0"), levels = c("y > 0", "y = 0"))) + hosp_layer() +
  scale_fill_manual(values = c("#f2f2f2", "#4d4d4d"))
save(m2, "map_2_zeros.pdf", 6.4, 6.3)

kind <- factor(ifelse(y_m > 0, "y > 0", ifelse(s_m == 1, "structural zero: no demand reaches an ED", "sampling zero: active, low intensity")),
               levels = c("y > 0", "sampling zero: active, low intensity", "structural zero: no demand reaches an ED"))
m3 <- map_base(kind) + hosp_layer() +
  scale_fill_manual(values = c("#f2f2f2", "#56B4E9", "#D55E00")) + guides(fill = guide_legend(nrow = 3))
d_any <- d_h
iA <- which(y_m == 0 & s_m == 0)[which.min(d_any[y_m == 0 & s_m == 0])]   # sampling zero closest to an ED
iB <- which(s_m == 1)[which.min(d_any[s_m == 1])]                         # structural zero closest to an ED
iC <- which(s_m == 1)[which.max(d_any[s_m == 1])]                         # structural zero farthest from any ED
mk_df <- data.frame(x = pts[c(iA, iB, iC), 1], y = pts[c(iA, iB, iC), 2], lab = c("A", "B", "C"))
m3 <- m3 + geom_label(data = mk_df, aes(x, y, label = lab), inherit.aes = FALSE, fill = "white",
                      fontface = "bold", size = 5.5, label.padding = unit(0.18, "lines"))
save(m3, "map_3_kinds.pdf", 6.4, 6.7)


# Excess zeros on the map data: observed frequencies vs Poisson and NB with the same moments
K <- 15
obs <- tabulate(pmin(y_m, K) + 1, nbins = K + 1) / n_m
m <- mean(y_m); v <- var(y_m); kappa_mm <- m^2 / (v - m)
labs_k <- factor(c(0:(K - 1), paste0(K, "+")), levels = c(0:(K - 1), paste0(K, "+")))
df1 <- data.frame(lab = labs_k, observed = obs)
tail_p <- function(p) c(p[1:K], 1 - sum(p[1:K]))
df1l <- rbind(data.frame(lab = labs_k, model = "Poisson", p = tail_p(dpois(0:(K - 1), m))),
              data.frame(lab = labs_k, model = "NB", p = tail_p(dnbinom(0:(K - 1), size = kappa_mm, mu = m))))
p1 <- ggplot(df1, aes(lab, observed)) +
  geom_col(fill = "grey82", width = 0.8) +
  geom_line(data = df1l, aes(lab, p, colour = model, group = model), linewidth = 1.2) +
  geom_point(data = df1l, aes(lab, p, colour = model, shape = model), size = 2.6) +
  scale_colour_manual(values = c(NB = col_zinb, Poisson = col_sarnb),
                      labels = c(NB = bquote("NB (" * kappa == .(sprintf("%.2f", kappa_mm)) * "; true 2)"), Poisson = "Poisson")) +
  scale_shape_manual(values = c(NB = 17, Poisson = 16),
                     labels = c(NB = bquote("NB (" * kappa == .(sprintf("%.2f", kappa_mm)) * "; true 2)"), Poisson = "Poisson")) +
  labs(x = "count y", y = "share of municipalities") +
  theme_deck()
save(p1, "fig_zeros.pdf", 7.5, 4.8)
cat(sprintf("fig_zeros: zeros %.3f, Poisson P0 %.4f, NB kappa_mm %.3f, NB P0 %.3f\n", mean(y_m == 0), dpois(0, m), kappa_mm, dnbinom(0, size = kappa_mm, mu = m)))

# A clinic in one municipality: change in expected ED visits across the map (true parameters)
j_c <- which.min((pts[,1] - 0.62)^2 + (pts[,2] - 0.42)^2)
S_m <- solve(diag(n_m) - rho_m * Wm)
d_eta <- -0.3 * S_m[, j_c]
Ey0 <- (1 - pi_m) * exp(eta_m); Ey1 <- (1 - pi_m) * exp(eta_m + d_eta)
pct <- 100 * (exp(d_eta) - 1)
eff_cls <- cut(pct, c(-Inf, -20, -5, -1, -0.1, Inf),
               labels = c("more than -20%", "-20% to -5%", "-5% to -1%", "-1% to -0.1%", "no change"))
cl <- data.frame(x = pts[j_c, 1], y = pts[j_c, 2])
mE <- map_base(eff_cls) + hosp_layer() +
  geom_point(data = cl, aes(x, y), inherit.aes = FALSE, shape = 23, fill = "white", size = 4.5, stroke = 1.2) +
  scale_fill_manual(values = c("#08306b", "#2171b5", "#6baed6", "#c6dbef", "#f7f7f7"), drop = FALSE) +
  guides(fill = guide_legend(nrow = 2))
save(mE, "fig_effects_map.pdf", 6.4, 6.6)
dE <- Ey1 - Ey0
cat(sprintf("clinic at %d: direct %.2f visits (%.1f%%), indirect %.2f visits over %d municipalities, total %.2f; share indirect %.2f; first-ring mean pct %.1f\n",
            j_c, dE[j_c], pct[j_c], sum(dE[-j_c]), sum(abs(pct[-j_c]) > 0.1), sum(dE), sum(dE[-j_c]) / sum(dE),
            mean(pct[Wm[j_c, ] > 0])))

eta_lim <- as.numeric(quantile(eta_m, c(0.02, 0.98)))
eta_scale <- function() scale_fill_gradient2(low = "#2c7bb6", mid = "#f7f7f7", high = "#d7191c",
                                             midpoint = median(eta_m), limits = eta_lim, oob = scales::squish, name = NULL)
m4 <- map_base(eta_m) + eta_scale() +
  theme(legend.position = "right", legend.key.height = unit(0.9, "cm"))
save(m4, "map_4_eta.pdf", 6.8, 6.0)

# =============================================================================
# 9. Estimation walk-through: the real estimator (notebooks/zisar_tmb.cpp,
#    TMB / Laplace), run step by step on the synthetic map with the same
#    starting values, bounds and optimiser as zisar_fit(). Cached in
#    figures/est_tmb_cache.rds (delete it to re-run).
# =============================================================================
cache <- file.path(out_dir, "est_tmb_cache.rds")
X_m <- cbind(1, logpop, dist_h); Z_m <- cbind(1, dist_h)
if (file.exists(cache)) {
  E <- readRDS(cache); list2env(E, environment())
} else {
  suppressPackageStartupMessages(library(TMB))
  tmb_dir <- normalizePath(file.path("..", "notebooks"))
  dll <- file.path(tmb_dir, dynlib("zisar_tmb"))
  if (!file.exists(dll)) { owd <- setwd(tmb_dir); compile("zisar_tmb.cpp"); setwd(owd) }
  dyn.load(dll)

  eig <- eigen(Wm, only.values = TRUE)$values
  tmb_data <- list(y = as.numeric(y_m), X = X_m, Z = Z_m,
                   W = as(as(Matrix(Wm, sparse = TRUE), "CsparseMatrix"), "TsparseMatrix"),
                   eigenvalues_W_re = Re(eig), eigenvalues_W_im = Im(eig))
  # starting values exactly as in zisar_fit()
  glm0 <- glm(y_m ~ logpop + dist_h, family = poisson())
  beta0 <- unname(coef(glm0))
  gamma0 <- c(qlogis(max(0.01, min(0.99, mean(y_m == 0)))), 0)
  eta_start <- as.numeric(X_m %*% beta0)
  pars0 <- list(beta = beta0, gamma = gamma0, rho_raw = 0, log_kappa = 0, log_sigma = log(0.5), eta = eta_start)
  lower <- c(rep(-Inf, 5), -4, log(0.01), log(0.001)); upper <- c(rep(Inf, 5), 4, log(500), log(10))

  obj <- MakeADFun(tmb_data, pars0, random = "eta", DLL = "zisar_tmb", silent = TRUE)

  # Step 1-2: one inner loop at theta(0), traced
  newtonOption(obj, trace = 1)
  inner_log <- capture.output(obj$fn(obj$par), type = "output")
  newtonOption(obj, trace = 0)
  mgc <- as.numeric(sub(".*mgc:\\s*([^ ]+).*", "\\1", grep("mgc:", inner_log, value = TRUE)))
  eta_theta0 <- obj$env$last.par[obj$env$random]

  # Step 3: outer loop, traced
  tr_rows <- list()
  fn_tr <- function(p) {
    v <- obj$fn(p)
    tr_rows[[length(tr_rows) + 1]] <<- unname(c(v, tanh(p[6]), exp(p[7]), exp(p[8]), p[2]))
    v
  }
  t0 <- Sys.time()
  opt <- nlminb(obj$par, fn_tr, obj$gr, lower = lower, upper = upper,
                control = list(eval.max = 1000, iter.max = 500, rel.tol = 1e-10))
  sdr <- sdreport(obj)
  t_fit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  tr <- do.call(rbind, tr_rows); colnames(tr) <- c("nll", "rho", "kappa", "sigma", "b1")
  par_best <- obj$env$last.par.best
  eta_hat <- par_best[obj$env$random]
  H <- obj$env$spHess(par_best, random = TRUE)
  H_nnz <- Matrix::nnzero(H) / prod(dim(H))
  Hs_idx <- Matrix::summary(as(H, "TsparseMatrix"))[, 1:2]

  sf <- summary(sdr, "fixed"); sr <- summary(sdr, "report")
  est <- data.frame(parameter = c("beta0", "beta_logpop", "beta_dist", "rho", "gamma0", "gamma_dist", "kappa", "sigma"),
                    truth = c(0.8, 0.7, -0.2, rho_m, -1.6, 1.4, 2, 0.5),
                    estimate = c(sf[1:3, 1], sr["rho", 1], sf[4:5, 1], sr["kappa", 1], sr["sigma", 1]),
                    se = c(sf[1:3, 2], sr["rho", 2], sf[4:5, 2], sr["kappa", 2], sr["sigma", 2]))
  ll_LA <- -opt$objective

  # LR test of rho = 0: same template with rho fixed at 0
  pars_r <- pars0
  obj0 <- MakeADFun(tmb_data, pars_r, random = "eta", DLL = "zisar_tmb", silent = TRUE, map = list(rho_raw = factor(NA)))
  opt0 <- nlminb(obj0$par, obj0$fn, obj0$gr, lower = lower[-6], upper = upper[-6],
                 control = list(eval.max = 1000, iter.max = 500, rel.tol = 1e-10))
  ll_0 <- -opt0$objective
  lr_rho <- 2 * (ll_LA - ll_0)

  # aspatial benchmarks on the same data, for the talk
  nb0 <- MASS::glm.nb(y_m ~ logpop + dist_h)
  zinb0 <- tryCatch(pscl::zeroinfl(y_m ~ logpop + dist_h | dist_h, dist = "negbin"), error = function(e) NULL)

  E <- list(beta0 = beta0, gamma0 = gamma0, eta_start = eta_start, mgc = mgc, inner_log = inner_log,
            eta_theta0 = eta_theta0, tr = tr, opt = opt, t_fit = t_fit, eta_hat = eta_hat,
            H_nnz = H_nnz, Hs_idx = Hs_idx, est = est, ll_LA = ll_LA, ll_0 = ll_0, lr_rho = lr_rho,
            nb_coef = coef(nb0), nb_theta = nb0$theta,
            zinb_coef = if (!is.null(zinb0)) coef(zinb0) else NULL, zinb_theta = if (!is.null(zinb0)) zinb0$theta else NA)
  saveRDS(E, cache)
}
list2env(E, environment())

sink(file.path(out_dir, "est_results.txt"))
cat("Synthetic region: n =", n_m, " zeros =", sum(y_m == 0), sprintf("(%.0f%%)", 100 * mean(y_m == 0)),
    " mean y =", round(mean(y_m), 2), " var y =", round(var(y_m), 1), "\n")
cat("Estimator: notebooks/zisar_tmb.cpp via TMB (Laplace, AD), nlminb, starting values as in zisar_fit()\n")
cat("Start (Poisson GLM): beta =", round(beta0, 3), " gamma =", round(gamma0, 3), "\n")
cat("Aspatial NB on the same data: beta =", round(nb_coef, 3), " kappa =", round(nb_theta, 3), "\n")
if (!is.null(zinb_coef)) cat("Aspatial ZINB on the same data: coef =", round(zinb_coef, 3), " kappa =", round(zinb_theta, 3), "\n")
cat("Inner Newton at theta0: iterations =", length(mgc), " max|gradient| =", signif(mgc, 3), "\n")
cat("Outer nlminb: iterations =", opt$iterations, " fn evals =", opt$evaluations[1], " gr evals =", opt$evaluations[2],
    " convergence =", opt$convergence, " (", opt$message, ")\n")
cat("Wall time (fit + sdreport):", round(t_fit, 2), "s\n")
cat("Hessian of f in eta at the optimum: non-zero share =", round(100 * H_nnz, 2), "%\n")
cat("Laplace logLik ZISAR =", round(ll_LA, 2), "  rho = 0 =", round(ll_0, 2), "  LR(rho=0) =", round(lr_rho, 2),
    "  p =", signif(pchisq(lr_rho, 1, lower.tail = FALSE), 2), "\n")
print(est, digits = 3)
cat("corr(eta_hat, eta_true) =", round(cor(eta_hat, eta_m), 3), "\n")
cat("Outer trace (all fn evaluations):\n"); print(round(tr, 3))
sink()

# est_1: inputs and the neighbourhood of one municipality
focus <- which.min((pts[,1] - 0.62)^2 + (pts[,2] - 0.42)^2); nbrs <- which(Wm[focus, ] > 0)
seg <- data.frame(x = pts[focus,1], y = pts[focus,2], xend = pts[nbrs,1], yend = pts[nbrs,2])
e1 <- map_base(y_cls) + hosp_layer() +
  scale_fill_manual(values = seq_cols, drop = FALSE, guide = "none") +
  geom_segment(data = seg, aes(x, y, xend = xend, yend = yend), inherit.aes = FALSE, colour = navy, linewidth = 1) +
  geom_point(data = data.frame(x = pts[nbrs,1], y = pts[nbrs,2]), aes(x, y), inherit.aes = FALSE, shape = 21, fill = col_sarnb, size = 3.6, stroke = 0.8) +
  geom_point(data = data.frame(x = pts[focus,1], y = pts[focus,2]), aes(x, y), inherit.aes = FALSE, shape = 21, fill = navy, size = 5, stroke = 0.8) +
  NULL
save(e1, "est_1_data.pdf", 6.4, 6.3)

two_maps <- function(v1, v2, t1, t2) {
  pd <- rbind(cbind(poly, fill = v1[poly$id], panel = t1), cbind(poly, fill = v2[poly$id], panel = t2))
  pd$panel <- factor(pd$panel, levels = c(t1, t2))
  ggplot(pd, aes(x, y, group = id, fill = fill)) + geom_polygon(colour = "grey35", linewidth = 0.2) +
    facet_wrap(~panel, strip.position = "bottom") + coord_equal(expand = FALSE) + eta_scale() +
    theme_map(16) + theme(legend.position = "right", legend.key.height = unit(0.9, "cm"),
                          strip.text = element_text(face = "plain", colour = "grey35", size = 13))
}
save(map_base(eta_start) + eta_scale() + hosp_layer() +
       theme(legend.position = "right", legend.key.height = unit(0.9, "cm")),
     "est_2_start.pdf", 6.8, 6.0)

save(two_maps(eta_hat, eta_m, "estimated", "true (simulated)"),
     "est_5_result.pdf", 10.5, 5)

# est_4: outer loop trace
it <- seq_len(nrow(tr))
rd <- data.frame(it = it, rho = tr[, "rho"])
pR <- ggplot(rd, aes(it, rho)) + geom_hline(yintercept = rho_m, linetype = 2, colour = "grey45") +
  annotate("text", x = max(it), y = rho_m, label = "true value 0.5", vjust = -0.6, hjust = 1, colour = "grey35", size = 5) +
  geom_line(colour = col_zisar, linewidth = 1.3) + geom_point(colour = col_zisar, size = 2) +
  labs(x = "successive trial values (each one: rebuild the latent map, score the fit)", y = expression(hat(rho))) + theme_deck(17)
save(pR, "est_4_outer.pdf", 8, 4.6)

cat("All figures written to", normalizePath(out_dir), "\n")

# =============================================================================
# 10. Power study: Moran's I vs likelihood tests on rho (notebooks/power_moran.r)
# =============================================================================
pw_file <- file.path("..", "notebooks", "power_moran.rds")
if (file.exists(pw_file)) {
  pw <- readRDS(pw_file)
  rate <- function(p) mean(p < 0.05, na.rm = TRUE)
  pr <- do.call(rbind, lapply(split(pw, list(pw$n, pw$rho)), function(d) data.frame(
    n = d$n[1], rho = d$rho[1],
    test = c("Moran on y", "Moran on ZINB residuals", "ZISAR Wald on rho", "ZISAR LR on rho"),
    power = c(rate(d$moran_y), rate(d$moran_zinb), rate(d$wald_rho), rate(d$lr_rho)))))
  pr$test <- factor(pr$test, levels = c("Moran on ZINB residuals", "Moran on y", "ZISAR Wald on rho", "ZISAR LR on rho"))
  pr$panel <- factor(paste0("n = ", pr$n), levels = c("n = 300", "n = 900"))
  pP <- ggplot(pr, aes(rho, power, colour = test, shape = test, linetype = test)) +
    geom_hline(yintercept = 0.05, linetype = 3, colour = "grey50") +
    geom_line(linewidth = 1.2) + geom_point(size = 3) + facet_wrap(~panel) +
    scale_colour_manual(values = c(col_sarnb, col_zinb, col_zisar, "#7f2704")) +
    scale_shape_manual(values = c(15, 17, 16, 18)) +
    scale_linetype_manual(values = c(2, 2, 1, 1)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = expression("true " * rho), y = "rejection rate at 5%",
         title = "How often each test detects spatial dependence (200 replications per point)") +
    theme_deck() + guides(colour = guide_legend(nrow = 2))
  save(pP, "fig_power_moran.pdf", 10, 4.8)
}
