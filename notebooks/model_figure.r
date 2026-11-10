# =============================================================================
#  Figure of the model section: the two components of ZISAR on a synthetic
#  region of 300 units (same design as the workshop slides, slides/figures.R §8)
#  (a) observed counts, (b) the two kinds of zeros, (c) the latent intensity
#  estimated by ZISAR, (d) the posterior probability that a zero is structural
#  Output: ../tex/figures/model_components.pdf
# =============================================================================

source("zisar.r")
suppressPackageStartupMessages({ library(Matrix); library(ggplot2); library(deldir); library(patchwork) })

# ---- synthetic region: Voronoi polygons around 5 urban clusters --------------
set.seed(31)
n_m <- 300
centres <- cbind(runif(5, 0.15, 0.85), runif(5, 0.15, 0.85))
urban <- do.call(rbind, lapply(1:5, function(k) cbind(rnorm(25, centres[k, 1], 0.05), rnorm(25, centres[k, 2], 0.05))))
rural <- cbind(runif(n_m - nrow(urban)), runif(n_m - nrow(urban)))
pts <- rbind(urban, rural); pts <- pmin(pmax(pts, 0.005), 0.995)
hubs <- centres[1:4, ]
d_h  <- apply(pts, 1, function(p) min(sqrt((hubs[, 1] - p[1])^2 + (hubs[, 2] - p[2])^2)))
d_c  <- apply(pts, 1, function(p) min(sqrt((centres[, 1] - p[1])^2 + (centres[, 2] - p[2])^2)))
logpop <- as.numeric(scale(-d_c + rnorm(n_m, 0, 0.25)))     # count covariate: size
remote <- as.numeric(scale(d_h))                            # remoteness: both processes
D <- as.matrix(dist(pts)); diag(D) <- Inf
Wm <- t(apply(D, 1, function(r) { w <- numeric(n_m); w[order(r)[1:6]] <- 1 / 6; w }))   # 6 nearest neighbours
rho_m <- 0.5
u_m <- rnorm(n_m, 0, 0.5)
eta_m <- as.numeric(solve(diag(n_m) - rho_m * Wm, 0.8 + 0.7 * logpop - 0.2 * remote + u_m))
pi_m <- plogis(-1.6 + 1.4 * remote)
s_m <- rbinom(n_m, 1, pi_m)
y_m <- ifelse(s_m == 1, 0L, rnbinom(n_m, size = 2, mu = exp(eta_m)))

# ---- ZISAR fit ----------------------------------------------------------------
X <- cbind(1, logpop, remote); Z <- cbind(1, remote)
fit <- zisar_fit(y_m, X, Z, Matrix(Wm, sparse = TRUE), silent = TRUE)
pi_hat <- zisar_predict(fit, type = "zero")
mu_hat <- zisar_predict(fit, type = "count")
kap <- fit$params$kappa
post <- ifelse(y_m == 0, pi_hat / (pi_hat + (1 - pi_hat) * dnbinom(0, size = kap, mu = mu_hat)), NA)
z0 <- y_m == 0
auc <- function(s, lab) { r <- rank(s); n1 <- sum(lab); (sum(r[lab]) - n1 * (n1 + 1) / 2) / (n1 * sum(!lab)) }
cat(sprintf("zeros %d (structural %d, sampling %d) | rho_hat %.2f (true 0.5) | posterior: structural %.2f, sampling %.2f, AUC %.2f | cor(eta_hat, eta) %.2f\n",
            sum(z0), sum(s_m == 1), sum(z0 & s_m == 0), fit$params$rho,
            mean(post[s_m == 1]), mean(post[z0 & s_m == 0]), auc(post[z0], s_m[z0] == 1), cor(log(mu_hat), eta_m)))

# ---- maps -----------------------------------------------------------------------
vor <- deldir(pts[, 1], pts[, 2], rw = c(0, 1, 0, 1))
tiles <- tile.list(vor)
poly <- do.call(rbind, lapply(seq_along(tiles), function(i) data.frame(id = i, x = tiles[[i]]$x, y = tiles[[i]]$y)))

blue <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
col_samp <- "#2a78d6"; col_struct <- "#eb6834"; col_pos <- "#eeeeec"
th <- theme_void(base_size = 9) +
  theme(plot.title = element_text(size = 9.5, face = "bold", hjust = 0, margin = margin(b = 3)),
        legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 8),
        legend.key.size = unit(0.32, "cm"), legend.key.width = unit(0.55, "cm"),
        plot.margin = margin(4, 6, 4, 6))
map_base <- function(fillvar, title) {
  pd <- poly; pd$fill <- fillvar[pd$id]
  ggplot(pd, aes(x, y, group = id, fill = fill)) +
    geom_polygon(colour = "white", linewidth = 0.2) + coord_equal(expand = FALSE) +
    labs(title = title) + th
}

y_cls <- cut(y_m, c(-0.5, 0.5, 2.5, 5.5, 10.5, 20.5, Inf), labels = c("0", "1-2", "3-5", "6-10", "11-20", "> 20"))
pa <- map_base(y_cls, "(a) Observed counts") +
  scale_fill_manual(values = c(col_pos, blue[c(1, 3, 4, 5, 7)]), drop = FALSE) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

kind <- factor(ifelse(y_m > 0, "count > 0", ifelse(s_m == 1, "structural zero", "sampling zero")),
               levels = c("count > 0", "sampling zero", "structural zero"))
pb <- map_base(kind, "(b) Two kinds of zeros (true)") +
  scale_fill_manual(values = c(col_pos, col_samp, col_struct)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

eta_hat <- log(mu_hat)
lim <- as.numeric(quantile(eta_hat, c(0.02, 0.98)))
pc <- map_base(pmin(pmax(eta_hat, lim[1]), lim[2]), "(c) Estimated latent intensity (log scale)") +
  scale_fill_gradientn(colours = blue, breaks = pretty(lim, 4)) +
  guides(fill = guide_colourbar(barheight = unit(0.25, "cm"), barwidth = unit(3.5, "cm")))

post_fill <- ifelse(z0, post, NA)
pd_ <- map_base(post_fill, "(d) Estimated P(zero is structural)") +
  scale_fill_gradient(low = "#fde3d5", high = "#9c3a12", limits = c(0, 1), breaks = c(0, 0.5, 1),
                      na.value = col_pos) +
  guides(fill = guide_colourbar(barheight = unit(0.25, "cm"), barwidth = unit(3.5, "cm")))

fig <- (pa | pb) / (pc | pd_)
ggsave("../tex/figures/model_components.pdf", fig, width = 6.3, height = 6.9, device = pdf)
