# =============================================================================
#  Monte Carlo figure and appendix tables from mc_results_1000.rds
#  (written by mc_run_1000.r: R = 1000 per scenario, rho in {0, .15, .3, .45, .6})
#  Output: ../tex/figures/mc_results.pdf and mc_1000_tables.txt (LaTeX rows)
# =============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

mc <- readRDS("mc_results_1000.rds")
BETA <- c(0.5, 0.8, -0.5); KAPPA <- 2
NS <- c(100, 225, 400); RHOS <- sort(unique(mc$rho_true))
MODELS <- c("ZISAR", "SAR-NB", "ZINB")
cols   <- c(ZISAR = "#2a78d6", `SAR-NB` = "#eb6834", ZINB = "#1baf7a")
shapes <- c(ZISAR = 16, `SAR-NB` = 17, ZINB = 15)
ltys   <- c(ZISAR = "solid", `SAR-NB` = "22", ZINB = "13")

# ZISAR summaries use converged fits with kappa below 100 (in practice the
# excluded fits sit at the upper bound kappa = 500); the competitors use
# their converged fits
ok <- list(ZISAR    = mc$zisar_conv %in% TRUE & mc$zisar_kappa < 100,
           `SAR-NB` = mc$sarnb_conv %in% TRUE,
           ZINB     = mc$zinb_conv %in% TRUE)
pref <- c(ZISAR = "zisar", `SAR-NB` = "sarnb", ZINB = "zinb")

long <- do.call(rbind, lapply(MODELS, function(m) {
  d <- mc[ok[[m]], ]
  rho <- if (m == "ZINB") NA_real_ else d[[paste0(pref[m], "_rho")]]
  data.frame(model = m, n = d$n, rho_true = d$rho_true, rho = rho,
             b1 = d[[paste0(pref[m], "_b1")]], b2 = d[[paste0(pref[m], "_b2")]],
             kappa = d[[paste0(pref[m], "_kappa")]], ll = d[[paste0(pref[m], "_ll")]])
}))
long$model <- factor(long$model, levels = MODELS)

sm <- do.call(rbind, lapply(split(long, list(long$model, long$n, long$rho_true), drop = TRUE), function(d) {
  r0 <- d$rho_true[1]
  data.frame(model = d$model[1], n = d$n[1], rho_true = r0, reps = nrow(d),
             rho_mean = mean(d$rho), rho_lo = unname(quantile(d$rho, 0.1, na.rm = TRUE)),
             rho_hi = unname(quantile(d$rho, 0.9, na.rm = TRUE)),
             rho_bias = mean(d$rho - r0), rho_rmse = sqrt(mean((d$rho - r0)^2)),
             b1_bias = mean(d$b1 - BETA[2]), b1_rmse = sqrt(mean((d$b1 - BETA[2])^2)),
             b2_bias = mean(d$b2 - BETA[3]), b2_rmse = sqrt(mean((d$b2 - BETA[3])^2)),
             kappa_medbias = median(d$kappa) - KAPPA, ll_med = median(d$ll))
}))
sm$panel <- factor(paste0("n = ", sm$n), levels = paste0("n = ", NS))

# ---- figure ------------------------------------------------------------------
fam <- ""  # Helvetica: base-14 font, no embedding needed
th <- theme_minimal(base_size = 9.5, base_family = fam) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
        strip.text = element_text(size = 9.5, face = "bold", hjust = 0),
        axis.title = element_text(size = 9.5), legend.position = "bottom",
        legend.title = element_blank(), legend.key.width = unit(1.6, "lines"),
        panel.spacing.x = unit(1, "lines"), plot.tag = element_text(size = 10, face = "bold"))
sx <- scale_x_continuous(breaks = RHOS, labels = c("0", ".15", ".3", ".45", ".6"))
enc <- list(scale_colour_manual(values = cols, drop = FALSE), scale_fill_manual(values = cols, drop = FALSE),
            scale_shape_manual(values = shapes, drop = FALSE), scale_linetype_manual(values = ltys, drop = FALSE))

d1 <- subset(sm, model != "ZINB")
p1 <- ggplot(d1, aes(rho_true, rho_mean, colour = model, fill = model, shape = model, linetype = model)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_ribbon(aes(ymin = rho_lo, ymax = rho_hi), alpha = 0.15, colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
  facet_wrap(~panel) + sx + enc +
  labs(x = NULL, y = expression("estimated " * rho), tag = "(a)") + th + theme(legend.position = "none")

# RMSE of beta_1 relative to ZISAR in the same scenario (ZISAR = 1): the three
# estimators see the same datasets, so the ratio removes the noise they share
ref <- subset(sm, model == "ZISAR")[, c("n", "rho_true", "b1_rmse")]
names(ref)[3] <- "b1_ref"
d2 <- merge(sm, ref, by = c("n", "rho_true"))
d2$b1_rel <- d2$b1_rmse / d2$b1_ref
p2 <- ggplot(d2, aes(rho_true, b1_rel, colour = model, shape = model, linetype = model)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
  facet_wrap(~panel) + sx + enc +
  labs(x = expression("true " * rho), y = expression("relative RMSE of " * beta[1]), tag = "(b)") + th +
  theme(strip.text = element_blank())

fig <- p1 / p2
ggsave("../tex/figures/mc_results.pdf", fig, width = 6.3, height = 4.9, device = pdf)

# ---- appendix tables: LaTeX rows for rho in {0, 0.3, 0.6} ---------------------
fmt <- function(x, d = 3) ifelse(is.na(x), "--", sprintf(paste0("$%.", d, "f$"), x))
row_of <- function(m, v, d = 3) {
  vals <- unlist(lapply(c(0, 0.3, 0.6), function(r) sapply(NS, function(nn)
    sm[sm$model == m & sm$n == nn & abs(sm$rho_true - r) < 1e-9, v])))
  paste0("  & ", sprintf("%-8s", m), " & ", paste(fmt(vals, d), collapse = " & "), " \\\\")
}
out <- c()
for (v in c("b1_bias", "b1_rmse", "b2_bias", "b2_rmse")) {
  out <- c(out, paste("%", v), sapply(MODELS, row_of, v = v))
}
for (v in c("rho_bias", "rho_rmse")) out <- c(out, paste("%", v), sapply(c("SAR-NB", "ZISAR"), row_of, v = v))
out <- c(out, "% kappa_medbias", sapply(MODELS, row_of, v = "kappa_medbias"))
out <- c(out, "% ll_med", sapply(MODELS, row_of, v = "ll_med", d = 0))

# replications kept for ZISAR, convergence, boundary share and time
cv <- do.call(rbind, lapply(split(mc, list(mc$n, mc$rho_true)), function(d) {
  cz <- d$zisar_conv %in% TRUE
  data.frame(n = d$n[1], rho = d$rho_true[1], zeros = mean(d$zeros), conv = mean(cz),
             bound = mean(d$zisar_kappa[cz] >= 100), kept = sum(cz & d$zisar_kappa < 100),
             time = mean(d$zisar_time, na.rm = TRUE))
}))
out <- c(out, "", capture.output(print(cv[order(cv$rho, cv$n), ], row.names = FALSE, digits = 3)))
writeLines(out, "mc_1000_tables.txt")
cat(out, sep = "\n")
