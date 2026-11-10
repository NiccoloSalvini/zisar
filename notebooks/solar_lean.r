# =============================================================================
#  ZISAR application: domestic solar PV adoption across English LSOAs
#  (specification of the paper: public area statistics only)
#
#  y_i : domestic PV installations of up to 10 kW commissioned under the
#        Feed-in Tariff after the tariff reform of 8 February 2016 and before
#        the scheme closed to new entrants (31 March 2019), in unit i
#        (n = 32,322 harmonised 2001/2011 LSOAs, see solar_data.r)
#  X   : count-process covariates (2011 Census, IoD 2019, VOA 2015, DESNZ
#        2010, northing), including the covariates of Z
#  Z   : structural-zero covariates: housing form and age, density, tenure,
#        communal and vacant dwellings
#  W   : queen contiguity, row-standardised
#
#  Six specifications: NB, ZINB, SAR-NB (instrumented lag), ZINB-RE (rho = 0),
#  SAR-NB-RE (pi = 0) and ZISAR; direct, indirect and total effects of ZISAR.
#  The ZISAR fit is the "nobar" specification of solar_application.r and is
#  read from its cache when present.
#
#  Run solar_data.r first. Writes solar_lean_results.txt and the figure
#  ../tex/figures/solar_lean_maps.png
# =============================================================================

source("zisar.r")
suppressPackageStartupMessages({
  library(pscl); library(spdep); library(sf); library(ggplot2); library(patchwork)
})
set.seed(2026)

solar <- readRDS("../data/solar/solar_lsoa.rds")
Ws    <- readRDS("../data/solar/W_lsoa.rds")
W     <- Ws$queen
n     <- nrow(solar)
y     <- solar$pv_dom_late

std <- function(x) (x - mean(x)) / sd(x)
cov_x <- with(solar, data.frame(
  households = std(log(households)),    # log number of households
  detached   = std(sh_detached),        # detached-house share
  social     = std(sh_social),          # social-rented share
  income_dep = std(income_dep),         # IoD income deprivation rate
  managerial = std(sh_nssec12),         # NS-SEC 1-2: managerial and professional
  age65      = std(sh_age65),           # share aged 65+
  children   = std(sh_children),        # households with dependent children
  offgas     = std(sh_offgas),          # households without gas central heating
  band_eh    = std(sh_band_eh),         # Council Tax bands E-H (VOA 2015)
  elec       = std(log(elec_kwh)),      # log domestic kWh per meter, 2010
  north      = std(north / 1000)))      # northing (km): solar resource gradient
cov_z <- with(solar, data.frame(
  flat       = std(sh_flat),            # share of household spaces in flats
  density    = std(log(density)),       # log residents per km^2
  private    = std(sh_private),         # private-rented share
  communal   = std(sh_communal),        # residents in communal establishments
  empty      = std(sh_empty),           # household spaces with no usual residents
  pre1919    = std(sh_pre1919)))        # dwellings built before 1919 (VOA 2015)
# Housing form also shifts the intensity among eligible units: the Z covariates
# enter X too (same names, so zisar_effects adds the inflation channel)
X <- cbind(intercept = 1, as.matrix(cov_x), as.matrix(cov_z))
Z <- cbind(intercept = 1, as.matrix(cov_z))
p <- ncol(X); q <- ncol(Z)
df <- data.frame(y = y, cov_x, cov_z)
f_count <- paste(c(colnames(cov_x), colnames(cov_z)), collapse = " + ")
f_zero  <- paste(colnames(cov_z), collapse = " + ")
lw <- nb2listw(Ws$nb_queen, style = "W")

# -----------------------------------------------------------------------------
# 1. Benchmarks without random effects
# -----------------------------------------------------------------------------
fit_nb <- MASS::glm.nb(as.formula(paste("y ~", f_count)), data = df)
fit_zinb <- zeroinfl(as.formula(paste("y ~", f_count, "|", f_zero)),
                     data = df, dist = "negbin")
# SAR-NB: NB GLM with the spatial lag W y instrumented by [X, WX, W^2X]
WX  <- as.matrix(W %*% X[, -1]); W2X <- as.matrix(W %*% WX)
H   <- cbind(X, WX, W2X)
Wy  <- as.numeric(W %*% y)
df$Wy_hat <- as.numeric(H %*% solve(crossprod(H), crossprod(H, Wy)))
fit_sarnb <- MASS::glm.nb(as.formula(paste("y ~", f_count, "+ Wy_hat")), data = df)

# -----------------------------------------------------------------------------
# 2. ZISAR and its two nested restrictions (sparse log-determinant), cached
# -----------------------------------------------------------------------------
fits_file <- "../data/solar/solar_lean_fits.rds"
fits <- if (file.exists(fits_file)) readRDS(fits_file) else list()
fit_one <- function(fixed = NULL, init = NULL) {
  t0 <- Sys.time()
  f <- zisar_fit(y, X, Z, W, fixed = fixed, init = init)
  f$elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  f$tmb_obj <- NULL                     # external pointers do not survive the fork
  f
}
# Starting values from a fitted model, with kappa near the upper bound where
# the fits of this outcome end
warm <- function(f) list(beta = f$params$beta, gamma = f$params$gamma,
                         rho_raw = atanh(f$params$rho), log_kappa = log(499),
                         log_sigma = log(f$params$sigma), eta = f$eta_hat)
.zisar_compile()
if (is.null(fits$zisar)) {
  old <- "../data/solar/solar_fits.rds"
  fits$zisar <- if (file.exists(old)) readRDS(old)$nobar
  if (is.null(fits$zisar)) {
    fits$zisar <- fit_one()
    if (fits$zisar$convergence != 0) fits$zisar <- fit_one(init = warm(fits$zisar))
  }
  stopifnot(isTRUE(all.equal(unname(fits$zisar$X), unname(X))),
            isTRUE(all.equal(unname(fits$zisar$Z), unname(Z))),
            identical(as.numeric(fits$zisar$y), as.numeric(y)))
  saveRDS(fits, fits_file)
}
nested <- list(
  zinbre  = list(rho_raw = 0),                       # rho = 0
  sarnbre = list(gamma = c(-20, rep(0, q - 1))))     # pi = 0
todo <- setdiff(names(nested), names(fits))
if (length(todo)) {
  fits[todo] <- parallel::mclapply(todo, function(nm) {
    f <- fit_one(fixed = nested[[nm]], init = warm(fits$zisar))
    # a false-convergence stop is restarted from its stopping point
    if (f$convergence != 0) { e <- f$elapsed; f <- fit_one(fixed = nested[[nm]], init = warm(f)); f$elapsed <- f$elapsed + e }
    f
  }, mc.cores = min(2, length(todo)), mc.preschedule = FALSE)
  stopifnot(all(sapply(fits[todo], function(f) is.list(f) && is.finite(f$loglik))))
  saveRDS(fits, fits_file)
}
fit_zisar <- fits$zisar; fit_zinbre <- fits$zinbre; fit_sarnbre <- fits$sarnbre
lr_rho <- 2 * (fit_zisar$loglik - fit_zinbre$loglik)
lr_pi  <- 2 * (fit_zisar$loglik - fit_sarnbre$loglik)

# -----------------------------------------------------------------------------
# 3. Effects decomposition
# -----------------------------------------------------------------------------
if (is.null(fits$eff)) {
  t_eff <- system.time(fits$eff <- zisar_effects(fit_zisar, B = 500))["elapsed"]
  fits$eff_elapsed <- as.numeric(t_eff)
  saveRDS(fits, fits_file)
}
eff <- fits$eff

# -----------------------------------------------------------------------------
# 4. Expected number of zeros (marginal: the random-effects models integrate
#    the latent field by simulation from the fitted distribution)
# -----------------------------------------------------------------------------
zero_nb <- function(mu, size) sum(dnbinom(0, size = size, mu = mu))
exp_nb    <- zero_nb(fitted(fit_nb), fit_nb$theta)
exp_sarnb <- zero_nb(fitted(fit_sarnb), fit_sarnb$theta)
pz_zinb <- predict(fit_zinb, type = "zero"); mu_zinb <- predict(fit_zinb, type = "count")
exp_zinb <- sum(pz_zinb + (1 - pz_zinb) * dnbinom(0, size = fit_zinb$theta, mu = mu_zinb))
marginal_zeros <- function(fit, R = 200) {
  A <- Diagonal(fit$n) - fit$params$rho * as(fit$W, "CsparseMatrix")
  luA <- lu(A)
  Xb <- as.numeric(fit$X %*% fit$params$beta)
  pi_i <- plogis(as.numeric(fit$Z %*% fit$params$gamma))
  kap <- fit$params$kappa
  acc <- 0
  for (r in seq_len(R)) {
    eta <- as.numeric(solve(luA, Xb + rnorm(fit$n, 0, fit$params$sigma)))
    acc <- acc + sum(pi_i + (1 - pi_i) * dnbinom(0, size = kap, mu = exp(eta)))
  }
  acc / R
}
set.seed(2026)   # same draws whether the fits above were cached or not
exp_re <- sapply(list(`ZINB-RE` = fit_zinbre, `SAR-NB-RE` = fit_sarnbre, ZISAR = fit_zisar),
                 marginal_zeros)

# Structural zeros of ZISAR
pi_hat <- zisar_predict(fit_zisar, type = "zero")
mu_hat <- zisar_predict(fit_zisar, type = "count")
kap <- fit_zisar$params$kappa
post_struct <- ifelse(y == 0, pi_hat / (pi_hat + (1 - pi_hat) * (kap / (kap + mu_hat))^kap), 0)
london <- substr(solar$lad19cd, 1, 3) == "E09"
# Validation against adoption before the reform, which the model does not use:
# mean posterior probability by number of pre-reform installations, and AUC of
# the posterior as a classifier of the zeros with none over the whole scheme
pre <- solar$pv_dom - y
val_a <- post_struct[y == 0 & pre == 0]; val_b <- post_struct[y == 0 & pre > 0]
validation <- c(none = mean(val_a), `1-4` = mean(post_struct[y == 0 & pre %in% 1:4]),
                `5+` = mean(post_struct[y == 0 & pre >= 5]),
                AUC = mean(outer(val_a, val_b, ">") + 0.5 * outer(val_a, val_b, "==")))
flat_dec <- cut(solar$sh_flat, quantile(solar$sh_flat, 0:10 / 10), include.lowest = TRUE, labels = FALSE)

# Moran's I of the counts and of the Pearson residuals of the benchmarks
mi <- function(x) { m <- moran.test(x, lw, randomisation = FALSE)
                    c(I = unname(m$estimate[1]), z = unname(m$statistic)) }
moran <- rbind(y = mi(y), NB = mi(residuals(fit_nb, type = "pearson")),
               ZINB = mi(residuals(fit_zinb, type = "pearson")),
               `SAR-NB` = mi(residuals(fit_sarnb, type = "pearson")))

# -----------------------------------------------------------------------------
# 5. Results
# -----------------------------------------------------------------------------
ll <- function(f) as.numeric(logLik(f)); k_of <- function(f) attr(logLik(f), "df")
info <- rbind(
  NB          = c(ll(fit_nb), k_of(fit_nb), exp_nb),
  ZINB        = c(ll(fit_zinb), k_of(fit_zinb), exp_zinb),
  `SAR-NB`    = c(ll(fit_sarnb), k_of(fit_sarnb), exp_sarnb),
  `ZINB-RE`   = c(fit_zinbre$loglik, length(fit_zinbre$coefficients), exp_re["ZINB-RE"]),
  `SAR-NB-RE` = c(fit_sarnbre$loglik, length(fit_sarnbre$coefficients), exp_re["SAR-NB-RE"]),
  ZISAR       = c(fit_zisar$loglik, length(fit_zisar$coefficients), exp_re["ZISAR"]))
info <- cbind(logLik = info[, 1], k = info[, 2], AIC = 2 * info[, 2] - 2 * info[, 1],
              BIC = log(n) * info[, 2] - 2 * info[, 1], exp_zeros = info[, 3])
named <- function(f) { ct <- f$coef_table
  free_b <- sum(grepl("^beta_", rownames(ct))); free_g <- sum(grepl("^gamma_", rownames(ct)))
  rownames(ct)[grepl("^beta_", rownames(ct))] <- colnames(X)[seq_len(free_b)]
  rownames(ct)[grepl("^gamma_", rownames(ct))] <- paste0("zero:", colnames(Z)[seq_len(free_g)])
  ct }

sink("solar_lean_results.txt")
cat("Domestic PV up to 10 kW installed under the FiT, 8 Feb 2016 - 31 Mar 2019, per unit, England\n")
cat("n =", n, " installations =", sum(y), " zeros =", sum(y == 0),
    sprintf("(%.2f%%)", 100 * mean(y == 0)), " ones =", sum(y == 1),
    " mean =", round(mean(y), 2), " var =", round(var(y), 1),
    " median =", median(y), " max =", max(y), "\n")
cat("W queen: mean neighbours", round(mean(card(lw$neighbours)), 2), "\n")
cat("\n--- Descriptives (raw covariates) ---\n")
print(round(t(sapply(solar[, c("pv_dom_late", "households", "sh_detached", "sh_social", "income_dep",
  "sh_nssec12", "sh_age65", "sh_children", "sh_offgas", "sh_band_eh", "elec_kwh", "north",
  "sh_flat", "density", "sh_private", "sh_communal", "sh_empty", "sh_pre1919")],
  function(v) c(mean = mean(v), sd = sd(v), min = min(v), median = median(v), max = max(v)))), 3))
cat("\n--- NB ---\n"); print(round(summary(fit_nb)$coefficients, 4)); cat("theta =", fit_nb$theta, "SE", fit_nb$SE.theta, "\n")
cat("\n--- ZINB ---\n"); print(lapply(summary(fit_zinb)$coefficients, round, 4)); cat("theta =", fit_zinb$theta, "\n")
cat("\n--- SAR-NB ---\n"); print(round(summary(fit_sarnb)$coefficients, 4)); cat("theta =", fit_sarnb$theta, "\n")
cat("\n--- ZISAR (convergence", fit_zisar$convergence, ", max_grad", signif(fit_zisar$max_grad, 3),
    ", elapsed", round(fit_zisar$elapsed, 1), "s) ---\n")
print(round(named(fit_zisar), 4)); print(fit_zisar$transformed)
cat("\n--- ZINB-RE (rho = 0; convergence", fit_zinbre$convergence, ", max_grad", signif(fit_zinbre$max_grad, 3), ") ---\n")
print(round(named(fit_zinbre), 4)); print(unlist(fit_zinbre$params[c("kappa", "sigma")]))
cat("\n--- SAR-NB-RE (pi = 0; convergence", fit_sarnbre$convergence, ", max_grad", signif(fit_sarnbre$max_grad, 3), ") ---\n")
print(round(named(fit_sarnbre), 4)); print(fit_sarnbre$transformed)
cat("\n--- Fit of the six specifications (observed zeros", sum(y == 0), ") ---\n"); print(round(info))
cat(sprintf("\nLR rho = 0 (ZISAR vs ZINB-RE): %.1f; Wald z on atanh(rho): %.1f\n", lr_rho, fit_zisar$coef_table["rho_raw", "z"]))
cat(sprintf("LR pi = 0 (ZISAR vs SAR-NB-RE): %.1f on %d restrictions\n", lr_pi, q))
cat("\n--- Effects on E[y] per 1 SD (average over units; 'zero' = inflation channel, part of direct) ---\n")
for (e in eff[-1]) cat(sprintf("%-11s direct %7.4f (%.4f) indirect %7.4f (%.4f) total %7.4f (%.4f) zero %7.4f\n",
  e$name, e$avg_direct, e$se_direct, e$avg_indirect, e$se_indirect, e$avg_total, e$se_total, e$avg_zero))
ind_share <- sapply(eff[-1], function(e) e$avg_indirect / e$avg_total)
cat("indirect share of the total effect, covariates of the count process only:",
    round(range(ind_share[!colnames(X)[-1] %in% colnames(Z)]), 3), "\n")
cat("mean of 1/(1 - rho):", round(1 / (1 - fit_zisar$params$rho), 3), " elapsed", round(fits$eff_elapsed, 1), "s\n")
cat("\n--- Structural zeros (ZISAR) ---\n")
cat("mean pi_hat:", round(mean(pi_hat), 4), " expected structural zeros:", round(sum(pi_hat), 1),
    " share of observed zeros:", round(sum(pi_hat) / sum(y == 0), 3),
    " posterior sum among observed zeros:", round(sum(post_struct), 1),
    " zeros with P(structural) > 0.5:", sum(post_struct > 0.5), "\n")
cat("London: units", sum(london), " zeros", sum(y[london] == 0), " expected structural", round(sum(pi_hat[london]), 1),
    " (share of all structural:", round(sum(pi_hat[london]) / sum(pi_hat), 3), ")\n")
cat("pi_hat by quintile of the flat share:", round(tapply(pi_hat, cut(solar$sh_flat, quantile(solar$sh_flat, 0:5 / 5), include.lowest = TRUE), mean), 4), "\n")
cat("units with pi_hat > 0.1:", sum(pi_hat > 0.1), " > 0.5:", sum(pi_hat > 0.5),
    " of which in London:", sum(pi_hat > 0.5 & london), "\n")
cat("pi_hat quantiles (50, 90, 99%):", signif(quantile(pi_hat, c(.5, .9, .99)), 3), "\n")
cat("share of units with y = 0: top decile of flats", round(mean(y[flat_dec == 10] == 0), 3),
    " below the median", round(mean(y[flat_dec <= 5] == 0), 3),
    " London", round(mean(y[london] == 0), 3), " elsewhere", round(mean(y[!london] == 0), 3), "\n")
cat("zero units by pre-reform installations: none", sum(y == 0 & pre == 0), " 1-4", sum(y == 0 & pre %in% 1:4),
    " 5+", sum(y == 0 & pre >= 5), sprintf("(%.1f%% of the zero units)", 100 * mean(pre[y == 0] >= 5)), "\n")
cat("posterior P(structural | y = 0) by pre-reform installations, and AUC for the zeros with none:\n")
print(round(validation, 3))
cat("\n--- Moran's I (queen) ---\n"); print(round(moran, 4))
sink()

# -----------------------------------------------------------------------------
# 6. Maps: adoption rate (England), observed zeros and structural-zero
#    probability (Greater London, where most structural zeros are)
# -----------------------------------------------------------------------------
map <- readRDS("../data/solar/lsoa_map.rds")
idx <- match(map$unit, solar$unit)
map$rate <- 1000 * y[idx] / solar$households[idx]
map$pi   <- pi_hat[idx]
map$zero <- y[idx] == 0
lon <- london[idx]
theme_map <- theme_void(base_size = 8) +
  theme(legend.position = "bottom", legend.key.width = unit(0.6, "cm"),
        legend.key.height = unit(0.22, "cm"), legend.title.position = "top",
        legend.title = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5))
g1 <- ggplot(map) + geom_sf(aes(fill = pmin(rate, 10)), colour = NA) +
  scale_fill_viridis_c(option = "magma", name = "installations per 1,000 households",
                       breaks = c(0, 5, 10), labels = c("0", "5", "10+")) +
  ggtitle("(a) England") + theme_map
g2 <- ggplot(map[lon, ]) + geom_sf(aes(fill = zero), colour = NA) +
  scale_fill_manual(values = c(`FALSE` = "grey85", `TRUE` = "#b2182b"),
                    labels = c("at least one", "none"), name = "installations") +
  ggtitle("(b) London: observed zeros") + theme_map
g3 <- ggplot(map[lon, ]) + geom_sf(aes(fill = pi), colour = NA) +
  scale_fill_viridis_c(option = "viridis", trans = "sqrt", limits = c(0, 1),
                       breaks = c(0.01, 0.1, 0.3, 0.6, 1),
                       name = "structural-zero probability") +
  ggtitle(expression(paste("(c) London: ", hat(pi)[i]))) + theme_map
dir.create("../tex/figures", showWarnings = FALSE)
ggsave("../tex/figures/solar_lean_maps.png", g1 + g2 + g3, width = 7.2, height = 3.6, dpi = 300, bg = "white")
