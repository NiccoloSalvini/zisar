# =============================================================================
#  ZISAR application: domestic solar PV adoption across English LSOAs
#
#  y_i : domestic PV installations commissioned under the Feed-in Tariff after
#        the tariff reform of 8 February 2016 and before the scheme closed to
#        new entrants (31 March 2019), in LSOA i (n = 32,322 harmonised
#        2001/2011 LSOAs, see solar_data.r)
#  X   : count-process covariates (2011 Census, IoD 2019, VOA 2015, DESNZ
#        2010, location), including the housing covariates of Z
#  Z   : structural-zero covariates (housing form and age, density, tenure,
#        communal and vacant dwellings)
#  W   : queen contiguity, row-standardised (k-NN in the sensitivity check)
#  The cumulative 2010-2019 count is used as a robustness check.
#
#  Run solar_data.r first. Writes solar_results.txt and solar_results.rds,
#  and the figure ../tex/figures/solar_maps.png
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
  social     = std(sh_social),          # social-rented share (landlord programmes)
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
# Housing form also shifts the intensity among eligible LSOAs: the Z covariates
# enter X too (same names, so zisar_effects adds the inflation channel)
X <- cbind(intercept = 1, as.matrix(cov_x), as.matrix(cov_z))
Z <- cbind(intercept = 1, as.matrix(cov_z))
p <- ncol(X); q <- ncol(Z)
df <- data.frame(y = y, cov_x, cov_z)
f_count <- paste(c(colnames(cov_x), colnames(cov_z)), collapse = " + ")
f_zero  <- paste(colnames(cov_z), collapse = " + ")

lw <- nb2listw(Ws$nb_queen, style = "W")
timing <- list()

# -----------------------------------------------------------------------------
# 1. Benchmarks
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
# 2. ZISAR, its nested restrictions, the alternative W and the cumulative
#    outcome (same template, sparse log-determinant). One fit takes about
#    40 min on a single core at n = 32,322, so the six fits run as parallel
#    processes and are cached in solar_fits.rds.
# -----------------------------------------------------------------------------
fits_file <- "../data/solar/solar_fits.rds"
if (file.exists(fits_file)) {
  fits <- readRDS(fits_file)
} else {
  specs <- list(
    zisar   = list(y = y, W = W),
    # rho = 0: aspatial ZINB with Gaussian random effects (LR test of rho)
    zinbre  = list(y = y, W = W, fixed = list(rho_raw = 0)),
    # pi = 0: spatial NB with random effects, no zero-inflation
    sarnbre = list(y = y, W = W, fixed = list(gamma = c(-20, rep(0, q - 1)))),
    k6      = list(y = y, W = Ws$k6),
    k10     = list(y = y, W = Ws$k10),
    cum     = list(y = solar$pv_dom, W = W))
  .zisar_compile()
  fits <- parallel::mclapply(specs, function(s) {
    t0 <- Sys.time()
    f <- zisar_fit(s$y, X, Z, s$W, fixed = s$fixed)
    f$elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    f$tmb_obj <- NULL                     # external pointers do not survive the fork
    f
  }, mc.cores = length(specs), mc.preschedule = FALSE)
  stopifnot(all(sapply(fits, inherits, "zisar")))
  # A fit of the post-reform outcome that stops with a false-convergence flag
  # is restarted from its stopping point with kappa near the upper bound,
  # where all the converged fits of this application end
  for (nm in setdiff(names(fits), "cum")) if (fits[[nm]]$convergence != 0) {
    f0 <- fits[[nm]]
    t0 <- Sys.time()
    f <- zisar_fit(specs[[nm]]$y, X, Z, specs[[nm]]$W, fixed = specs[[nm]]$fixed,
                   init = list(beta = f0$params$beta, gamma = f0$params$gamma,
                               rho_raw = atanh(f0$params$rho), log_kappa = log(499),
                               log_sigma = log(f0$params$sigma), eta = f0$eta_hat))
    f$elapsed <- f0$elapsed + as.numeric(difftime(Sys.time(), t0, units = "secs"))
    f$restarted <- TRUE
    f$tmb_obj <- NULL
    fits[[nm]] <- f
  }
  saveRDS(fits, fits_file)
}
fit_zisar   <- fits$zisar
fit_zinbre  <- fits$zinbre
fit_sarnbre <- fits$sarnbre
timing$zisar <- fit_zisar$elapsed

lr_rho <- 2 * (fit_zisar$loglik - fit_zinbre$loglik)
wald_rho <- fit_zisar$coef_table["rho_raw", ]

# -----------------------------------------------------------------------------
# 3. Effects decomposition
# -----------------------------------------------------------------------------
timing$effects <- system.time(eff <- zisar_effects(fit_zisar, B = 500))["elapsed"]

# -----------------------------------------------------------------------------
# 4. Fit of the count distribution: observed vs expected frequencies
#    (ZISAR expected frequencies integrate the latent field by simulation
#     from the fitted SAR prior, i.e. they are marginal, like the others)
# -----------------------------------------------------------------------------
K <- 0:5
obs_freq <- sapply(K, function(k) sum(y == k))

nb_pmf <- function(k, mu, size) dnbinom(k, size = size, mu = mu)
exp_nb <- sapply(K, function(k) sum(nb_pmf(k, fitted(fit_nb), fit_nb$theta)))
exp_zinb <- colSums(predict(fit_zinb, type = "prob")[, K + 1])
exp_sarnb <- sapply(K, function(k) sum(nb_pmf(k, fitted(fit_sarnb), fit_sarnb$theta)))

zisar_marginal_freq <- function(fit, K, R = 200) {
  A <- Diagonal(fit$n) - fit$params$rho * as(fit$W, "CsparseMatrix")
  luA <- lu(A)
  Xb <- as.numeric(fit$X %*% fit$params$beta)
  pi_i <- plogis(as.numeric(fit$Z %*% fit$params$gamma))
  kap <- fit$params$kappa
  acc <- numeric(length(K))
  for (r in seq_len(R)) {
    eta <- as.numeric(solve(luA, Xb + rnorm(fit$n, 0, fit$params$sigma)))
    for (j in seq_along(K))
      acc[j] <- acc[j] + sum((K[j] == 0) * pi_i +
                             (1 - pi_i) * dnbinom(K[j], size = kap, mu = exp(eta)))
  }
  acc / R
}
exp_zisar <- zisar_marginal_freq(fit_zisar, K)
freq_tab <- rbind(observed = obs_freq, NB = exp_nb, ZINB = exp_zinb,
                  `SAR-NB` = exp_sarnb, ZISAR = exp_zisar)
colnames(freq_tab) <- K

# Posterior probability that an observed zero is structural (ZISAR, at the mode)
pi_hat <- zisar_predict(fit_zisar, type = "zero")
mu_hat <- zisar_predict(fit_zisar, type = "count")
p0_nb  <- (fit_zisar$params$kappa / (fit_zisar$params$kappa + mu_hat))^fit_zisar$params$kappa
post_struct <- ifelse(y == 0, pi_hat / (pi_hat + (1 - pi_hat) * p0_nb), 0)

# -----------------------------------------------------------------------------
# 5. Residual spatial autocorrelation (Pearson residuals)
# -----------------------------------------------------------------------------
pearson_zisar <- function(fit) {
  m <- zisar_predict(fit, type = "count"); pz <- zisar_predict(fit, type = "zero")
  v <- (1 - pz) * m * (1 + m / fit$params$kappa + pz * m)
  (fit$y - (1 - pz) * m) / sqrt(v)
}
moran <- list(
  y      = moran.test(y, lw),
  NB     = moran.test(residuals(fit_nb, type = "pearson"), lw),
  ZINB   = moran.test(residuals(fit_zinb, type = "pearson"), lw),
  SARNB  = moran.test(residuals(fit_sarnb, type = "pearson"), lw),
  ZISAR  = moran.test(pearson_zisar(fit_zisar), lw))

# -----------------------------------------------------------------------------
# 6. Sensitivity to W
# -----------------------------------------------------------------------------
sens <- fits[c("k6", "k10")]

# -----------------------------------------------------------------------------
# 6b. Robustness: cumulative count over the whole scheme (2010-2019)
# -----------------------------------------------------------------------------
y_cum  <- solar$pv_dom
df_cum <- transform(df, y = y_cum)
rob <- list(
  nb    = MASS::glm.nb(as.formula(paste("y ~", f_count)), data = df_cum),
  zinb  = zeroinfl(as.formula(paste("y ~", f_count, "|", f_zero)),
                   data = df_cum, dist = "negbin"),
  zisar = fits$cum)

# With counts in the hundreds the optimiser stops with a false-convergence flag.
# Restart twice from the stopping point, the second time with kappa held at its
# value there; the random-number stream is restored so that the draws used
# below are unaffected.
if (is.null(fits$cum_restarts)) {
  seed <- .Random.seed
  start_from <- function(f) list(beta = f$params$beta, gamma = f$params$gamma,
    rho_raw = atanh(f$params$rho), log_kappa = log(f$params$kappa),
    log_sigma = log(f$params$sigma), eta = f$eta_hat)
  r1 <- zisar_fit(y_cum, X, Z, W, init = start_from(fits$cum))
  i2 <- start_from(r1); i2$log_kappa <- NULL
  r2 <- zisar_fit(y_cum, X, Z, W, init = i2,
                  fixed = list(log_kappa = log(r1$params$kappa)))
  fits$cum_restarts <- lapply(list(r1, r2), `[`,
    c("params", "coef_table", "loglik", "AIC", "convergence"))
  saveRDS(fits, fits_file)
  assign(".Random.seed", seed, envir = globalenv())
}

# -----------------------------------------------------------------------------
# 7. Results
# -----------------------------------------------------------------------------
saveRDS(list(fit_nb = fit_nb, fit_zinb = fit_zinb, fit_sarnb = fit_sarnb,
             zisar = fit_zisar[setdiff(names(fit_zisar), c("tmb_obj", "W", "X", "Z"))],
             zinbre = fit_zinbre[c("params", "coef_table", "loglik", "AIC", "BIC")],
             sarnbre = fit_sarnbre[c("params", "coef_table", "loglik", "AIC", "BIC")],
             sens = lapply(sens, `[`, c("params", "coef_table", "transformed", "loglik", "AIC")),
             eff = eff, freq_tab = freq_tab, moran = moran, timing = timing,
             pi_hat = pi_hat, mu_hat = mu_hat, post_struct = post_struct,
             rob = list(nb = rob$nb, zinb = rob$zinb,
                        zisar = rob$zisar[c("params", "coef_table", "transformed",
                                            "loglik", "AIC", "BIC", "convergence")])),
        "solar_results.rds")

ll <- function(f) as.numeric(logLik(f))
k_of <- function(f) attr(logLik(f), "df")
info <- rbind(
  NB      = c(ll(fit_nb),  k_of(fit_nb)),
  ZINB    = c(ll(fit_zinb), k_of(fit_zinb)),
  `SAR-NB` = c(ll(fit_sarnb), k_of(fit_sarnb)),
  `ZINB-RE` = c(fit_zinbre$loglik, length(fit_zinbre$coefficients)),
  `SAR-NB-RE` = c(fit_sarnbre$loglik, length(fit_sarnbre$coefficients)),
  ZISAR   = c(fit_zisar$loglik, length(fit_zisar$coefficients)))
info <- cbind(info, AIC = 2 * info[, 2] - 2 * info[, 1],
              BIC = log(n) * info[, 2] - 2 * info[, 1])
colnames(info)[1:2] <- c("logLik", "k")

sink("solar_results.txt")
cat("Domestic PV installed under the FiT, 8 Feb 2016 - 31 Mar 2019, per LSOA, England\n")
cat("n =", n, " installations =", sum(y), " zeros =", sum(y == 0),
    sprintf("(%.2f%%)", 100 * mean(y == 0)), " ones =", sum(y == 1),
    " mean =", round(mean(y), 2), " var =", round(var(y), 1),
    " median =", median(y), " max =", max(y), "\n")
cat("W queen: mean neighbours", round(mean(card(lw$neighbours)), 2),
    " islands relinked:", length(Ws$islands), "\n\n")
cat("--- Descriptives (raw covariates) ---\n")
print(round(t(sapply(solar[, c("pv_dom_late", "pv_dom", "households", "sh_detached",
  "sh_social", "income_dep", "sh_nssec12", "sh_age65", "sh_children", "sh_offgas", "sh_band_eh", "elec_kwh", "north",
  "sh_flat", "density", "sh_private", "sh_communal", "sh_empty", "sh_pre1919")], function(v) c(mean = mean(v), sd = sd(v),
  min = min(v), median = median(v), max = max(v)))), 3))
cat("\n--- Correlations of standardised covariates ---\n")
print(round(cor(cbind(cov_x, cov_z)), 2))
cat("\n--- NB ---\n"); print(summary(fit_nb)$coefficients); cat("theta =", fit_nb$theta, "SE", fit_nb$SE.theta, "\n")
cat("\n--- ZINB ---\n"); print(summary(fit_zinb)$coefficients); cat("theta =", fit_zinb$theta, "\n")
cat("\n--- SAR-NB ---\n"); print(summary(fit_sarnb)$coefficients); cat("theta =", fit_sarnb$theta, "\n")
cat("\n--- ZISAR (convergence", fit_zisar$convergence, ", elapsed", round(timing$zisar, 1), "s) ---\n")
print(fit_zisar$coef_table); print(fit_zisar$transformed)
cat("\n--- ZINB-RE (rho = 0) ---\n"); print(fit_zinbre$coef_table); print(unlist(fit_zinbre$params[c("kappa", "sigma")]))
cat("\n--- SAR-NB-RE (pi = 0) ---\n"); print(fit_sarnbre$coef_table); print(fit_sarnbre$transformed)
cat("\n--- Information criteria ---\n"); print(round(info, 1))
cat(sprintf("\nLR test rho = 0: %.2f (p = %.3g); Wald z on atanh(rho): %.2f\n",
            lr_rho, pchisq(lr_rho, 1, lower.tail = FALSE), wald_rho$z))
cat("\n--- Effects on E[y] per 1 SD (avg over LSOAs; 'zero' = inflation channel, in direct) ---\n")
for (e in eff[-1]) cat(sprintf("%-11s direct %7.4f (%.4f) indirect %7.4f (%.4f) total %7.4f (%.4f) zero %7.4f\n",
  e$name, e$avg_direct, e$se_direct, e$avg_indirect, e$se_indirect,
  e$avg_total, e$se_total, e$avg_zero))
cat("elapsed", round(timing$effects, 1), "s\n")
cat("\n--- Observed vs expected frequencies ---\n"); print(round(freq_tab, 1))
cat("\n--- Structural zeros ---\n")
cat("mean pi_hat:", round(mean(pi_hat), 4), " sum pi_hat:", round(sum(pi_hat), 1),
    " expected structural zeros among observed zeros:", round(sum(post_struct), 1),
    " zeros with P(structural) > 0.5:", sum(post_struct > 0.5), "\n")
cat("\n--- Moran's I (queen) ---\n")
for (m in names(moran)) cat(sprintf("%-6s I = %.4f  (p = %.3g)\n", m,
                                    moran[[m]]$estimate[1], moran[[m]]$p.value))
cat("\n--- Sensitivity to W ---\n")
st <- sapply(c(list(queen = fit_zisar), sens), function(f)
  c(rho = f$params$rho, f$params$beta[-1], kappa = f$params$kappa,
    sigma = f$params$sigma, f$params$gamma[-1], logLik = f$loglik, AIC = f$AIC))
rownames(st) <- c("rho", colnames(X)[-1], "kappa", "sigma", colnames(Z)[-1], "logLik", "AIC")
print(round(st, 3))
cat("\n--- Robustness: cumulative 2010-2019 count ---\n")
cat("zeros =", sum(y_cum == 0), " mean =", round(mean(y_cum), 2), " var =", round(var(y_cum), 1), "\n")
cat("AIC: NB", round(AIC(rob$nb), 1), " ZINB", round(AIC(rob$zinb), 1),
    " ZISAR", round(rob$zisar$AIC, 1), " (convergence", rob$zisar$convergence, ")\n")
cat("expected zeros: NB", round(sum(dnbinom(0, size = rob$nb$theta, mu = fitted(rob$nb))), 1),
    " ZINB", round(sum(predict(rob$zinb, type = "prob")[, 1]), 1),
    " ZISAR", round(zisar_marginal_freq(rob$zisar, 0), 1), "\n")
print(rob$zisar$coef_table); print(rob$zisar$transformed)
print(round(summary(rob$zinb)$coefficients$zero, 3))
cat("\nRestarts from the stopping point (the second with kappa fixed):\n")
for (r in fits$cum_restarts) cat(sprintf(
  "  convergence %d  logLik %.1f  rho %.3f  kappa %.1f  sigma %.3f  gamma %s\n",
  r$convergence, r$loglik, r$params$rho, r$params$kappa, r$params$sigma,
  paste(sprintf("%.3f", r$params$gamma), collapse = " ")))
sink()

# -----------------------------------------------------------------------------
# 8. Maps: adoption rate and structural-zero probability (England, London)
# -----------------------------------------------------------------------------
map <- readRDS("../data/solar/lsoa_map.rds")
idx <- match(map$unit, solar$unit)
map$rate <- 1000 * y[idx] / solar$households[idx]
map$pi   <- pi_hat[idx]
map$zero <- y[idx] == 0
london <- substr(solar$lad19cd[idx], 1, 3) == "E09"
theme_map <- theme_void(base_size = 8) +
  theme(legend.position = "bottom", legend.key.width = unit(0.6, "cm"),
        legend.key.height = unit(0.22, "cm"), legend.title.position = "top",
        legend.title = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5))
g1 <- ggplot(map) + geom_sf(aes(fill = pmin(rate, 10)), colour = NA) +
  scale_fill_viridis_c(option = "magma", name = "installations per 1,000 households",
                       breaks = c(0, 5, 10), labels = c("0", "5", "10+")) +
  ggtitle("(a) England") + theme_map
g2 <- ggplot(map[london, ]) + geom_sf(aes(fill = zero), colour = NA) +
  scale_fill_manual(values = c(`FALSE` = "grey85", `TRUE` = "#b2182b"),
                    labels = c("at least one", "none"), name = "installations") +
  ggtitle("(b) London: observed zeros") + theme_map
g3 <- ggplot(map[london, ]) + geom_sf(aes(fill = pi), colour = NA) +
  scale_fill_viridis_c(option = "viridis", trans = "sqrt", limits = c(0, 1),
                       breaks = c(0.01, 0.1, 0.3, 0.6, 1),
                       name = "structural-zero probability") +
  ggtitle(expression(paste("(c) London: ", hat(pi)[i]))) + theme_map
dir.create("../tex/figures", showWarnings = FALSE)
ggsave("../tex/figures/solar_maps.png", g1 + g2 + g3, width = 7.2, height = 3.6, dpi = 300)
