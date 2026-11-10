# =============================================================================
#  Second empirical application: shoplifting in the 32,322 LSOA units of
#  England, July 2016 - June 2019, where the inflation field of the general
#  ZISAR is identifiable.
#
#  1. Restricted ZISAR with retail employment in both processes, and the
#     benchmarks of the Monte Carlo study: NB, ZINB, SAR-NB (instrumented
#     lag), ZINB-RE (rho = 0) and SAR-NB-RE (pi = 0)
#  2. Score-type test of a spatial field in the zero regime (England, Greater
#     London, outside London), parametric bootstrap
#  3. General ZISAR with the reduced-rank field (10 km knots), and the same
#     with police-force fixed effects in both processes
#  4. Checks on data simulated from the restricted fit (same X, Z, W): no
#     field; a field the basis can represent; a field constant within police
#     force areas; a unit-level SAR field. The unit-level field of the general
#     model on the no-field data shows the breakdown of the Laplace
#     approximation that the reduced-rank field avoids
#  5. Figure: shoplifting rate, predicted field, force means
#
#  Fits are cached in crime_fits.rds. Writes crime_results.rds,
#  crime_results.txt and ../tex/figures/crime_maps.png
# =============================================================================

source("zisar.r")
library(dplyr)
.zisar_compile()
set.seed(2026)

crime <- readRDS("../data/crime/shoplifting_lsoa.rds")
Ws    <- readRDS("../data/solar/W_lsoa.rds")
W     <- Ws$queen
n     <- nrow(crime)
y     <- crime$shoplifting
london <- crime$london

# Count process: retail and other employment (log of one plus), residents,
# density, income deprivation, social renting. Inflation process: an
# indicator of no recorded retail employment, retail employment (log of one
# plus) and density. (Retail bands in Z instead give the same picture, with
# the coefficients of the bands above 50 employees diverging to -Inf: no
# structural zero is attributed there.)
X <- cbind(intercept = 1,
           retail_jobs = log1p(crime$retail_jobs), other_jobs = log1p(crime$other_jobs),
           pop = log(crime$pop), density = log(crime$density),
           income_dep = crime$income_dep, social = crime$sh_social)
Z <- cbind(intercept = 1, no_retail = as.numeric(crime$retail_jobs == 0),
           retail_jobs = log1p(crime$retail_jobs), density = log(crime$density))
FE <- model.matrix(~ force, transform(crime, force = relevel(force, "Metropolitan Police Service")))[, -1]
colnames(FE) <- sub("^force", "", colnames(FE))
p <- ncol(X); q <- ncol(Z)

fits_file <- "crime_fits.rds"
fits <- if (file.exists(fits_file)) readRDS(fits_file) else list()
keep <- function(f) f[setdiff(names(f), "tmb_obj")]
start_from <- function(f, extra_x = 0, extra_z = 0) {
  pr <- f$params
  list(beta = c(pr$beta, rep(0, extra_x)), gamma = c(pr$gamma, rep(0, extra_z)),
       rho_raw = atanh(pr$rho), log_kappa = log(min(pr$kappa, 499)), log_sigma = log(pr$sigma))
}

# -----------------------------------------------------------------------------
# 1. Restricted ZISAR, and with police-force fixed effects
# -----------------------------------------------------------------------------
if (is.null(fits$zisar)) {
  t0 <- Sys.time()
  fits$zisar <- keep(zisar_fit(y, X, Z, W))
  fits$zisar$elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  saveRDS(fits, fits_file)
}
if (is.null(fits$force)) {
  t0 <- Sys.time()
  fits$force <- keep(zisar_fit(y, cbind(X, FE), cbind(Z, FE), W,
                               init = start_from(fits$zisar, ncol(FE), ncol(FE))))
  fits$force$elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  saveRDS(fits, fits_file)
}
fit <- fits$zisar; fit_f <- fits$force

# -----------------------------------------------------------------------------
# 1b. Benchmarks, as in the Monte Carlo study: NB, ZINB, SAR-NB (NB GLM with
#     the spatial lag W y instrumented by [X, WX, W^2X]), and the two nested
#     restrictions of ZISAR fitted with the template, ZINB-RE (rho = 0) and
#     SAR-NB-RE (pi = 0: inflation intercept fixed at -20, slopes at zero).
#     retail_jobs and density appear in both X and Z with the same values
# -----------------------------------------------------------------------------
df <- data.frame(y = y, X[, -1]); df$no_retail <- Z[, "no_retail"]
f_count <- paste(colnames(X)[-1], collapse = " + ")
f_zero  <- paste(colnames(Z)[-1], collapse = " + ")
if (is.null(fits$bench)) {
  fit_nb   <- MASS::glm.nb(as.formula(paste("y ~", f_count)), data = df)
  fit_zinb <- pscl::zeroinfl(as.formula(paste("y ~", f_count, "|", f_zero)), data = df, dist = "negbin")
  WX <- as.matrix(W %*% X[, -1]); H <- cbind(X, WX, as.matrix(W %*% WX))
  df$Wy_hat <- as.numeric(H %*% solve(crossprod(H), crossprod(H, as.numeric(W %*% y))))
  fit_sarnb <- MASS::glm.nb(as.formula(paste("y ~", f_count, "+ Wy_hat")), data = df)
  re <- parallel::mclapply(list(zinbre  = list(rho_raw = 0),
                                sarnbre = list(gamma = c(-20, rep(0, q - 1)))),
                           function(fx) keep(zisar_fit(y, X, Z, W, fixed = fx, init = start_from(fit))),
                           mc.cores = 2, mc.preschedule = FALSE)
  stopifnot(!sapply(re, inherits, "try-error"), sapply(re, function(f) is.finite(f$loglik)))
  fits$bench <- c(list(nb = fit_nb, zinb = fit_zinb, sarnb = fit_sarnb), re)
  saveRDS(fits, fits_file)
}
bench <- fits$bench

# -----------------------------------------------------------------------------
# 2. Score-type test of the zero regime, as in solar_application.r: r_i =
#    P(s_i = 1 | y_i) - pi_i under the restricted fit, Moran's I of r, null
#    by parametric bootstrap with the conditional mode of eta re-solved
# -----------------------------------------------------------------------------
rs <- function(M) Diagonal(x = 1 / pmax(rowSums(M), 1e-12)) %*% M
W_lon <- rs(W[london, london]); W_out <- rs(W[!london, !london])
moran_of <- function(x, Wm) { z <- x - mean(x); as.numeric(crossprod(z, Wm %*% z) / crossprod(z)) }
struct_test <- function(fit, B = 200, seed = 9000, cores = 4) {
  pr <- fit$params; Wc <- as(fit$W, "CsparseMatrix")
  pz <- plogis(as.numeric(fit$Z %*% pr$gamma)); Xb <- as.numeric(fit$X %*% pr$beta)
  data0 <- .zisar_data(fit$y, fit$X, fit$Z, fit$W, "sparse", 0L)
  eta_mode <- function(yy) {
    d <- data0; d$y <- yy
    par <- list(beta = pr$beta, gamma = pr$gamma, rho_raw = atanh(pr$rho), log_kappa = log(min(pr$kappa, 499)),
                log_sigma = log(pr$sigma), rho0_raw = 0, log_sigma0 = log(0.5), eta = Xb,
                g = numeric(0), alpha = numeric(0))
    map <- lapply(par[c("beta", "gamma", "rho_raw", "log_kappa", "log_sigma", "rho0_raw", "log_sigma0")],
                  function(v) factor(rep(NA, length(v))))
    ob <- MakeADFun(d, par, random = "eta", map = map, DLL = "zisar_tmb", silent = TRUE)
    ob$fn(ob$par); lp <- ob$env$last.par.best; as.numeric(lp[names(lp) == "eta"])
  }
  resid <- function(yy, eta) {
    p0 <- (pr$kappa / (pr$kappa + exp(eta)))^pr$kappa
    ifelse(yy == 0, pz / (pz + (1 - pz) * p0), 0) - pz
  }
  stat <- function(r) c(england = moran_of(r, Wc), london = moran_of(r[london], W_lon),
                        outside = moran_of(r[!london], W_out))
  obs <- stat(resid(fit$y, fit$eta_hat))
  sim <- do.call(rbind, parallel::mclapply(seq_len(B), function(b) {
    set.seed(seed + b)
    eta <- as.numeric(solve(Diagonal(n) - pr$rho * Wc, Xb + rnorm(n, 0, pr$sigma)))
    yy <- ifelse(rbinom(n, 1, pz) == 1, 0L, rnbinom(n, size = pr$kappa, mu = exp(eta)))
    stat(resid(yy, eta_mode(yy)))
  }, mc.cores = cores, mc.preschedule = FALSE))
  list(obs = obs, null_mean = colMeans(sim), null_sd = apply(sim, 2, sd),
       z = (obs - colMeans(sim)) / apply(sim, 2, sd),
       p = sapply(names(obs), function(k) (1 + sum(sim[, k] >= obs[k])) / (B + 1)), B = B)
}
if (is.null(fits$struct)) {
  fits$struct <- list(zisar = struct_test(fit), force = struct_test(fit_f))
  saveRDS(fits, fits_file)
}

# -----------------------------------------------------------------------------
# 3. General ZISAR with the reduced-rank inflation field: 10 km knot lattice,
#    bisquare functions of radius 15 km; profile in (rho0, sigma0) and fit
#    from the best grid point, without and with the force effects
# -----------------------------------------------------------------------------
basis <- zisar_basis(crime$east, crime$north, h = 10000)
if (is.null(fits$field)) {
  fits$field <- list(zisar = zisar_field_profile(fit, basis, cores = 4),
                     force = zisar_field_profile(fit_f, basis, cores = 4))
  saveRDS(fits, fits_file)
}

# -----------------------------------------------------------------------------
# 4. Simulated data with the real design. Fields in the inflation predictor:
#    none; smooth, g = B alpha with (rho0, sigma0) = (0.9, 0.5) on the knots;
#    step, constant within police force areas, N(0, 0.7^2); unit, the SAR field
#    of the general model with (0.9, 0.5). Each dataset gets a restricted fit
#    and the reduced-rank profile; the no-field data also get the unit-level
#    profile of the general model
# -----------------------------------------------------------------------------
simulate <- function(kind, seed) {
  set.seed(seed)
  pr <- fit$params; Wc <- as(W, "CsparseMatrix"); Wk <- as(basis$Wk, "CsparseMatrix")
  eta <- as.numeric(solve(Diagonal(n) - pr$rho * Wc, as.numeric(X %*% pr$beta) + rnorm(n, 0, pr$sigma)))
  g <- switch(kind,
    none   = rep(0, n),
    smooth = as.numeric(basis$B %*% solve(Diagonal(basis$K) - 0.9 * Wk, rnorm(basis$K, 0, 0.5))),
    step   = rnorm(nlevels(crime$force), 0, 0.7)[as.integer(crime$force)],
    unit   = as.numeric(solve(Diagonal(n) - 0.9 * Wc, rnorm(n, 0, 0.5))))
  s <- rbinom(n, 1, plogis(as.numeric(Z %*% pr$gamma) + g))
  list(y = ifelse(s == 1, 0L, rnbinom(n, size = pr$kappa, mu = exp(eta))), g = g, s = s)
}
if (is.null(fits$sim)) {
  kinds <- c(none = 101, smooth = 301, step = 401, unit = 201)
  fits$sim <- parallel::mclapply(names(kinds), function(kind) {
    d <- simulate(kind, kinds[kind])
    r <- zisar_fit(d$y, X, Z, W, init = start_from(fit))
    f <- zisar_field_profile(r, basis, cores = 1)
    out <- list(kind = kind, zeros = sum(d$y == 0), struct = sum(d$s), sd_g = sd(d$g),
                restricted = r[c("convergence", "max_grad", "loglik")], field = f,
                cor_g = if (kind == "none") NA else cor(f$g_hat, d$g))
    if (kind == "none") out$unit_field <- zisar_field_profile(r, cores = 1)
    out
  }, mc.cores = 4, mc.preschedule = FALSE)
  names(fits$sim) <- names(kinds)
  saveRDS(fits, fits_file)
}

# -----------------------------------------------------------------------------
# 5. Results
# -----------------------------------------------------------------------------
fld <- fits$field$zisar; fld_f <- fits$field$force
pz <- zisar_predict(fit, type = "zero"); mu <- zisar_predict(fit, type = "count")
p0 <- (fit$params$kappa / (fit$params$kappa + mu))^fit$params$kappa
post <- ifelse(y == 0, pz / (pz + (1 - pz) * p0), 0)
force_tab <- crime |> mutate(g = fld$g_hat, post = post) |> group_by(force) |>
  summarise(units = n(), zero = mean(shoplifting == 0), per_1000 = 1000 * sum(shoplifting) / sum(pop),
            g = mean(g), struct = sum(post), .groups = "drop") |> arrange(g)
r2_force <- summary(lm(fld$g_hat ~ crime$force))$r.squared
r2_force_f <- summary(lm(fld_f$g_hat ~ crime$force))$r.squared

# Fit of the six specifications. Expected zeros are marginal: for the
# random-effects models the latent field is integrated out by simulation from
# the fitted SAR distribution (200 draws)
marginal_zeros <- function(f, R = 200, seed = 7) {
  set.seed(seed)
  A <- Diagonal(f$n) - f$params$rho * as(f$W, "CsparseMatrix")
  Xb <- as.numeric(f$X %*% f$params$beta); pi_i <- plogis(as.numeric(f$Z %*% f$params$gamma)); kap <- f$params$kappa
  mean(sapply(seq_len(R), function(r) {
    eta <- as.numeric(solve(A, Xb + rnorm(f$n, 0, f$params$sigma)))
    sum(pi_i + (1 - pi_i) * (kap / (kap + exp(eta)))^kap)
  }))
}
ll <- function(f) as.numeric(logLik(f)); k_of <- function(f) attr(logLik(f), "df")
nb_zeros <- function(f) sum(dnbinom(0, size = f$theta, mu = fitted(f)))
pz_zinb <- predict(bench$zinb, type = "zero"); mu_zinb <- predict(bench$zinb, type = "count")
fit_tab <- rbind(
  NB          = c(ll(bench$nb), k_of(bench$nb), nb_zeros(bench$nb)),
  ZINB        = c(ll(bench$zinb), k_of(bench$zinb), sum(pz_zinb + (1 - pz_zinb) * dnbinom(0, size = bench$zinb$theta, mu = mu_zinb))),
  `SAR-NB`    = c(ll(bench$sarnb), k_of(bench$sarnb), nb_zeros(bench$sarnb)),
  `ZINB-RE`   = c(bench$zinbre$loglik, length(bench$zinbre$coefficients), marginal_zeros(bench$zinbre)),
  `SAR-NB-RE` = c(bench$sarnbre$loglik, length(bench$sarnbre$coefficients), marginal_zeros(bench$sarnbre)),
  ZISAR       = c(fit$loglik, length(fit$coefficients), marginal_zeros(fit)))
fit_tab <- cbind(logLik = fit_tab[, 1], k = fit_tab[, 2], AIC = 2 * fit_tab[, 2] - 2 * fit_tab[, 1],
                 BIC = log(n) * fit_tab[, 2] - 2 * fit_tab[, 1], exp_zeros = fit_tab[, 3])
# Residual dependence: Moran's I of the counts and of the Pearson residuals of
# the aspatial and instrumented models, and of the ZISAR innovations
lw <- spdep::nb2listw(Ws$nb_queen, style = "W")
mi <- function(x) { m <- spdep::moran.test(x, lw, randomisation = FALSE, zero.policy = TRUE)
                    c(I = unname(m$estimate[1]), z = unname(m$statistic)) }
u_hat <- as.numeric((Diagonal(n) - fit$params$rho * W) %*% fit$eta_hat - X %*% fit$params$beta)
moran_tab <- rbind(counts = mi(y), NB = mi(residuals(bench$nb, type = "pearson")),
                   ZINB = mi(residuals(bench$zinb, type = "pearson")), `SAR-NB` = mi(residuals(bench$sarnb, type = "pearson")),
                   `ZISAR u_hat` = mi(u_hat))

saveRDS(list(zisar = fit[c("params", "coef_table", "transformed", "loglik", "AIC", "BIC", "convergence", "max_grad", "elapsed")],
             force = fit_f[c("params", "coef_table", "transformed", "loglik", "AIC", "BIC", "convergence", "max_grad", "elapsed")],
             bench = list(fit_tab = fit_tab, moran_tab = moran_tab, nb = coef(summary(bench$nb)),
                          zinb = summary(bench$zinb)$coefficients, sarnb = coef(summary(bench$sarnb)),
                          zinbre = bench$zinbre[c("params", "coef_table", "loglik")], sarnbre = bench$sarnbre[c("params", "coef_table", "loglik")]),
             struct = fits$struct, field = lapply(fits$field, `[`, setdiff(names(fld), c("opt", "g_hat"))),
             g_hat = fld$g_hat, sim = fits$sim, force_tab = force_tab, r2_force = c(r2_force, r2_force_f),
             expected_struct = sum(pz), posterior_struct = sum(post)),
        "crime_results.rds")

sink("crime_results.txt")
cat("=== Shoplifting, England, Jul 2016 - Jun 2019: n =", n, " offences", sum(y), " zeros", sum(y == 0), "\n")
cat("\n--- Fit of the six specifications (expected zeros marginal; observed", sum(y == 0), ") ---\n")
print(round(fit_tab, 1))
cat(sprintf("LR rho = 0 (ZISAR vs ZINB-RE): %.2f; LR pi = 0 (ZISAR vs SAR-NB-RE): %.1f on %d restrictions\n",
            2 * (fit$loglik - bench$zinbre$loglik), 2 * (fit$loglik - bench$sarnbre$loglik), q))
cat(sprintf("expected structural zeros: ZINB %.0f, ZINB-RE %.0f, ZISAR %.0f\n", sum(pz_zinb),
            sum(zisar_predict(bench$zinbre, type = "zero")), sum(pz)))
cat("\n--- Moran's I (normal approximation) ---\n"); print(round(moran_tab, 3))
cat("\n--- NB ---\n"); print(round(coef(summary(bench$nb)), 3)); cat("theta =", bench$nb$theta, "\n")
cat("\n--- ZINB ---\n"); print(lapply(summary(bench$zinb)$coefficients, round, 3)); cat("theta =", bench$zinb$theta, "\n")
cat("\n--- SAR-NB (instrumented W y) ---\n"); print(round(coef(summary(bench$sarnb)), 3)); cat("theta =", bench$sarnb$theta, "\n")
cat("\n--- ZINB-RE (rho = 0): convergence", bench$zinbre$convergence, "---\n"); print(round(bench$zinbre$coef_table, 3)); print(round(bench$zinbre$transformed, 4))
cat("\n--- SAR-NB-RE (pi = 0): convergence", bench$sarnbre$convergence, "---\n"); print(round(bench$sarnbre$coef_table, 3)); print(round(bench$sarnbre$transformed, 4))
cat("\n--- Restricted ZISAR (convergence", fit$convergence, ", max gradient", signif(fit$max_grad, 3),
    ",", round(fit$elapsed, 1), "min) ---\n")
ct <- fit$coef_table; rownames(ct) <- c(colnames(X), paste0("zero:", colnames(Z)), "rho_raw", "log_kappa", "log_sigma")
print(round(ct, 3)); print(round(fit$transformed, 4))
cat(sprintf("loglik %.1f AIC %.1f BIC %.1f; expected structural zeros %.0f, posterior sum %.0f\n",
            fit$loglik, fit$AIC, fit$BIC, sum(pz), sum(post)))
cat("\n--- With police-force fixed effects (convergence", fit_f$convergence, ", max gradient", signif(fit_f$max_grad, 3), ") ---\n")
ct <- fit_f$coef_table; rownames(ct) <- c(colnames(X), colnames(FE), paste0("zero:", colnames(Z)), paste0("zero:", colnames(FE)),
                                          "rho_raw", "log_kappa", "log_sigma")
print(round(ct[c(seq_len(p), p + ncol(FE) + seq_len(q), nrow(ct) - 2:0), ], 3)); print(round(fit_f$transformed, 4))
cat(sprintf("loglik %.1f; LR force effects %.1f on %d df; expected structural zeros %.0f\n", fit_f$loglik,
            2 * (fit_f$loglik - fit$loglik), 2 * ncol(FE), sum(zisar_predict(fit_f, type = "zero"))))
cat("\n--- Score-type test of the zero regime (Moran's I of r, bootstrap B =", fits$struct$zisar$B, ") ---\n")
for (k in names(fits$struct)) { cat(k, "\n"); print(round(do.call(rbind, fits$struct[[k]][c("obs", "null_mean", "null_sd", "z", "p")]), 4)) }
cat("\n--- Reduced-rank field, 10 km knots (K =", basis$K, ") ---\n")
for (k in names(fits$field)) {
  x <- fits$field[[k]]
  cat(sprintf("%s: LR %.1f rho0 %.3f sigma0 %.3f rho %.3f kappa %.2f conv %d max_grad %.3g sd(g_hat) %.3f (%.1f min)\n",
              k, x$LR, x$rho0, x$sigma0, x$rho, x$kappa, x$convergence, x$max_grad, sd(x$g_hat), x$elapsed))
  print(xtabs(round(gain, 1) ~ rho0 + sigma0, x$grid))
}
cat(sprintf("R2 of g_hat on police force: %.3f without, %.3f with force effects\n", r2_force, r2_force_f))
cat("\nforce means of g_hat, zero share and rate per 1,000 residents:\n")
print(as.data.frame(force_tab |> mutate(across(c(zero, per_1000, g, struct), ~ round(.x, 3)))))
cat("\n--- Simulated data (real design): restricted fit, then reduced-rank profile ---\n")
for (k in names(fits$sim)) {
  s <- fits$sim[[k]]; f <- s$field
  cat(sprintf("%s: zeros %d structural %d sd(g) %.2f | restricted conv %d max_grad %.3g | field LR %.1f rho0 %.3f sigma0 %.3f conv %d max_grad %.3g cor(g_hat, g) %.3f\n",
              k, s$zeros, s$struct, s$sd_g, s$restricted$convergence, s$restricted$max_grad,
              f$LR, f$rho0, f$sigma0, f$convergence, f$max_grad, s$cor_g))
  print(xtabs(round(gain, 1) ~ rho0 + sigma0, f$grid))
  if (!is.null(s$unit_field)) {
    u <- s$unit_field
    cat(sprintf("  unit-level field on the same data: LR %.1f rho0 %.3f sigma0 %.3f conv %d (%s) max_grad %.3g\n",
                u$LR, u$rho0, u$sigma0, u$convergence, u$message, u$max_grad))
    print(xtabs(round(gain, 1) ~ rho0 + sigma0, u$grid))
  }
}
sink()
cat(readLines("crime_results.txt"), sep = "\n")

# -----------------------------------------------------------------------------
# 6. Figure: shoplifting per 1,000 residents, predicted field, force means
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({library(sf); library(ggplot2); library(patchwork)})
geo <- readRDS("../data/solar/lsoa_map.rds")
geo <- geo[match(crime$unit, geo$unit), ]
geo$rate <- 1000 * y / crime$pop
geo$g <- fld$g_hat
geo$g_force <- force_tab$g[match(crime$force, force_tab$force)]
forces <- geo |> group_by(force = crime$force) |> summarise(.groups = "drop")
theme_map <- theme_void(base_size = 8) + theme(legend.position = "bottom", legend.key.height = unit(2, "mm"),
                                              legend.key.width = unit(8, "mm"), plot.title = element_text(size = 8))
g1 <- ggplot(geo) + geom_sf(aes(fill = log10(1 + y)), colour = NA) +
  scale_fill_viridis_c(name = "offences", breaks = log10(1 + c(0, 1, 10, 100, 1000)), labels = c("0", "1", "10", "100", "1,000")) +
  labs(title = "(a) Shoplifting offences, Jul 2016 - Jun 2019") + theme_map
g2 <- ggplot(geo) + geom_sf(aes(fill = g), colour = NA) + geom_sf(data = forces, fill = NA, colour = "white", linewidth = 0.15) +
  scale_fill_gradient2(name = expression(hat(g)[i]), low = "#2166ac", mid = "grey95", high = "#b2182b") +
  labs(title = "(b) Predicted inflation field, 10 km knots") + theme_map
g3 <- ggplot(geo) + geom_sf(aes(fill = g_force), colour = NA) + geom_sf(data = forces, fill = NA, colour = "white", linewidth = 0.15) +
  scale_fill_gradient2(name = "force mean", low = "#2166ac", mid = "grey95", high = "#b2182b") +
  labs(title = "(c) Mean of the field by police force") + theme_map
ggsave("../tex/figures/crime_maps.png", g1 + g2 + g3, width = 7.2, height = 3.9, dpi = 300)
