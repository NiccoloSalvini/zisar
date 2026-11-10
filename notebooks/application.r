# =============================================================================
#  ZISAR Empirical Application: Emergency Department Visits
#  across Italian Municipalities
# =============================================================================

source("zisar.r")

library(sf)
library(spdep)
library(pscl)
library(MASS)
library(ggplot2)
library(viridis)

# =============================================================================
# DATA PREPARATION
# Note: This script uses a simulated realistic dataset mimicking the structure
# of Italian municipality-level ED data.  Replace the simulation block below
# with real data loading when the actual dataset is available.
# =============================================================================

set.seed(2024)
n_muni <- 300

coords <- data.frame(
  lon = runif(n_muni, 7.0, 9.5),
  lat = runif(n_muni, 44.0, 46.0)
)

# Covariates (standardised)
pop         <- rnorm(n_muni, 0, 1)
aging_index <- rnorm(n_muni, 0, 1)
dist_hosp   <- rnorm(n_muni, 0, 1)
income      <- rnorm(n_muni, 0, 1)
gp_density  <- rnorm(n_muni, 0, 1)

# Spatial weights: 6-nearest neighbours, row-standardised
coords_mat <- as.matrix(coords)
knn6 <- knearneigh(coords_mat, k = 6)
nb6  <- knn2nb(knn6)
W_listw <- nb2listw(nb6, style = "W")
W <- listw2mat(W_listw)
W_sparse <- Matrix(W, sparse = TRUE)

X <- cbind(1, pop, aging_index, dist_hosp, income)
Z <- cbind(1, dist_hosp, gp_density)

# True parameters
beta_true  <- c(1.2, 0.4, 0.3, 0.2, -0.15)
gamma_true <- c(-1.5, 0.6, -0.4)
rho_true   <- 0.35
kappa_true <- 2
sigma_true <- 0.4

# Generate from ZISAR-NB DGP (SAR on eta)
u <- rnorm(n_muni, 0, sigma_true)
I_n <- Diagonal(n_muni)
A <- I_n - rho_true * W_sparse
eta <- as.numeric(solve(A, X %*% beta_true + u))
mu <- exp(eta)

pi_i <- 1 / (1 + exp(-(Z %*% gamma_true)))
s_i  <- rbinom(n_muni, 1, pi_i)

y_count <- rnbinom(n_muni, size = kappa_true, mu = mu)
y <- ifelse(s_i == 1, 0L, y_count)

cat("Outcome summary:\n")
cat("  n =", n_muni, "\n")
cat("  zeros:", sum(y == 0), "(", round(100 * mean(y == 0), 1), "%)\n")
cat("  mean:", round(mean(y), 2), " | var:", round(var(y), 2), "\n")
cat("  max:", max(y), "\n\n")

# =============================================================================
# MODEL ESTIMATION
# =============================================================================

df <- data.frame(y = y, pop = pop, aging = aging_index,
                 dist = dist_hosp, income = income, gp = gp_density)

# ---- Model 1: Standard ZINB (no spatial) ----
cat("=== Fitting standard ZINB ===\n")
fit_zinb <- zeroinfl(y ~ pop + aging + dist + income | dist + gp,
                     data = df, dist = "negbin")
cat("ZINB log-lik:", logLik(fit_zinb), " | AIC:", AIC(fit_zinb), "\n\n")

# ---- Model 2: SAR-NB (2SLS plug-in, no ZI) ----
cat("=== Fitting SAR-NB ===\n")
WX  <- as.matrix(W_sparse %*% X)
W2X <- as.matrix(W_sparse %*% WX)
H   <- cbind(X, WX[, -1], W2X[, -1])  # W1 = 1: drop duplicated intercepts
Wy  <- as.numeric(W_sparse %*% y)
Wy_hat <- as.numeric(H %*% solve(crossprod(H)) %*% crossprod(H, Wy))

df$Wy <- Wy_hat
fit_sarnb <- glm.nb(y ~ pop + aging + dist + income + Wy, data = df)
cat("SAR-NB log-lik:", logLik(fit_sarnb), " | AIC:", AIC(fit_sarnb), "\n\n")

# ---- Model 3: ZISAR (proposed) ----
cat("=== Fitting ZISAR (TMB/Laplace) ===\n")
fit_zisar <- zisar_fit(y, X, Z, W_sparse, silent = FALSE)

cat("\n")
print(fit_zisar)
cat("\n")

# ---- Model comparison ----
cat("=== Model Comparison ===\n")
comparison <- data.frame(
  Model  = c("ZINB", "SAR-NB", "ZISAR"),
  LogLik = c(as.numeric(logLik(fit_zinb)),
             as.numeric(logLik(fit_sarnb)),
             fit_zisar$loglik),
  AIC    = c(AIC(fit_zinb), AIC(fit_sarnb), fit_zisar$AIC),
  BIC    = c(BIC(fit_zinb), BIC(fit_sarnb), fit_zisar$BIC)
)
print(comparison)
cat("\n")

# =============================================================================
# EFFECTS DECOMPOSITION
# =============================================================================

cat("=== Direct / Indirect / Total Effects ===\n")
effects <- zisar_effects(fit_zisar, B = 500)
print(effects)

# =============================================================================
# DIAGNOSTICS: Residual spatial autocorrelation
# =============================================================================

cat("=== Moran's I on residuals ===\n")

resid_zinb <- residuals(fit_zinb, type = "pearson")
moran_zinb <- moran.test(resid_zinb, W_listw)

resid_sarnb <- residuals(fit_sarnb, type = "pearson")
moran_sarnb <- moran.test(resid_sarnb, W_listw)

# Pearson residuals of the ZINB mixture, conditional on the estimated field:
# E[y] = (1 - pi) mu,  Var[y] = (1 - pi) mu (1 + mu / kappa + pi mu)
pred_zisar <- zisar_predict(fit_zisar, type = "response")
mu_zisar   <- zisar_predict(fit_zisar, type = "count")
pi_zisar   <- zisar_predict(fit_zisar, type = "zero")
var_zisar  <- (1 - pi_zisar) * mu_zisar *
              (1 + mu_zisar / fit_zisar$params$kappa + pi_zisar * mu_zisar)
resid_zisar <- (y - pred_zisar) / sqrt(var_zisar)
moran_zisar <- moran.test(resid_zisar, W_listw)

cat(sprintf("  ZINB:    Moran's I = %.4f (p = %.4f)\n",
            moran_zinb$estimate[1], moran_zinb$p.value))
cat(sprintf("  SAR-NB:  Moran's I = %.4f (p = %.4f)\n",
            moran_sarnb$estimate[1], moran_sarnb$p.value))
cat(sprintf("  ZISAR:   Moran's I = %.4f (p = %.4f)\n",
            moran_zisar$estimate[1], moran_zisar$p.value))

# =============================================================================
# SENSITIVITY TO W (k-nearest neighbours, k = 4, 6, 8) and results file
# =============================================================================

sens <- lapply(c(4, 6, 8), function(k) {
  lw <- nb2listw(knn2nb(knearneigh(coords_mat, k = k)), style = "W")
  Wk <- Matrix(listw2mat(lw), sparse = TRUE)
  f  <- if (k == 6) fit_zisar else zisar_fit(y, X, Z, Wk)
  c(k = k, rho = f$params$rho, beta_pop = f$params$beta[2],
    beta_aging = f$params$beta[3], kappa = f$params$kappa, AIC = f$AIC)
})
sens <- do.call(rbind, sens)

sink("application_results.txt")
cat("Synthetic Piedmont-like data, seed 2024; truth: beta =", beta_true,
    " gamma =", gamma_true, " rho =", rho_true, " kappa =", kappa_true,
    " sigma =", sigma_true, "\n")
cat("n =", n_muni, " zeros =", sum(y == 0), " mean =", round(mean(y), 2),
    " var =", round(var(y), 2), "\n\n")
cat("--- ZINB ---\n"); print(summary(fit_zinb)$coefficients); cat("theta =", fit_zinb$theta, "\n")
cat("--- SAR-NB ---\n"); print(summary(fit_sarnb)$coefficients); cat("theta =", fit_sarnb$theta, "\n")
cat("--- ZISAR ---\n"); print(fit_zisar$coef_table); print(fit_zisar$transformed)
cat("\n--- Comparison ---\n"); print(comparison)
cat("\n--- Effects (covariate 1 = intercept) ---\n")
for (e in effects) cat(sprintf("cov %d: direct %.3f (%.3f) indirect %.3f (%.3f) total %.3f (%.3f)\n",
  e$covariate, e$avg_direct, e$se_direct, e$avg_indirect, e$se_indirect, e$avg_total, e$se_total))
cat("\n--- Moran's I (6-NN) ---\n")
mt_y <- moran.test(y, W_listw)
cat(sprintf("y: I = %.4f (p = %.4f)\n", mt_y$estimate[1], mt_y$p.value))
cat(sprintf("ZINB: I = %.4f (p = %.4f)\nSAR-NB: I = %.4f (p = %.4f)\nZISAR: I = %.4f (p = %.4f)\n",
  moran_zinb$estimate[1], moran_zinb$p.value, moran_sarnb$estimate[1], moran_sarnb$p.value,
  moran_zisar$estimate[1], moran_zisar$p.value))
cat("\n--- Sensitivity to W ---\n"); print(round(sens, 3))
sink()

# =============================================================================
# VISUALIZATION
# =============================================================================

coords$y <- y
coords$pred_zisar <- pred_zisar
coords$resid_zisar <- resid_zisar

p1 <- ggplot(coords, aes(x = lon, y = lat, colour = y)) +
  geom_point(size = 2.5) +
  scale_colour_viridis(option = "C", name = "ED visits") +
  labs(title = "Observed ED visit counts", x = "Longitude", y = "Latitude") +
  theme_minimal()

p2 <- ggplot(coords, aes(x = lon, y = lat, colour = pred_zisar)) +
  geom_point(size = 2.5) +
  scale_colour_viridis(option = "C", name = "Predicted") +
  labs(title = "ZISAR predicted", x = "Longitude", y = "Latitude") +
  theme_minimal()

p3 <- ggplot(coords, aes(x = lon, y = lat, colour = resid_zisar)) +
  geom_point(size = 2.5) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, name = "Residual") +
  labs(title = "ZISAR residuals", x = "Longitude", y = "Latitude") +
  theme_minimal()

cat("\nPlots generated. Save with ggsave() if needed.\n")
