source("zisar.r")
library(pscl)
library(MASS)

cat("============================================\n")
cat("  ZISAR Quick Test: Parameter Recovery\n")
cat("============================================\n\n")

beta_true  <- c(0.5, 0.8, -0.5)
gamma_true <- c(-1.0, 0.4)
rho_true   <- 0.4
kappa_true <- 2
sigma_true <- 0.5

cat("TRUE PARAMETERS:\n")
cat("  beta  =", beta_true, "\n")
cat("  gamma =", gamma_true, "\n")
cat("  rho   =", rho_true, "\n")
cat("  kappa =", kappa_true, "\n")
cat("  sigma =", sigma_true, "\n\n")

run_test <- function(n_side, seed, label) {
  cat("============================================\n")
  cat(sprintf("  %s: n=%d (%dx%d lattice)\n", label, n_side^2, n_side, n_side))
  cat("============================================\n\n")

  dat <- rzisar(n_side = n_side, beta = beta_true, rho = rho_true,
                gamma = gamma_true, kappa = kappa_true, sigma = sigma_true,
                seed = seed)
  y <- dat$y; X <- dat$X; Z <- dat$Z; W <- dat$W; n <- length(y)

  cat("  zeros:", sum(y == 0), "/", n,
      "(", round(100*mean(y==0),1), "%)\n")
  cat("  mean(y) =", round(mean(y), 2),
      " var(y) =", round(var(y), 2), "\n\n")

  # --- ZINB ---
  cat("--- ZINB (no spatial) ---\n")
  df <- data.frame(y=y, x1=X[,2], x2=X[,3], z1=Z[,2])
  t0 <- proc.time()
  fit_zinb <- zeroinfl(y ~ x1 + x2 | z1, data=df, dist="negbin")
  cat("  time:", round((proc.time()-t0)[3], 2), "s\n")
  zinb_beta  <- coef(fit_zinb, "count")
  zinb_kappa <- fit_zinb$theta
  zinb_ll    <- as.numeric(logLik(fit_zinb))
  zinb_aic   <- AIC(fit_zinb)

  # --- SAR-NB (2SLS, no ZI) ---
  cat("--- SAR-NB (2SLS, no ZI) ---\n")
  Wy <- as.numeric(W %*% y)
  WX <- as.matrix(W %*% X); W2X <- as.matrix(W %*% WX)
  H <- cbind(X, WX[,-1], W2X[,-1])  # drop redundant intercept columns
  Wy_hat <- as.numeric(H %*% qr.solve(crossprod(H), crossprod(H, Wy)))
  df$Wy <- Wy_hat
  t0 <- proc.time()
  fit_sarnb <- suppressWarnings(glm.nb(y ~ x1 + x2 + Wy, data=df))
  cat("  time:", round((proc.time()-t0)[3], 2), "s\n")
  sarnb_beta  <- coef(fit_sarnb)[1:3]
  sarnb_rho   <- coef(fit_sarnb)["Wy"]
  sarnb_kappa <- fit_sarnb$theta
  sarnb_ll    <- as.numeric(logLik(fit_sarnb))
  sarnb_aic   <- AIC(fit_sarnb)

  # --- ZISAR ---
  cat("--- ZISAR (TMB/Laplace) ---\n")
  t0 <- proc.time()
  fit_zisar <- zisar_fit(y, X, Z, W, silent = TRUE)
  elapsed <- (proc.time()-t0)[3]
  cat("  time:", round(elapsed, 2), "s\n")
  cat("  convergence:", fit_zisar$convergence, "\n\n")

  # --- Table ---
  cat(sprintf("%-12s %8s %8s %8s %8s\n", "Parameter", "True", "ZINB", "SAR-NB", "ZISAR"))
  cat(paste(rep("-", 52), collapse=""), "\n")
  cat(sprintf("%-12s %8.3f %8.3f %8.3f %8.3f\n", "beta_0",
      beta_true[1], zinb_beta[1], sarnb_beta[1], fit_zisar$params$beta[1]))
  cat(sprintf("%-12s %8.3f %8.3f %8.3f %8.3f\n", "beta_1",
      beta_true[2], zinb_beta[2], sarnb_beta[2], fit_zisar$params$beta[2]))
  cat(sprintf("%-12s %8.3f %8.3f %8.3f %8.3f\n", "beta_2",
      beta_true[3], zinb_beta[3], sarnb_beta[3], fit_zisar$params$beta[3]))
  cat(sprintf("%-12s %8.3f %8s %8.3f %8.3f\n", "rho",
      rho_true, "--", sarnb_rho, fit_zisar$params$rho))
  cat(sprintf("%-12s %8.3f %8.3f %8.3f %8.3f\n", "kappa",
      kappa_true, zinb_kappa, sarnb_kappa, fit_zisar$params$kappa))
  cat(sprintf("%-12s %8.3f %8s %8s %8.3f\n", "sigma",
      sigma_true, "--", "--", fit_zisar$params$sigma))

  cat("\n")
  cat(sprintf("%-12s %8s %8.1f %8.1f %8.1f\n", "logLik", "",
      zinb_ll, sarnb_ll, fit_zisar$loglik))
  cat(sprintf("%-12s %8s %8.1f %8.1f %8.1f\n", "AIC", "",
      zinb_aic, sarnb_aic, fit_zisar$AIC))

  # --- Effects ---
  cat("\n--- ZISAR Direct/Indirect/Total Effects ---\n")
  eff <- zisar_effects(fit_zisar, B = 0)
  print(eff)
  cat("\n")

  invisible(list(zinb=fit_zinb, sarnb=fit_sarnb, zisar=fit_zisar))
}

res1 <- run_test(10, 42, "SMALL")
res2 <- run_test(15, 123, "MEDIUM")
res3 <- run_test(20, 7, "LARGE")

cat("Done!\n")
