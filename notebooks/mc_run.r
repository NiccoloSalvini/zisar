source("zisar.r")
library(pscl)
library(MASS)

BETA_TRUE  <- c(0.5, 0.8, -0.5)
GAMMA_TRUE <- c(-1.0, 0.4)
KAPPA_TRUE <- 2
SIGMA_TRUE <- 0.5
N_REP      <- 100

SCENARIOS <- expand.grid(
  n_side   = c(10, 15, 20),
  rho_true = c(0.0, 0.3, 0.6),
  stringsAsFactors = FALSE
)

run_one <- function(rep_id, n_side, rho_true) {
  seed <- rep_id * 1000 + n_side * 10 + round(rho_true * 10)
  dat <- rzisar(n_side, BETA_TRUE, rho_true, GAMMA_TRUE,
                KAPPA_TRUE, SIGMA_TRUE, seed = seed)
  y <- dat$y; X <- dat$X; Z <- dat$Z; W <- dat$W; n <- length(y)
  df <- data.frame(y=y, x1=X[,2], x2=X[,3], z1=Z[,2])

  out <- data.frame(rep=rep_id, n=n, rho_true=rho_true,
                    stringsAsFactors=FALSE)

  # --- ZINB ---
  r_zinb <- tryCatch({
    fit <- zeroinfl(y ~ x1 + x2 | z1, data=df, dist="negbin")
    data.frame(zinb_b0=coef(fit,"count")[1], zinb_b1=coef(fit,"count")[2],
               zinb_b2=coef(fit,"count")[3], zinb_kappa=fit$theta,
               zinb_ll=as.numeric(logLik(fit)), zinb_conv=fit$converged)
  }, error=function(e) data.frame(zinb_b0=NA,zinb_b1=NA,zinb_b2=NA,
                                   zinb_kappa=NA,zinb_ll=NA,zinb_conv=FALSE))

  # --- SAR-NB (2SLS) ---
  r_sarnb <- tryCatch({
    Wy <- as.numeric(W %*% y)
    WX <- as.matrix(W %*% X); W2X <- as.matrix(W %*% WX)
    H <- cbind(X, WX[,-1], W2X[,-1])
    Wy_hat <- as.numeric(H %*% qr.solve(crossprod(H), crossprod(H, Wy)))
    df$Wy <- Wy_hat
    fit <- suppressWarnings(glm.nb(y ~ x1 + x2 + Wy, data=df))
    data.frame(sarnb_b0=coef(fit)[1], sarnb_b1=coef(fit)[2],
               sarnb_b2=coef(fit)[3], sarnb_rho=coef(fit)[4],
               sarnb_kappa=fit$theta, sarnb_ll=as.numeric(logLik(fit)),
               sarnb_conv=fit$converged)
  }, error=function(e) data.frame(sarnb_b0=NA,sarnb_b1=NA,sarnb_b2=NA,
                                   sarnb_rho=NA,sarnb_kappa=NA,sarnb_ll=NA,
                                   sarnb_conv=FALSE))

  # --- ZISAR ---
  r_zisar <- tryCatch({
    t0 <- proc.time()
    fit <- zisar_fit(y, X, Z, W, silent=TRUE)
    tm <- (proc.time()-t0)[3]
    data.frame(zisar_b0=fit$params$beta[1], zisar_b1=fit$params$beta[2],
               zisar_b2=fit$params$beta[3], zisar_rho=fit$params$rho,
               zisar_kappa=fit$params$kappa, zisar_sigma=fit$params$sigma,
               zisar_ll=fit$loglik, zisar_time=tm,
               zisar_conv=(fit$convergence==0))
  }, error=function(e) data.frame(zisar_b0=NA,zisar_b1=NA,zisar_b2=NA,
                                   zisar_rho=NA,zisar_kappa=NA,zisar_sigma=NA,
                                   zisar_ll=NA,zisar_time=NA,zisar_conv=FALSE))

  cbind(out, r_zinb, r_sarnb, r_zisar)
}

cat("============================================================\n")
cat("  ZISAR Monte Carlo Simulation: R=100, 9 scenarios\n")
cat("============================================================\n\n")

all_res <- list()
idx <- 0
t_start <- proc.time()

for (s in seq_len(nrow(SCENARIOS))) {
  sc <- SCENARIOS[s,]
  n <- sc$n_side^2
  cat(sprintf("Scenario %d/9: n=%d, rho=%.1f  ", s, n, sc$rho_true))

  sc_res <- vector("list", N_REP)
  for (r in seq_len(N_REP)) {
    sc_res[[r]] <- run_one(r, sc$n_side, sc$rho_true)
    if (r %% 25 == 0) cat(".")
  }
  idx <- idx + 1
  all_res[[idx]] <- do.call(rbind, sc_res)
  cat(sprintf(" done (%.0fs)\n", (proc.time()-t_start)[3]))
}

results <- do.call(rbind, all_res)
saveRDS(results, "mc_results_100.rds")

cat(sprintf("\nTotal time: %.1f min\n\n", (proc.time()-t_start)[3]/60))

# ---- Summary ----
cat("============================================================\n")
cat("  RESULTS SUMMARY\n")
cat("============================================================\n\n")

for (nn in c(100, 225, 400)) {
  for (rr in c(0.0, 0.3, 0.6)) {
    sub <- results[results$n == nn & results$rho_true == rr, ]
    zconv <- sum(sub$zisar_conv, na.rm=TRUE)

    sub_z <- sub[sub$zisar_conv == TRUE & !is.na(sub$zisar_b0), ]
    sub_nb <- sub[sub$sarnb_conv == TRUE & !is.na(sub$sarnb_b0), ]
    sub_zi <- sub[sub$zinb_conv == TRUE & !is.na(sub$zinb_b0), ]

    cat(sprintf("--- n=%d, rho_true=%.1f (ZISAR converged: %d/%d) ---\n",
        nn, rr, zconv, N_REP))

    if (nrow(sub_z) < 5) { cat("  Too few converged.\n\n"); next }

    # Filter kappa outliers (boundary hits)
    sub_z_clean <- sub_z[sub_z$zisar_kappa < 100, ]
    n_clean <- nrow(sub_z_clean)

    cat(sprintf("%-10s %10s %10s %10s\n", "", "ZINB", "SAR-NB", "ZISAR"))
    cat(paste(rep("-", 44), collapse=""), "\n")

    # Bias beta1
    cat(sprintf("%-10s %10.3f %10.3f %10.3f\n", "Bias(b1)",
        mean(sub_zi$zinb_b1, na.rm=T) - BETA_TRUE[2],
        mean(sub_nb$sarnb_b1, na.rm=T) - BETA_TRUE[2],
        mean(sub_z_clean$zisar_b1, na.rm=T) - BETA_TRUE[2]))

    # RMSE beta1
    cat(sprintf("%-10s %10.3f %10.3f %10.3f\n", "RMSE(b1)",
        sqrt(mean((sub_zi$zinb_b1 - BETA_TRUE[2])^2, na.rm=T)),
        sqrt(mean((sub_nb$sarnb_b1 - BETA_TRUE[2])^2, na.rm=T)),
        sqrt(mean((sub_z_clean$zisar_b1 - BETA_TRUE[2])^2, na.rm=T))))

    # Bias beta2
    cat(sprintf("%-10s %10.3f %10.3f %10.3f\n", "Bias(b2)",
        mean(sub_zi$zinb_b2, na.rm=T) - BETA_TRUE[3],
        mean(sub_nb$sarnb_b2, na.rm=T) - BETA_TRUE[3],
        mean(sub_z_clean$zisar_b2, na.rm=T) - BETA_TRUE[3]))

    # RMSE beta2
    cat(sprintf("%-10s %10.3f %10.3f %10.3f\n", "RMSE(b2)",
        sqrt(mean((sub_zi$zinb_b2 - BETA_TRUE[3])^2, na.rm=T)),
        sqrt(mean((sub_nb$sarnb_b2 - BETA_TRUE[3])^2, na.rm=T)),
        sqrt(mean((sub_z_clean$zisar_b2 - BETA_TRUE[3])^2, na.rm=T))))

    # Bias/RMSE rho (SAR-NB and ZISAR only)
    cat(sprintf("%-10s %10s %10.3f %10.3f\n", "Bias(rho)", "--",
        mean(sub_nb$sarnb_rho, na.rm=T) - rr,
        mean(sub_z_clean$zisar_rho, na.rm=T) - rr))
    cat(sprintf("%-10s %10s %10.3f %10.3f\n", "RMSE(rho)", "--",
        sqrt(mean((sub_nb$sarnb_rho - rr)^2, na.rm=T)),
        sqrt(mean((sub_z_clean$zisar_rho - rr)^2, na.rm=T))))

    # Bias/RMSE kappa
    cat(sprintf("%-10s %10.3f %10.3f %10.3f\n", "Bias(k)",
        mean(sub_zi$zinb_kappa, na.rm=T) - KAPPA_TRUE,
        mean(sub_nb$sarnb_kappa, na.rm=T) - KAPPA_TRUE,
        mean(sub_z_clean$zisar_kappa, na.rm=T) - KAPPA_TRUE))

    # Median log-likelihood
    cat(sprintf("%-10s %10.1f %10.1f %10.1f\n", "Med(ll)",
        median(sub_zi$zinb_ll, na.rm=T),
        median(sub_nb$sarnb_ll, na.rm=T),
        median(sub_z_clean$zisar_ll, na.rm=T)))

    # Mean time
    cat(sprintf("%-10s %10s %10s %10.2f\n", "Time(s)", "<0.1", "<0.1",
        mean(sub_z_clean$zisar_time, na.rm=T)))

    cat(sprintf("  [%d/%d clean ZISAR fits (kappa<100)]\n\n", n_clean, nrow(sub_z)))
  }
}
