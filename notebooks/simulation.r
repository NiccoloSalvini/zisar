# =============================================================================
#  ZISAR Monte Carlo Simulation Study
# =============================================================================

source("zisar.r")

library(pscl)       # zeroinfl() for standard ZINB
library(MASS)       # glm.nb()
library(Matrix)
library(parallel)

# -----------------------------------------------------------------------------
# Simulation design
# -----------------------------------------------------------------------------
SCENARIOS <- expand.grid(
  n_side   = c(10, 20, 30),        # n = 100, 400, 900
  rho_true = c(0.0, 0.3, 0.6),
  zi_rate  = c("low", "high"),      # ~15% vs ~35% zeros
  stringsAsFactors = FALSE
)

BETA_TRUE  <- c(0.5, 0.8, -0.5)    # intercept + 2 covariates
KAPPA_TRUE <- 2                     # NB shape (higher = less overdispersion)
SIGMA_TRUE <- 0.5                   # RE std dev

gamma_low  <- c(-2.0, 0.3)         # ~15% structural zeros
gamma_high <- c(-0.5, 0.5)         # ~35% structural zeros

N_REP <- 200                        # Monte Carlo replications

# -----------------------------------------------------------------------------
# Single-replication worker
# -----------------------------------------------------------------------------
run_one_rep <- function(rep_id, sc) {

  gamma_true <- if (sc$zi_rate == "low") gamma_low else gamma_high

  dat <- rzisar(
    n_side = sc$n_side,
    beta   = BETA_TRUE,
    rho    = sc$rho_true,
    gamma  = gamma_true,
    kappa  = KAPPA_TRUE,
    sigma  = SIGMA_TRUE,
    seed   = rep_id * 1000 + sc$n_side * 100 +
             as.integer(sc$rho_true * 10) +
             ifelse(sc$zi_rate == "high", 50, 0)
  )

  n <- length(dat$y)
  X <- dat$X
  Z <- dat$Z
  W <- dat$W
  y <- dat$y

  results <- list()

  # ---- Model 1: Standard ZINB (no spatial) ----
  tryCatch({
    df <- data.frame(y = y, x1 = X[, 2], x2 = X[, 3], z1 = Z[, 2])
    t0 <- proc.time()
    fit_zinb <- zeroinfl(y ~ x1 + x2 | z1, data = df, dist = "negbin")
    time_zinb <- (proc.time() - t0)[3]

    results[["ZINB"]] <- data.frame(
      model     = "ZINB",
      beta0_hat = coef(fit_zinb, "count")[1],
      beta1_hat = coef(fit_zinb, "count")[2],
      beta2_hat = coef(fit_zinb, "count")[3],
      rho_hat   = NA,
      kappa_hat = fit_zinb$theta,
      loglik    = logLik(fit_zinb),
      time      = time_zinb,
      converged = fit_zinb$converged,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    results[["ZINB"]] <<- data.frame(
      model = "ZINB", beta0_hat = NA, beta1_hat = NA, beta2_hat = NA,
      rho_hat = NA, kappa_hat = NA, loglik = NA, time = NA, converged = FALSE,
      stringsAsFactors = FALSE
    )
  })

  # ---- Model 2: SAR-NB (spatial lag on eta via 2SLS plug-in, no ZI) ----
  tryCatch({
    Wy <- as.numeric(W %*% y)
    WX  <- as.matrix(W %*% X)
    W2X <- as.matrix(W %*% WX)
    H <- cbind(X, WX, W2X)
    Wy_hat <- as.numeric(H %*% solve(crossprod(H)) %*% crossprod(H, Wy))

    df_sarnb <- data.frame(y = y, x1 = X[, 2], x2 = X[, 3], Wy = Wy_hat)
    t0 <- proc.time()
    fit_sarnb <- glm.nb(y ~ x1 + x2 + Wy, data = df_sarnb)
    time_sarnb <- (proc.time() - t0)[3]

    results[["SAR_NB"]] <- data.frame(
      model     = "SAR_NB",
      beta0_hat = coef(fit_sarnb)[1],
      beta1_hat = coef(fit_sarnb)[2],
      beta2_hat = coef(fit_sarnb)[3],
      rho_hat   = coef(fit_sarnb)[4],
      kappa_hat = fit_sarnb$theta,
      loglik    = logLik(fit_sarnb),
      time      = time_sarnb,
      converged = fit_sarnb$converged,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    results[["SAR_NB"]] <<- data.frame(
      model = "SAR_NB", beta0_hat = NA, beta1_hat = NA, beta2_hat = NA,
      rho_hat = NA, kappa_hat = NA, loglik = NA, time = NA, converged = FALSE,
      stringsAsFactors = FALSE
    )
  })

  # ---- Model 3: ZINB + eigenvector spatial filtering ----
  tryCatch({
    eig_W <- eigen(as.matrix(W))
    ev_real <- Re(eig_W$vectors[, 1:min(10, n)])
    df_esf <- data.frame(y = y, x1 = X[, 2], x2 = X[, 3], z1 = Z[, 2])
    for (k in 1:ncol(ev_real)) df_esf[[paste0("ev", k)]] <- ev_real[, k]
    ev_formula <- paste(paste0("ev", 1:ncol(ev_real)), collapse = " + ")
    f_str <- paste0("y ~ x1 + x2 + ", ev_formula, " | z1")

    t0 <- proc.time()
    fit_esf <- zeroinfl(as.formula(f_str), data = df_esf, dist = "negbin")
    time_esf <- (proc.time() - t0)[3]

    results[["ZINB_ESF"]] <- data.frame(
      model     = "ZINB_ESF",
      beta0_hat = coef(fit_esf, "count")[1],
      beta1_hat = coef(fit_esf, "count")[2],
      beta2_hat = coef(fit_esf, "count")[3],
      rho_hat   = NA,
      kappa_hat = fit_esf$theta,
      loglik    = logLik(fit_esf),
      time      = time_esf,
      converged = fit_esf$converged,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    results[["ZINB_ESF"]] <<- data.frame(
      model = "ZINB_ESF", beta0_hat = NA, beta1_hat = NA, beta2_hat = NA,
      rho_hat = NA, kappa_hat = NA, loglik = NA, time = NA, converged = FALSE,
      stringsAsFactors = FALSE
    )
  })

  # ---- Model 4: ZISAR (proposed) ----
  tryCatch({
    t0 <- proc.time()
    fit_zisar <- zisar_fit(y, X, Z, W, silent = TRUE)
    time_zisar <- (proc.time() - t0)[3]

    results[["ZISAR"]] <- data.frame(
      model     = "ZISAR",
      beta0_hat = fit_zisar$params$beta[1],
      beta1_hat = fit_zisar$params$beta[2],
      beta2_hat = fit_zisar$params$beta[3],
      rho_hat   = fit_zisar$params$rho,
      kappa_hat = fit_zisar$params$kappa,
      loglik    = fit_zisar$loglik,
      time      = time_zisar,
      converged = (fit_zisar$convergence == 0),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    results[["ZISAR"]] <<- data.frame(
      model = "ZISAR", beta0_hat = NA, beta1_hat = NA, beta2_hat = NA,
      rho_hat = NA, kappa_hat = NA, loglik = NA, time = NA, converged = FALSE,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, results)
  out$rep <- rep_id
  out$n   <- sc$n_side^2
  out$rho_true <- sc$rho_true
  out$zi_rate  <- sc$zi_rate
  out
}

# -----------------------------------------------------------------------------
# Run full simulation
# -----------------------------------------------------------------------------
run_simulation <- function(n_cores = 1, save_path = "mc_results.rds") {

  all_results <- list()
  counter <- 0

  for (s in seq_len(nrow(SCENARIOS))) {
    sc <- SCENARIOS[s, ]
    cat(sprintf("\n=== Scenario %d/%d: n=%d, rho=%.1f, ZI=%s ===\n",
                s, nrow(SCENARIOS), sc$n_side^2, sc$rho_true, sc$zi_rate))

    sc_results <- lapply(seq_len(N_REP), function(r) {
      if (r %% 20 == 0) cat(sprintf("  rep %d/%d\n", r, N_REP))
      run_one_rep(r, sc)
    })

    counter <- counter + 1
    all_results[[counter]] <- do.call(rbind, sc_results)
  }

  results_df <- do.call(rbind, all_results)
  saveRDS(results_df, save_path)
  cat("\nResults saved to", save_path, "\n")
  results_df
}

# -----------------------------------------------------------------------------
# Summarise Monte Carlo results
# -----------------------------------------------------------------------------
summarise_mc <- function(results_df) {

  summary_list <- list()
  for (mdl in unique(results_df$model)) {
    for (nn in unique(results_df$n)) {
      for (rr in unique(results_df$rho_true)) {
        for (zz in unique(results_df$zi_rate)) {
          sub <- results_df[results_df$model == mdl &
                            results_df$n == nn &
                            results_df$rho_true == rr &
                            results_df$zi_rate == zz, ]
          sub <- sub[sub$converged == TRUE, ]
          if (nrow(sub) < 10) next

          bias_beta1 <- mean(sub$beta1_hat, na.rm = TRUE) - BETA_TRUE[2]
          rmse_beta1 <- sqrt(mean((sub$beta1_hat - BETA_TRUE[2])^2, na.rm = TRUE))

          bias_rho <- if (all(is.na(sub$rho_hat))) NA else
                      mean(sub$rho_hat, na.rm = TRUE) - rr
          rmse_rho <- if (all(is.na(sub$rho_hat))) NA else
                      sqrt(mean((sub$rho_hat - rr)^2, na.rm = TRUE))

          summary_list[[length(summary_list) + 1]] <- data.frame(
            model    = mdl,
            n        = nn,
            rho_true = rr,
            zi_rate  = zz,
            n_conv   = nrow(sub),
            bias_beta1 = bias_beta1,
            rmse_beta1 = rmse_beta1,
            bias_rho   = bias_rho,
            rmse_rho   = rmse_rho,
            avg_time   = mean(sub$time, na.rm = TRUE),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  do.call(rbind, summary_list)
}

# =============================================================================
# To run (uncomment):
# results <- run_simulation(n_cores = 1)
# summary_table <- summarise_mc(results)
# print(summary_table)
# =============================================================================
