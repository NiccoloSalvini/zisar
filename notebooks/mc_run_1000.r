# =============================================================================
#  Monte Carlo, R = 1000: same DGP and competitors as mc_run.r, on a finer rho
#  grid {0, 0.15, 0.3, 0.45, 0.6} (it contains the paper's {0, 0.3, 0.6}).
#  Output: mc_results_1000.rds (one row per replication)
#  Smoke run: MC_SMOKE=1 Rscript mc_run_1000.r
# =============================================================================

source("zisar.r")
suppressPackageStartupMessages({ library(pscl); library(MASS); library(parallel) })
.zisar_compile()
TMB::openmp(1)

BETA_TRUE  <- c(0.5, 0.8, -0.5)
GAMMA_TRUE <- c(-1.0, 0.4)
KAPPA_TRUE <- 2
SIGMA_TRUE <- 0.5
N_REP      <- if (nzchar(Sys.getenv("MC_SMOKE"))) 2 else 1000

run_one <- function(rep_id, n_side, rho_true) {
  seed <- rep_id * 1000 + n_side * 10 + round(rho_true * 100) %/% 5
  dat <- rzisar(n_side, BETA_TRUE, rho_true, GAMMA_TRUE,
                KAPPA_TRUE, SIGMA_TRUE, seed = seed)
  y <- dat$y; X <- dat$X; Z <- dat$Z; W <- dat$W; n <- length(y)
  df <- data.frame(y = y, x1 = X[, 2], x2 = X[, 3], z1 = Z[, 2])
  out <- data.frame(rep = rep_id, n = n, rho_true = rho_true, zeros = mean(y == 0))

  r_zinb <- tryCatch({
    fit <- suppressWarnings(zeroinfl(y ~ x1 + x2 | z1, data = df, dist = "negbin"))
    data.frame(zinb_b0 = coef(fit, "count")[1], zinb_b1 = coef(fit, "count")[2],
               zinb_b2 = coef(fit, "count")[3], zinb_kappa = fit$theta,
               zinb_ll = as.numeric(logLik(fit)), zinb_conv = fit$converged)
  }, error = function(e) data.frame(zinb_b0 = NA, zinb_b1 = NA, zinb_b2 = NA,
                                    zinb_kappa = NA, zinb_ll = NA, zinb_conv = FALSE))

  r_sarnb <- tryCatch({
    Wy <- as.numeric(W %*% y)
    WX <- as.matrix(W %*% X); W2X <- as.matrix(W %*% WX)
    H <- cbind(X, WX[, -1], W2X[, -1])
    df$Wy <- as.numeric(H %*% qr.solve(crossprod(H), crossprod(H, Wy)))
    fit <- suppressWarnings(glm.nb(y ~ x1 + x2 + Wy, data = df))
    data.frame(sarnb_b0 = coef(fit)[1], sarnb_b1 = coef(fit)[2],
               sarnb_b2 = coef(fit)[3], sarnb_rho = coef(fit)[4],
               sarnb_kappa = fit$theta, sarnb_ll = as.numeric(logLik(fit)),
               sarnb_conv = fit$converged)
  }, error = function(e) data.frame(sarnb_b0 = NA, sarnb_b1 = NA, sarnb_b2 = NA,
                                    sarnb_rho = NA, sarnb_kappa = NA, sarnb_ll = NA,
                                    sarnb_conv = FALSE))

  r_zisar <- tryCatch({
    t0 <- proc.time()
    fit <- zisar_fit(y, X, Z, W, silent = TRUE)
    tm <- (proc.time() - t0)[3]
    data.frame(zisar_b0 = fit$params$beta[1], zisar_b1 = fit$params$beta[2],
               zisar_b2 = fit$params$beta[3], zisar_rho = fit$params$rho,
               zisar_kappa = fit$params$kappa, zisar_sigma = fit$params$sigma,
               zisar_ll = fit$loglik, zisar_time = tm,
               zisar_conv = (fit$convergence == 0))
  }, error = function(e) data.frame(zisar_b0 = NA, zisar_b1 = NA, zisar_b2 = NA,
                                    zisar_rho = NA, zisar_kappa = NA, zisar_sigma = NA,
                                    zisar_ll = NA, zisar_time = NA, zisar_conv = FALSE))

  cbind(out, r_zinb, r_sarnb, r_zisar)
}

grid <- expand.grid(rep = seq_len(N_REP), rho_true = c(0, 0.15, 0.3, 0.45, 0.6),
                    n_side = c(10, 15, 20))
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(grid)), function(i)
  tryCatch(run_one(grid$rep[i], grid$n_side[i], grid$rho_true[i]), error = function(e) NULL),
  mc.cores = max(1, detectCores() - 4), mc.preschedule = TRUE)
res <- do.call(rbind, res)
rownames(res) <- NULL
saveRDS(res, if (N_REP == 2) "mc_results_smoke.rds" else "mc_results_1000.rds")
cat(nrow(res), "replications in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
