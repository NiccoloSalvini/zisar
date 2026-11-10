# =============================================================================
#  Solar application: reduced-rank inflation field on the main specification.
#  Profiles of the Laplace log-likelihood in (rho0, sigma0) with the field on
#  10 km and 5 km knot lattices, and the same 10 km profile on one dataset
#  simulated from the restricted fit without a field (a check that the
#  estimator does not find a field where there is none).
#  Reads data/solar/solar_fits.rds (solar_application.r); writes
#  data/solar/solar_rr.rds
# =============================================================================
source("zisar.r")
.zisar_compile()
set.seed(2026)
fits  <- readRDS("../data/solar/solar_fits.rds")
fit   <- fits$zisar
solar <- readRDS("../data/solar/solar_lsoa.rds")
pr <- fit$params; n <- fit$n; Wc <- as(fit$W, "CsparseMatrix")

basis10 <- zisar_basis(solar$east, solar$north, h = 10000)
basis5  <- zisar_basis(solar$east, solar$north, h = 5000)
cat("knots: 10 km", basis10$K, " 5 km", basis5$K, "\n")

jobs <- list(
  rr10 = function() zisar_field_profile(fit, basis10, cores = 2),
  rr5  = function() zisar_field_profile(fit, basis5, cores = 2),
  null10 = function() {
    eta <- as.numeric(solve(Diagonal(n) - pr$rho * Wc, as.numeric(fit$X %*% pr$beta) + rnorm(n, 0, pr$sigma)))
    s <- rbinom(n, 1, plogis(as.numeric(fit$Z %*% pr$gamma)))
    y <- ifelse(s == 1, 0L, rnbinom(n, size = pr$kappa, mu = exp(eta)))
    r <- zisar_fit(y, fit$X, fit$Z, fit$W,
                   init = list(beta = pr$beta, gamma = pr$gamma, rho_raw = atanh(pr$rho),
                               log_kappa = log(min(pr$kappa, 499)), log_sigma = log(pr$sigma)))
    c(zisar_field_profile(r, basis10, cores = 2),
      list(restricted = r[c("convergence", "max_grad", "loglik")], zeros = sum(y == 0), struct = sum(s)))
  })
res <- parallel::mclapply(jobs, function(f) f(), mc.cores = 3, mc.preschedule = FALSE)
res$K <- c(k10 = basis10$K, k5 = basis5$K)
saveRDS(res, "../data/solar/solar_rr.rds")
for (k in c("rr10", "rr5", "null10")) {
  x <- res[[k]]
  if (inherits(x, "try-error")) { print(x); next }
  cat(sprintf("\n== %s: LR %.2f rho0 %.3f sigma0 %.3f rho %.3f kappa %.1f conv %d max_grad %.3g sd(g_hat) %.3f (%.1f min)\n",
              k, x$LR, x$rho0, x$sigma0, x$rho, x$kappa, x$convergence, x$max_grad, sd(x$g_hat), x$elapsed))
  print(xtabs(round(gain, 1) ~ rho0 + sigma0, x$grid))
}
