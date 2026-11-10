# Appendix 'Spatial inflation predictor: simulation checks': general ZISAR
# (SAR field g in the inflation predictor) estimated by Laplace approximation,
# 24 replications for each of n = 900, 3600 and 10,000. Writes mc_general.rds
suppressPackageStartupMessages({library(TMB); library(Matrix); library(parallel); library(spdep)})
source("zisar.r"); .zisar_compile()
tru <- c(b0 = 1, b1 = 0.5, b2 = -0.3, g0 = -1, g1 = 0.8, rho = 0.4, sigma = 0.5, rho0 = 0.6, sigma0 = 1)
lattice_W <- function(k) {
  nb <- cell2nb(k, k, type = "queen")
  i <- rep(seq_along(nb), lengths(nb)); j <- unlist(nb)
  W <- sparseMatrix(i, j, x = 1, dims = c(k^2, k^2))
  Diagonal(x = 1 / rowSums(W)) %*% W
}
sim <- function(k, W, seed) {
  set.seed(seed); n <- k^2; I_n <- Diagonal(n)
  X <- cbind(1, matrix(rnorm(n * 2), n, 2)); Z <- cbind(1, rnorm(n))
  eta <- as.numeric(solve(I_n - tru["rho"] * W, X %*% tru[1:3] + rnorm(n, 0, tru["sigma"])))
  g <- as.numeric(solve(I_n - tru["rho0"] * W, rnorm(n, 0, tru["sigma0"])))
  s <- rbinom(n, 1, plogis(as.numeric(Z %*% tru[4:5]) + g))
  list(y = ifelse(s == 1, 0L, rnbinom(n, size = 2, mu = exp(eta))), X = X, Z = Z)
}
res <- list(); for (k in c(30, 60, 100)) {
  W <- lattice_W(k); R <- 24
  out <- mclapply(1:R, function(r) {
    d <- sim(k, W, 3000 + r)
    f <- tryCatch(zisar_fit(d$y, d$X, d$Z, W, zero_spatial = TRUE), error = function(e) NULL)
    if (is.null(f)) return(rep(NA, 12))
    c(f$params$beta, f$params$gamma, f$params$rho, f$params$sigma, f$params$rho0, f$params$sigma0,
      f$loglik, f$max_grad, f$convergence)
  }, mc.cores = 12, mc.preschedule = FALSE)
  G <- do.call(rbind, out); colnames(G) <- c(names(tru), "loglik", "max_grad", "conv")
  res[[as.character(k^2)]] <- G; ok <- G[, "max_grad"] < 0.01 & G[, "loglik"] < 0
  cat("\nn =", k^2, " R =", R, " ok:", sum(ok, na.rm = TRUE), "\n")
  print(round(rbind(true = tru, mean = colMeans(G[ok, 1:9, drop = FALSE]),
                    sd = apply(G[ok, 1:9, drop = FALSE], 2, sd)), 3))
}
saveRDS(res, "mc_general.rds")
for (nm in names(res)) { G <- res[[nm]]
  cat("\nn =", nm, " conv==0:", sum(G[, "conv"] == 0, na.rm = TRUE), " max_grad quantiles:",
      signif(quantile(G[, "max_grad"], c(.25, .5, .75, 1), na.rm = TRUE), 2), " loglik>0:", sum(G[, "loglik"] > 0, na.rm = TRUE), "\n")
  print(round(rbind(mean_all = colMeans(G[, 1:9], na.rm = TRUE), sd_all = apply(G[, 1:9], 2, sd, na.rm = TRUE)), 3)) }
