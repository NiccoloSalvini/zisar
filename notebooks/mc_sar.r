# Appendix 'Spatial inflation predictor: simulation checks': SAR filter on the
# inflation predictor, logit(pi) = (I - rho0 W)^{-1} Z gamma, 24 replications
# for each of n = 900, 3600 and 10,000. Writes mc_sar.rds
suppressPackageStartupMessages({library(TMB); library(Matrix); library(parallel); library(spdep)})
source("zisar.r"); .zisar_compile()
tru <- c(b0 = 1, b1 = 0.5, b2 = -0.3, g0 = -1, g1 = 0.8, rho = 0.4, sigma = 0.5, rho0 = 0.6)
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
  zeta <- as.numeric(solve(I_n - tru["rho0"] * W, Z %*% tru[4:5]))
  s <- rbinom(n, 1, plogis(zeta))
  list(y = ifelse(s == 1, 0L, rnbinom(n, size = 2, mu = exp(eta))), X = X, Z = Z, s = s)
}
res <- list(); for (k in c(30, 60, 100)) {
  W <- lattice_W(k); R <- 24
  t0 <- Sys.time()
  out <- mclapply(1:R, function(r) {
    d <- sim(k, W, 5000 + r)
    f <- tryCatch(zisar_fit(d$y, d$X, d$Z, W, zero_spatial = "sar"), error = function(e) NULL)
    if (is.null(f)) return(rep(NA, 13))
    c(f$params$beta, f$params$gamma, f$params$rho, f$params$sigma, f$params$rho0,
      f$loglik, f$max_grad, f$convergence, mean(d$y == 0), mean(d$s))
  }, mc.cores = 11, mc.preschedule = FALSE)
  G <- do.call(rbind, out); colnames(G) <- c(names(tru), "loglik", "max_grad", "conv", "zero", "struct")
  res[[as.character(k^2)]] <- G
  cat("\nn =", k^2, " R =", R, " conv==0:", sum(G[, "conv"] == 0, na.rm = TRUE),
      " max_grad median/max:", signif(median(G[, "max_grad"], na.rm = TRUE), 2), signif(max(G[, "max_grad"], na.rm = TRUE), 2),
      " zeros", round(mean(G[, "zero"]), 3), " structural", round(mean(G[, "struct"]), 3),
      " time", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
  print(round(rbind(true = tru, mean = colMeans(G[, 1:8], na.rm = TRUE),
                    sd = apply(G[, 1:8], 2, sd, na.rm = TRUE)), 3))
}
saveRDS(res, "mc_sar.rds")
