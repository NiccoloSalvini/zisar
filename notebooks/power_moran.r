# =============================================================================
#  Power study: does a residual Moran's I test detect spatial dependence in
#  zero-inflated counts, compared with likelihood-based tests on rho?
#
#  DGP: the illustrative application (application.r) -- random centroids in the
#  Piedmont box, 6-NN row-standardised W, iid standard-normal covariates,
#  beta = (1.2, 0.4, 0.3, 0.2, -0.15), gamma = (-1.5, 0.6, -0.4), kappa = 2,
#  sigma = 0.4 -- with rho in {0, 0.2, 0.35, 0.5, 0.7} and n in {300, 900}.
#
#  Tests at the 5% level, per replication:
#    moran_y     Moran's I on the raw counts
#    moran_zinb  Moran's I on the Pearson residuals of an aspatial ZINB
#    wald_rho    ZISAR Wald test on rho_raw = atanh(rho)
#    lr_rho      ZISAR likelihood-ratio test of rho = 0
#
#  Output: power_moran.rds (all replications), power_moran.txt (rejection rates)
# =============================================================================

source("zisar.r")
suppressPackageStartupMessages({ library(spdep); library(pscl); library(parallel) })
.zisar_compile()
TMB::openmp(1)

BETA <- c(1.2, 0.4, 0.3, 0.2, -0.15); GAMMA <- c(-1.5, 0.6, -0.4)
KAPPA <- 2; SIGMA <- 0.4
RHOS <- c(0, 0.2, 0.35, 0.5, 0.7); NS <- c(300, 900); R <- 200
ALPHA <- 0.05

one_rep <- function(rep_id, n, rho) {
  set.seed(100000 * n + 1000 * round(100 * rho) + rep_id)
  coords <- cbind(runif(n, 7.0, 9.5), runif(n, 44.0, 46.0))
  lw <- nb2listw(knn2nb(knearneigh(coords, k = 6)), style = "W")
  W  <- Matrix(listw2mat(lw), sparse = TRUE)
  X  <- cbind(1, matrix(rnorm(n * 4), n))
  Z  <- cbind(1, X[, 4], rnorm(n))                  # dist_hosp enters both processes
  u   <- rnorm(n, 0, SIGMA)
  eta <- as.numeric(solve(Diagonal(n) - rho * W, X %*% BETA + u))
  pi_i <- plogis(as.numeric(Z %*% GAMMA))
  y <- ifelse(rbinom(n, 1, pi_i) == 1, 0L, rnbinom(n, size = KAPPA, mu = exp(eta)))

  out <- c(n = n, rho = rho, rep = rep_id, zeros = mean(y == 0),
           moran_y = NA, moran_zinb = NA, wald_rho = NA, lr_rho = NA,
           rho_hat = NA, I_y = NA, I_zinb = NA)

  out["I_y"] <- (mt <- moran.test(y, lw))$estimate[1]
  out["moran_y"] <- mt$p.value

  df <- data.frame(y = y, x1 = X[, 2], x2 = X[, 3], x3 = X[, 4], x4 = X[, 5], z2 = Z[, 3])
  zinb <- tryCatch(suppressWarnings(zeroinfl(y ~ x1 + x2 + x3 + x4 | x3 + z2, data = df, dist = "negbin")),
                   error = function(e) NULL)
  if (!is.null(zinb)) {
    mt <- moran.test(residuals(zinb, type = "pearson"), lw)
    out["I_zinb"] <- mt$estimate[1]; out["moran_zinb"] <- mt$p.value
  }

  fit <- tryCatch(zisar_fit(y, X, Z, W), error = function(e) NULL)
  if (!is.null(fit) && fit$convergence == 0) {
    out["rho_hat"] <- fit$params$rho
    ct <- fit$coef_table
    if (is.finite(ct["rho_raw", "p_value"])) out["wald_rho"] <- ct["rho_raw", "p_value"]
    # restricted fit with rho fixed at 0, same template and starting values
    obj0 <- fit$tmb_obj
    pars <- obj0$env$parList(fit$opt$par)
    pars$rho_raw <- 0
    obj_r <- MakeADFun(obj0$env$data, pars, random = "eta", DLL = "zisar_tmb",
                       silent = TRUE, map = list(rho_raw = factor(NA)))
    k <- length(obj_r$par)
    opt_r <- tryCatch(nlminb(obj_r$par, obj_r$fn, obj_r$gr,
                             lower = c(rep(-Inf, k - 2), log(0.01), log(0.001)),
                             upper = c(rep(Inf, k - 2), log(500), log(10)),
                             control = list(eval.max = 1000, iter.max = 500, rel.tol = 1e-10)),
                      error = function(e) NULL)
    if (!is.null(opt_r)) {
      lr <- max(0, 2 * (fit$loglik - (-opt_r$objective)))
      out["lr_rho"] <- pchisq(lr, 1, lower.tail = FALSE)
    }
  }
  out
}

grid <- expand.grid(rep = seq_len(R), rho = RHOS, n = NS)
if (nzchar(Sys.getenv("POWER_SMOKE"))) grid <- grid[grid$rep <= 2, ]
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(grid)), function(i) tryCatch(one_rep(grid$rep[i], grid$n[i], grid$rho[i]),
                                                         error = function(e) NULL),
                mc.cores = max(1, detectCores() - 4), mc.preschedule = FALSE)
res <- as.data.frame(do.call(rbind, res))
saveRDS(res, "power_moran.rds")

rate <- function(p) mean(p < ALPHA, na.rm = TRUE)
tab <- do.call(rbind, lapply(split(res, list(res$n, res$rho)), function(d)
  data.frame(n = d$n[1], rho = d$rho[1], reps = nrow(d), zeros = round(mean(d$zeros), 2),
             mean_rho_hat = round(mean(d$rho_hat, na.rm = TRUE), 3),
             moran_y = rate(d$moran_y), moran_zinb = rate(d$moran_zinb),
             wald_rho = rate(d$wald_rho), lr_rho = rate(d$lr_rho),
             mean_I_y = round(mean(d$I_y), 3), mean_I_zinb = round(mean(d$I_zinb, na.rm = TRUE), 3),
             zisar_ok = sum(is.finite(d$wald_rho)))))
tab <- tab[order(tab$n, tab$rho), ]
sink("power_moran.txt")
cat("Rejection rates at 5% (", R, "replications per cell); elapsed",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n\n")
print(tab, row.names = FALSE)
sink()
cat(readLines("power_moran.txt"), sep = "\n")
