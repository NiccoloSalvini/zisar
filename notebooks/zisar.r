# =============================================================================
#  ZISAR: Zero-Inflated Spatial Autoregressive Model
#  R interface to the TMB implementation
# =============================================================================

library(TMB)
library(Matrix)
library(MASS)

# Compile TMB template (only needed once per session)
.zisar_compile <- function(cpp_file = "zisar_tmb.cpp") {
  if (!file.exists(dynlib("zisar_tmb"))) {
    compile(cpp_file)
  }
  dyn.load(dynlib("zisar_tmb"))
}

# -----------------------------------------------------------------------------
# zisar_fit: ML estimation via Laplace approximation (TMB)
#
# Arguments:
#   y       - n-vector of counts
#   X       - n x p covariate matrix for count process (include intercept)
#   Z       - n x q covariate matrix for zero-inflation (include intercept)
#   W       - n x n spatial weights matrix (row-standardised, sparse or dense)
#   control - list of nlminb() control parameters
#   silent  - suppress TMB output
#
# Returns: S3 object of class "zisar"
# -----------------------------------------------------------------------------
zisar_fit <- function(y, X, Z, W,
                      control = list(eval.max = 1000, iter.max = 500,
                                     rel.tol = 1e-10),
                      silent = TRUE) {

  .zisar_compile()

  n <- length(y)
  p <- ncol(X)
  q <- ncol(Z)
  W_sp <- as(as(W, "CsparseMatrix"), "TsparseMatrix")

  # Pre-compute eigenvalues of W
  eig <- eigen(as.matrix(W), only.values = TRUE)$values
  eig_re <- Re(eig)
  eig_im <- Im(eig)

  # TMB data
  data <- list(
    y = as.numeric(y),
    X = as.matrix(X),
    Z = as.matrix(Z),
    W = W_sp,
    eigenvalues_W_re = eig_re,
    eigenvalues_W_im = eig_im
  )

  # Starting values from aspatial ZINB (or simple fallback)
  init_beta  <- tryCatch({
    fit0 <- glm(y ~ X[, -1], family = poisson())
    coef(fit0)
  }, error = function(e) rep(0, p))
  init_beta[is.na(init_beta)] <- 0

  init_gamma <- tryCatch({
    prop_zero <- mean(y == 0)
    c(qlogis(max(0.01, min(0.99, prop_zero))), rep(0, q - 1))
  }, error = function(e) rep(0, q))

  parameters <- list(
    beta      = init_beta,
    gamma     = init_gamma,
    rho_raw   = 0,
    log_kappa = log(1),
    log_sigma = log(0.5),
    eta       = as.numeric(X %*% init_beta)
  )

  # Build TMB objective
  obj <- MakeADFun(
    data       = data,
    parameters = parameters,
    random     = "eta",
    DLL        = "zisar_tmb",
    silent     = silent
  )

  # Parameter bounds: beta/gamma unbounded, rho_raw in [-4,4],
  # log_kappa in [log(0.01), log(500)], log_sigma in [log(0.001), log(10)]
  n_fixed <- p + q + 3
  lower <- c(rep(-Inf, p + q), -4, log(0.01), log(0.001))
  upper <- c(rep( Inf, p + q),  4, log(500),  log(10))

  # Optimise outer parameters
  opt <- nlminb(obj$par, obj$fn, obj$gr, lower = lower, upper = upper,
                control = control)

  # Retry with perturbed start if convergence failed
  if (opt$convergence != 0) {
    for (attempt in 1:3) {
      par_new <- opt$par + rnorm(length(opt$par), 0, 0.1)
      par_new <- pmin(pmax(par_new, lower), upper)
      opt2 <- nlminb(par_new, obj$fn, obj$gr, lower = lower, upper = upper,
                     control = control)
      if (opt2$objective < opt$objective) opt <- opt2
      if (opt$convergence == 0) break
    }
  }

  # Standard errors via sdreport
  sdr <- tryCatch(sdreport(obj), error = function(e) NULL)

  # Extract estimated parameters
  par_est <- opt$par
  idx <- 1
  beta_hat  <- par_est[idx:(idx + p - 1)]; idx <- idx + p
  gamma_hat <- par_est[idx:(idx + q - 1)]; idx <- idx + q
  rho_raw   <- par_est[idx]; idx <- idx + 1
  log_kappa <- par_est[idx]; idx <- idx + 1
  log_sigma <- par_est[idx]

  rho_hat   <- tanh(rho_raw)
  kappa_hat <- exp(log_kappa)
  sigma_hat <- exp(log_sigma)

  # Estimated latent field
  eta_hat <- tryCatch({
    re <- obj$env$last.par.best
    n_fixed <- length(opt$par)
    re[(n_fixed + 1):(n_fixed + n)]
  }, error = function(e) as.numeric(X %*% beta_hat))

  # Build coefficient table
  if (!is.null(sdr)) {
    summ_fixed <- summary(sdr, "fixed")
    summ_report <- summary(sdr, "report")

    par_names <- c(paste0("beta_", 0:(p-1)),
                   paste0("gamma_", 0:(q-1)),
                   "rho_raw", "log_kappa", "log_sigma")
    coef_table <- data.frame(
      estimate = summ_fixed[, 1],
      se       = summ_fixed[, 2],
      z        = summ_fixed[, 1] / summ_fixed[, 2],
      p_value  = 2 * pnorm(-abs(summ_fixed[, 1] / summ_fixed[, 2])),
      row.names = par_names
    )

    transformed_names <- rownames(summ_report)
    transformed_table <- data.frame(
      estimate = summ_report[, 1],
      se       = summ_report[, 2],
      row.names = transformed_names
    )
  } else {
    coef_table <- data.frame(
      estimate = par_est, se = NA, z = NA, p_value = NA
    )
    transformed_table <- NULL
  }

  # Variance-covariance of fixed effects
  vcov_mat <- tryCatch({
    H <- optimHess(opt$par, obj$fn, obj$gr)
    solve(H)
  }, error = function(e) matrix(NA, length(opt$par), length(opt$par)))

  structure(
    list(
      coefficients = par_est,
      params = list(
        beta  = beta_hat,
        gamma = gamma_hat,
        rho   = rho_hat,
        kappa = kappa_hat,
        sigma = sigma_hat
      ),
      eta_hat       = eta_hat,
      coef_table    = coef_table,
      transformed   = transformed_table,
      vcov          = vcov_mat,
      loglik        = -opt$objective,
      AIC           = 2 * length(opt$par) + 2 * opt$objective,
      BIC           = log(n) * length(opt$par) + 2 * opt$objective,
      convergence   = opt$convergence,
      n             = n,
      p             = p,
      q             = q,
      W             = W_sp,
      X             = X,
      Z             = Z,
      y             = y,
      eig_re        = eig_re,
      eig_im        = eig_im,
      tmb_obj       = obj,
      opt           = opt,
      sdr           = sdr
    ),
    class = "zisar"
  )
}

# -----------------------------------------------------------------------------
# print / summary methods
# -----------------------------------------------------------------------------
print.zisar <- function(x, ...) {
  cat("ZISAR: Zero-Inflated Spatial Autoregressive Model (TMB/Laplace)\n")
  cat("  n =", x$n, "\n")
  cat("  Log-likelihood:", round(x$loglik, 2), "\n")
  cat("  AIC:", round(x$AIC, 2), " | BIC:", round(x$BIC, 2), "\n\n")

  p <- x$params
  cat("  rho (spatial lag):", round(p$rho, 4), "\n")
  cat("  kappa (NB shape):", round(p$kappa, 4), "\n")
  cat("  sigma (RE std dev):", round(p$sigma, 4), "\n\n")

  cat("Coefficients:\n")
  print(round(x$coef_table, 4))
  if (!is.null(x$transformed)) {
    cat("\nTransformed parameters:\n")
    print(round(x$transformed, 4))
  }
  invisible(x)
}

summary.zisar <- function(object, ...) {
  cat("=== ZISAR Model Summary ===\n\n")
  cat("Count:     log(mu_i) = eta_i,  where eta = (I-rho*W)^{-1}(X*beta + u)\n")
  cat("Inflation: logit(pi_i) = Z'gamma\n")
  cat("Count:     y_i | s_i=0 ~ NB(exp(eta_i), kappa)\n\n")
  print(object)
  invisible(object)
}

# -----------------------------------------------------------------------------
# zisar_effects: Direct, indirect, total effects
#
# fit   - fitted zisar object
# B     - number of simulation draws for SEs (0 = point estimates only)
# -----------------------------------------------------------------------------
zisar_effects <- function(fit, B = 500) {
  n <- fit$n
  p <- fit$p
  W <- fit$W
  rho <- fit$params$rho
  beta <- fit$params$beta
  gamma <- fit$params$gamma
  eta_hat <- fit$eta_hat

  # Spatial multiplier S(rho) = (I - rho*W)^{-1}
  I_n <- Diagonal(n)
  A <- I_n - rho * W
  S_rho <- solve(A)

  # Fitted values at the mode
  mu_hat <- exp(eta_hat)
  pi_hat <- 1 / (1 + exp(-as.numeric(fit$Z %*% gamma)))

  # D = diag((1 - pi) * mu)
  d_vec <- (1 - pi_hat) * mu_hat

  effects <- list()
  for (r in seq_len(p)) {
    M_r <- Diagonal(x = d_vec) %*% S_rho * beta[r]
    M_dense <- as.matrix(M_r)

    direct_i   <- diag(M_dense)
    total_i    <- rowSums(M_dense)
    indirect_i <- total_i - direct_i

    effects[[r]] <- list(
      covariate    = r,
      avg_direct   = mean(direct_i),
      avg_indirect = mean(indirect_i),
      avg_total    = mean(total_i)
    )
  }

  # Simulation-based SEs
  if (B > 0 && all(is.finite(fit$vcov))) {
    n_par <- length(fit$coefficients)
    sim_effects <- array(0, dim = c(B, p, 3))

    theta_draws <- MASS::mvrnorm(B, fit$coefficients, fit$vcov)

    for (b in seq_len(B)) {
      tryCatch({
        idx <- 1
        beta_b  <- theta_draws[b, idx:(idx + p - 1)]; idx <- idx + p
        gamma_b <- theta_draws[b, idx:(idx + fit$q - 1)]; idx <- idx + fit$q
        rho_b   <- tanh(theta_draws[b, idx]); idx <- idx + 1
        # log_kappa and log_sigma not needed for effects

        if (abs(rho_b) >= 1) next

        S_b <- solve(I_n - rho_b * W)
        # Approximate eta with new beta and rho
        eta_b <- as.numeric(S_b %*% (fit$X %*% beta_b))
        mu_b  <- exp(eta_b)
        pi_b  <- 1 / (1 + exp(-as.numeric(fit$Z %*% gamma_b)))
        d_b   <- (1 - pi_b) * mu_b

        for (r in seq_len(p)) {
          M_b <- as.matrix(Diagonal(x = d_b) %*% S_b * beta_b[r])
          sim_effects[b, r, 1] <- mean(diag(M_b))
          sim_effects[b, r, 3] <- mean(rowSums(M_b))
          sim_effects[b, r, 2] <- sim_effects[b, r, 3] - sim_effects[b, r, 1]
        }
      }, error = function(e) NULL)
    }

    for (r in seq_len(p)) {
      valid <- sim_effects[, r, 1] != 0
      if (sum(valid) > 10) {
        effects[[r]]$se_direct   <- sd(sim_effects[valid, r, 1])
        effects[[r]]$se_indirect <- sd(sim_effects[valid, r, 2])
        effects[[r]]$se_total    <- sd(sim_effects[valid, r, 3])
        effects[[r]]$ci_direct   <- quantile(sim_effects[valid, r, 1], c(0.025, 0.975))
        effects[[r]]$ci_indirect <- quantile(sim_effects[valid, r, 2], c(0.025, 0.975))
        effects[[r]]$ci_total    <- quantile(sim_effects[valid, r, 3], c(0.025, 0.975))
      }
    }
  }

  structure(effects, class = "zisar_effects")
}

print.zisar_effects <- function(x, ...) {
  cat("ZISAR: Direct, Indirect, and Total Effects\n")
  cat(paste(rep("-", 65), collapse = ""), "\n")
  for (r in seq_along(x)) {
    eff <- x[[r]]
    cat(sprintf("Covariate %d:\n", eff$covariate))
    cat(sprintf("  Direct:   %8.4f", eff$avg_direct))
    if (!is.null(eff$se_direct)) cat(sprintf("  (SE: %.4f)", eff$se_direct))
    cat("\n")
    cat(sprintf("  Indirect: %8.4f", eff$avg_indirect))
    if (!is.null(eff$se_indirect)) cat(sprintf("  (SE: %.4f)", eff$se_indirect))
    cat("\n")
    cat(sprintf("  Total:    %8.4f", eff$avg_total))
    if (!is.null(eff$se_total)) cat(sprintf("  (SE: %.4f)", eff$se_total))
    cat("\n\n")
  }
  invisible(x)
}

# -----------------------------------------------------------------------------
# zisar_predict: Predicted values
# type = "response" (E[y]), "count" (mu), "zero" (pi), "eta" (latent field)
# -----------------------------------------------------------------------------
zisar_predict <- function(fit, newX = NULL, newZ = NULL,
                          type = c("response", "count", "zero", "eta")) {
  type <- match.arg(type)
  params <- fit$params

  eta <- fit$eta_hat
  mu  <- exp(eta)
  pi_hat <- 1 / (1 + exp(-as.numeric(fit$Z %*% params$gamma)))

  switch(type,
    response = (1 - pi_hat) * mu,
    count    = mu,
    zero     = pi_hat,
    eta      = eta
  )
}

# -----------------------------------------------------------------------------
# rzisar: Data-generating process for Monte Carlo
#
# n_side  - lattice side (n = n_side^2)
# beta    - count regression coefficients (including intercept)
# rho     - SAR coefficient on latent field
# gamma   - inflation regression coefficients (including intercept)
# kappa   - NB shape parameter
# sigma   - RE std dev
# seed    - random seed
# -----------------------------------------------------------------------------
rzisar <- function(n_side, beta, rho, gamma, kappa = 2, sigma = 0.5,
                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- n_side^2

  # Queen contiguity on regular lattice
  coords <- expand.grid(row = 1:n_side, col = 1:n_side)
  W_raw <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        dr <- abs(coords$row[i] - coords$row[j])
        dc <- abs(coords$col[i] - coords$col[j])
        if (dr <= 1 && dc <= 1) W_raw[i, j] <- 1
      }
    }
  }
  rs <- rowSums(W_raw)
  rs[rs == 0] <- 1
  W <- W_raw / rs
  W_sp <- Matrix(W, sparse = TRUE)

  p <- length(beta)
  q <- length(gamma)

  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  Z <- cbind(1, matrix(rnorm(n * (q - 1)), n, q - 1))

  # Generate latent field: eta = (I - rho*W)^{-1}(X*beta + u)
  u <- rnorm(n, 0, sigma)
  I_n <- Diagonal(n)
  A <- I_n - rho * W_sp
  eta <- as.numeric(solve(A, X %*% beta + u))
  mu <- exp(eta)

  # Zero-inflation
  pi_i <- 1 / (1 + exp(-as.numeric(Z %*% gamma)))
  s_i  <- rbinom(n, 1, pi_i)

  # NB counts
  y_count <- rnbinom(n, size = kappa, mu = mu)
  y <- ifelse(s_i == 1, 0L, y_count)

  list(
    y = as.integer(y),
    X = X,
    Z = Z,
    W = W_sp,
    coords = coords,
    params = list(beta = beta, rho = rho, gamma = gamma,
                  kappa = kappa, sigma = sigma),
    eta_true = eta,
    mu_true  = mu,
    pi_true  = as.numeric(pi_i),
    s_true   = s_i
  )
}
