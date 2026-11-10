# =============================================================================
#  ZISAR: Zero-Inflated Spatial Autoregressive Model
#  R interface to the TMB implementation
# =============================================================================

library(TMB)
library(Matrix)
library(MASS)

# Compile TMB template (only needed once per session)
.zisar_compile <- function(cpp_file = "zisar_tmb.cpp") {
  if (!file.exists(dynlib("zisar_tmb")) ||
      file.mtime(cpp_file) > file.mtime(dynlib("zisar_tmb"))) {
    if (is.loaded("MakeADFunObject", PACKAGE = "zisar_tmb"))
      dyn.unload(dynlib("zisar_tmb"))
    # TMBad: needed for the AD sparse log-determinant (logdet = "sparse")
    compile(cpp_file, framework = "TMBad")
  }
  dyn.load(dynlib("zisar_tmb"))
}

# TMB data list; logdet is "eigen" or "sparse"
.zisar_data <- function(y, X, Z, W, logdet, zero_spatial, basis = NULL) {
  n <- length(y)
  W_sp <- as(as(W, "CsparseMatrix"), "TsparseMatrix")
  empty <- as(Matrix(0, 1, 1, sparse = TRUE), "TsparseMatrix")

  if (logdet == "eigen") {
    # Pre-compute eigenvalues of W
    eig <- eigen(as.matrix(W), only.values = TRUE)$values
    eig_re <- Re(eig)
    eig_im <- Im(eig)
    I_n <- WpWt <- WtW <- as(Matrix(0, 1, 1, sparse = TRUE), "TsparseMatrix")
  } else {
    eig_re <- eig_im <- numeric(0)
    Wc   <- as(W, "CsparseMatrix")
    I_n  <- as(Diagonal(n), "TsparseMatrix")
    WpWt <- as(Wc + t(Wc), "TsparseMatrix")
    WtW  <- as(crossprod(Wc), "TsparseMatrix")
  }

  list(
    y = as.numeric(y),
    X = as.matrix(X),
    Z = as.matrix(Z),
    W = W_sp,
    eigenvalues_W_re = eig_re,
    eigenvalues_W_im = eig_im,
    logdet_sparse = as.integer(logdet == "sparse"),
    I_n = I_n,
    WpWt = WpWt,
    WtW = WtW,
    zero_spatial = as.integer(zero_spatial),
    # reduced-rank inflation field (zero_spatial = 3): basis B (n x K), knot
    # weights W_K and their eigenvalues; dummies otherwise
    B = if (is.null(basis)) empty else basis$B,
    Wk = if (is.null(basis)) empty else basis$Wk,
    eig_k_re = if (is.null(basis)) numeric(0) else basis$eig_re,
    eig_k_im = if (is.null(basis)) numeric(0) else basis$eig_im
  )
}

# -----------------------------------------------------------------------------
# zisar_basis: basis of a reduced-rank inflation field, g = B alpha
#
# east, north - unit coordinates (same units as h)
# h           - spacing of the square knot lattice
# r_fac       - radius of the bisquare functions, as a multiple of h
#
# Knots lie on a square lattice of spacing h and are kept where a unit centroid
# is within h. Column k of B is the bisquare function (1 - d^2/r^2)^2 of the
# distance d from knot k, zero beyond r = r_fac * h, and the rows of B are
# normalised to sum to one, so that g_i is a local weighted mean of the knot
# values alpha. The knot weights W_K are rook contiguity on the lattice, row
# standardised, with isolated knots linked to the nearest knot; their
# eigenvalues give the log-determinant of (I_K - rho0 W_K).
# Returns a list for the basis argument of zisar_fit and zisar_objective.
# -----------------------------------------------------------------------------
zisar_basis <- function(east, north, h, r_fac = 1.5) {
  x0 <- floor(min(east) / h) * h
  y0 <- floor(min(north) / h) * h
  nx <- ceiling((max(east) - x0) / h) + 1
  ny <- ceiling((max(north) - y0) / h) + 1
  kn <- expand.grid(i = seq_len(nx) - 1, j = seq_len(ny) - 1)
  kn$x <- x0 + kn$i * h
  kn$y <- y0 + kn$j * h
  # units indexed by their nearest lattice cell
  cells <- split(seq_along(east), paste(round((east - x0) / h), round((north - y0) / h)))
  units_near <- function(i, j, m) {
    keys <- as.vector(outer(i + (-m:m), j + (-m:m), paste))
    unlist(cells[keys], use.names = FALSE)
  }
  keep <- vapply(seq_len(nrow(kn)), function(k) {
    u <- units_near(kn$i[k], kn$j[k], 1)
    length(u) > 0 && min((east[u] - kn$x[k])^2 + (north[u] - kn$y[k])^2) <= h^2
  }, logical(1))
  kn <- kn[keep, ]
  K <- nrow(kn)
  r <- r_fac * h
  m_r <- ceiling(r_fac + 0.5)
  trip <- do.call(rbind, lapply(seq_len(K), function(k) {
    u <- units_near(kn$i[k], kn$j[k], m_r)
    d2 <- (east[u] - kn$x[k])^2 + (north[u] - kn$y[k])^2
    s <- d2 < r^2
    if (!any(s)) return(NULL)
    cbind(u[s], k, (1 - d2[s] / r^2)^2)
  }))
  B <- sparseMatrix(i = trip[, 1], j = trip[, 2], x = trip[, 3], dims = c(length(east), K))
  if (any(rowSums(B) == 0)) stop("some units are not covered by the basis: increase r_fac")
  B <- Diagonal(x = 1 / rowSums(B)) %*% B
  key <- paste(kn$i, kn$j)
  nbs <- lapply(seq_len(K), function(k) {
    v <- match(paste(kn$i[k] + c(1, -1, 0, 0), kn$j[k] + c(0, 0, 1, -1)), key)
    v[!is.na(v)]
  })
  for (k in which(lengths(nbs) == 0)) {
    d <- (kn$x - kn$x[k])^2 + (kn$y - kn$y[k])^2
    d[k] <- Inf
    j <- which.min(d)
    nbs[[k]] <- j
    nbs[[j]] <- c(nbs[[j]], k)
  }
  A <- sparseMatrix(i = rep(seq_len(K), lengths(nbs)), j = unlist(nbs), x = 1, dims = c(K, K))
  A <- (A + t(A) > 0) * 1
  Wk <- Diagonal(x = 1 / rowSums(A)) %*% A
  ev <- eigen(as.matrix(Wk), only.values = TRUE)$values
  list(B = as(as(B, "CsparseMatrix"), "TsparseMatrix"),
       Wk = as(as(Wk, "CsparseMatrix"), "TsparseMatrix"),
       eig_re = Re(ev), eig_im = Im(ev), knots = kn[, c("x", "y")], K = K, h = h, r = r)
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
#   logdet  - "eigen": log|I - rho*W| from the eigenvalues of W (dense, O(n^3)
#             once); "sparse": from a sparse Cholesky factor of
#             (I - rho*W)'(I - rho*W) at each evaluation; "auto": sparse
#             when n > 5000
#   fixed   - optional named list of parameter blocks held at given values,
#             e.g. list(rho_raw = 0) for the aspatial ZINB-RE, or
#             list(gamma = c(-20, rep(0, q - 1))) for SAR-NB-RE (pi ~ 0)
#   init    - optional named list of starting values (beta, gamma, rho_raw,
#             log_kappa, log_sigma, eta), e.g. from a previous fit
#   zero_spatial - inflation predictor zeta, logit(pi) = zeta:
#             "none"  - restricted model, zeta = Z gamma
#             "field" - zeta = Z gamma + g with a SAR random field
#                       (I - rho0 W) g ~ N(0, sigma0^2 I); adds rho0_raw,
#                       log_sigma0 (fixed) and g (random); TRUE is accepted
#             "sar"   - SAR on the predictor, (I - rho0 W) zeta = Z gamma;
#                       adds rho0_raw only
#             "reduced" - zeta = Z gamma + B alpha with a reduced-rank field on
#                       a knot lattice, (I - rho0 W_K) alpha ~ N(0, sigma0^2 I_K),
#                       rho0 in [0, 1); needs basis from zisar_basis; adds
#                       rho0_raw, log_sigma0 (fixed) and alpha (random)
#   basis   - list from zisar_basis (zero_spatial = "reduced" only)
#
# Returns: S3 object of class "zisar"
# -----------------------------------------------------------------------------
zisar_fit <- function(y, X, Z, W,
                      control = list(eval.max = 1000, iter.max = 500,
                                     rel.tol = 1e-10),
                      silent = TRUE,
                      logdet = c("auto", "eigen", "sparse"),
                      fixed = NULL, init = NULL,
                      zero_spatial = c("none", "field", "sar", "reduced"),
                      basis = NULL) {

  .zisar_compile()
  if (is.logical(zero_spatial)) zero_spatial <- if (zero_spatial) "field" else "none"
  zero_spatial <- match.arg(zero_spatial)
  zmode <- match(zero_spatial, c("none", "field", "sar", "reduced")) - 1L
  if (zmode == 3L && is.null(basis)) stop("zero_spatial = \"reduced\" needs a basis from zisar_basis()")

  n <- length(y)
  p <- ncol(X)
  q <- ncol(Z)
  W_sp <- as(as(W, "CsparseMatrix"), "TsparseMatrix")
  logdet <- match.arg(logdet)
  if (logdet == "auto") logdet <- if (n > 5000) "sparse" else "eigen"
  data <- .zisar_data(y, X, Z, W, logdet, zmode, basis)
  eig_re <- data$eigenvalues_W_re
  eig_im <- data$eigenvalues_W_im

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
    rho0_raw   = 0,
    log_sigma0 = log(0.5),
    eta       = as.numeric(X %*% init_beta),
    g         = if (zmode == 1) rep(0, n) else numeric(0),
    alpha     = if (zmode == 3) rep(0, basis$K) else numeric(0)
  )
  parameters[names(init)] <- init
  parameters[names(fixed)] <- fixed
  map <- lapply(fixed, function(v) factor(rep(NA, length(v))))
  if (zmode == 0) map$rho0_raw <- factor(NA)
  if (!zmode %in% c(1L, 3L)) map$log_sigma0 <- factor(NA)

  # Build TMB objective
  obj <- MakeADFun(
    data       = data,
    parameters = parameters,
    random     = c("eta", if (zmode == 1) "g", if (zmode == 3) "alpha"),
    map        = map,
    DLL        = "zisar_tmb",
    silent     = silent
  )

  # Parameter bounds: beta/gamma unbounded, rho_raw (rho0_raw) in [-4,4],
  # log_kappa in [log(0.01), log(500)], log_sigma (log_sigma0) in
  # [log(0.001), log(10)]; the knot field has rho0 in [0, 1)
  zextra <- c("rho0_raw", "log_sigma0")[seq_len(c(0, 2, 1, 2)[zmode + 1])]
  block <- c(rep("beta", p), rep("gamma", q), "rho_raw", "log_kappa", "log_sigma",
             zextra)
  free  <- !(block %in% names(fixed))
  lower <- c(rep(-Inf, p + q), -4, log(0.01), log(0.001),
             c(if (zmode == 3) 0 else -4, log(0.001))[seq_along(zextra)])[free]
  upper <- c(rep( Inf, p + q),  4, log(500),  log(10),
             c(4, log(10))[seq_along(zextra)])[free]

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

  # Largest gradient component away from the bounds: a Laplace breakdown
  # shows up as a huge gradient at a spurious optimum
  grad <- tryCatch(as.numeric(obj$gr(opt$par)), error = function(e) NA_real_)
  inside <- opt$par > lower + 1e-6 & opt$par < upper - 1e-6
  max_grad <- if (anyNA(grad)) Inf else max(abs(grad[inside]), 0)

  # Hessian of the Laplace objective (fixed effects), shared with sdreport
  H_fixed <- tryCatch(optimHess(opt$par, obj$fn, obj$gr),
                      error = function(e) NULL)

  # Standard errors via sdreport
  sdr <- tryCatch(sdreport(obj, par.fixed = opt$par, hessian.fixed = H_fixed),
                  error = function(e) NULL)

  # Extract estimated parameters (fixed blocks included)
  par_est <- opt$par
  pl <- obj$env$parList(opt$par)
  beta_hat  <- pl$beta
  gamma_hat <- pl$gamma

  rho_hat   <- tanh(pl$rho_raw)
  kappa_hat <- exp(pl$log_kappa)
  sigma_hat <- exp(pl$log_sigma)
  rho0_hat   <- if (zmode > 0) tanh(pl$rho0_raw) else NULL
  sigma0_hat <- if (zmode %in% c(1L, 3L)) exp(pl$log_sigma0) else NULL

  # Estimated latent fields (g = 0 in the restricted model)
  lpb <- obj$env$last.par.best
  eta_hat <- tryCatch(as.numeric(lpb[names(lpb) == "eta"]),
                      error = function(e) as.numeric(X %*% beta_hat))
  alpha_hat <- if (zmode == 3) as.numeric(lpb[names(lpb) == "alpha"]) else NULL
  g_hat <- if (zmode == 1) as.numeric(lpb[names(lpb) == "g"]) else
           if (zmode == 3) as.numeric(basis$B %*% alpha_hat) else rep(0, n)
  # Inflation predictor zeta, logit(pi) = zeta
  zeta_hat <- as.numeric(Z %*% gamma_hat)
  if (zmode %in% c(1L, 3L)) zeta_hat <- zeta_hat + g_hat
  if (zmode == 2) zeta_hat <- as.numeric(solve(Diagonal(n) - rho0_hat * W_sp, zeta_hat))

  # Build coefficient table
  if (!is.null(sdr)) {
    summ_fixed <- summary(sdr, "fixed")
    summ_report <- summary(sdr, "report")

    par_names <- c(paste0("beta_", 0:(p-1)),
                   paste0("gamma_", 0:(q-1)),
                   "rho_raw", "log_kappa", "log_sigma",
                   zextra)[free]
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
    H <- H_fixed
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
        sigma = sigma_hat,
        rho0  = rho0_hat,
        sigma0 = sigma0_hat,
        alpha = alpha_hat
      ),
      eta_hat       = eta_hat,
      g_hat         = g_hat,
      zeta_hat      = zeta_hat,
      zero_spatial  = zero_spatial,
      basis         = basis,
      max_grad      = max_grad,
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
# zisar_objective: TMB objective of the general model at the estimates of a
# restricted fit, with the inflation field at (rho0, sigma0), for profiling
# the Laplace likelihood in (rho0, sigma0) and checking the approximation.
# kappa starts just inside the box when the fit has it at the bound.
# With a basis from zisar_basis the field is the reduced-rank one of
# zero_spatial = "reduced" (random effects alpha instead of g, rho0 >= 0).
# Attributes lower and upper hold the optimisation box of zisar_fit.
# -----------------------------------------------------------------------------
zisar_objective <- function(fit, rho0 = 0, sigma0 = 0.5, silent = TRUE, basis = NULL) {
  .zisar_compile()
  n <- fit$n
  data <- .zisar_data(fit$y, fit$X, fit$Z, fit$W,
                      if (n > 5000) "sparse" else "eigen",
                      if (is.null(basis)) 1L else 3L, basis)
  pr <- fit$params
  parameters <- list(beta = pr$beta, gamma = pr$gamma, rho_raw = atanh(pr$rho),
                     log_kappa = log(min(pr$kappa, 499)), log_sigma = log(pr$sigma),
                     rho0_raw = atanh(rho0), log_sigma0 = log(sigma0),
                     eta = fit$eta_hat,
                     g = if (is.null(basis)) rep(0, n) else numeric(0),
                     alpha = if (is.null(basis)) numeric(0) else rep(0, basis$K))
  obj <- MakeADFun(data, parameters, random = c("eta", if (is.null(basis)) "g" else "alpha"),
                   DLL = "zisar_tmb", silent = silent)
  pq <- fit$p + fit$q
  attr(obj, "lower") <- c(rep(-Inf, pq), -4, log(0.01), log(0.001),
                          if (is.null(basis)) -4 else 0, log(0.001))
  attr(obj, "upper") <- c(rep( Inf, pq),  4, log(500),  log(10),   4, log(10))
  obj
}

# -----------------------------------------------------------------------------
# zisar_field_profile: the inflation field added to a restricted fit. Laplace
# log-likelihood gain over the restricted fit on a grid of (rho0, sigma0), the
# other parameters held at the restricted estimates, then a full optimisation
# from the best grid point. With basis (zisar_basis) the field is the
# reduced-rank one, otherwise the unit-level SAR field. cores: forked
# processes for the grid.
# Returns the grid with the gains, the nlminb result, LR against the restricted
# fit, the field parameters, the largest gradient component, the predicted
# field g_hat and the elapsed minutes.
# -----------------------------------------------------------------------------
zisar_field_profile <- function(fit, basis = NULL,
                                grid = expand.grid(rho0 = c(0, 0.5, 0.8, 0.9, 0.95),
                                                   sigma0 = c(0.1, 0.3, 0.6, 1, 2)),
                                cores = 1, silent = TRUE) {
  t0 <- Sys.time()
  ob <- zisar_objective(fit, silent = silent, basis = basis)
  L <- length(ob$par)
  at <- function(rho0, sigma0) replace(ob$par, c(L - 1, L), c(atanh(rho0), log(sigma0)))
  grid$gain <- -unlist(parallel::mclapply(seq_len(nrow(grid)), function(i)
    ob$fn(at(grid$rho0[i], grid$sigma0[i])), mc.cores = cores, mc.preschedule = FALSE)) - fit$loglik
  best <- grid[which.max(grid$gain), ]
  opt <- nlminb(at(best$rho0, best$sigma0), ob$fn, ob$gr,
                lower = attr(ob, "lower"), upper = attr(ob, "upper"),
                control = list(eval.max = 1000, iter.max = 500, rel.tol = 1e-10))
  gr <- ob$gr(opt$par)
  lp <- ob$env$last.par.best
  pq <- fit$p + fit$q
  g_hat <- if (is.null(basis)) as.numeric(lp[names(lp) == "g"]) else
           as.numeric(basis$B %*% as.numeric(lp[names(lp) == "alpha"]))
  list(grid = grid, start = best, opt = opt, par = opt$par,
       LR = 2 * (-opt$objective - fit$loglik), loglik = -opt$objective,
       rho0 = tanh(opt$par[L - 1]), sigma0 = exp(opt$par[L]),
       rho = tanh(opt$par[pq + 1]), kappa = exp(opt$par[pq + 2]), sigma = exp(opt$par[pq + 3]),
       beta = opt$par[seq_len(fit$p)], gamma = opt$par[fit$p + seq_len(fit$q)],
       convergence = opt$convergence, message = opt$message, max_grad = max(abs(gr)),
       g_hat = g_hat, basis_K = if (is.null(basis)) NA else basis$K,
       elapsed = as.numeric(difftime(Sys.time(), t0, units = "mins")))
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
  cat("  sigma (RE std dev):", round(p$sigma, 4), "\n")
  if (!is.null(p$rho0)) cat("  rho0 (inflation SAR):", round(p$rho0, 4), "\n")
  if (!is.null(p$sigma0)) cat("  sigma0 (inflation RE std dev):", round(p$sigma0, 4), "\n")
  cat("\n")

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
  cat(switch(.zero_mode(object),
             none  = "Inflation: logit(pi_i) = Z'gamma\n",
             field = "Inflation: logit(pi_i) = Z'gamma + g_i, where g = (I-rho0*W)^{-1} e0\n",
             reduced = "Inflation: logit(pi_i) = Z'gamma + (B alpha)_i, where alpha = (I-rho0*W_K)^{-1} e0\n",
             sar   = "Inflation: logit(pi) = (I-rho0*W)^{-1} Z gamma\n"))
  cat("Count:     y_i | s_i=0 ~ NB(exp(eta_i), kappa)\n\n")
  print(object)
  invisible(object)
}

# -----------------------------------------------------------------------------
# zisar_effects: Direct, indirect, total effects
#
# fit   - fitted zisar object
# B     - number of simulation draws for SEs (0 = point estimates only)
#
# The effect of x_r on E[y] = (1 - pi) * exp(eta) is diag(d) S(rho) beta_r,
# d = (1 - pi) * mu, so the averages need only diag(S) and S 1: no dense
# n x n matrix is formed. A covariate that also enters Z (same column name
# in X and Z) adds the inflation channel -diag(mu pi (1 - pi)) S0 gamma_r, with
# S0 = I, or S0 = (I - rho0 W)^{-1} when zero_spatial = "sar", in which case
# the inflation channel has indirect effects too.
# The simulation draws keep the predicted random effects
# u_hat = (I - rho W) eta_hat - X beta, and g_hat when zero_spatial = "field",
# at their estimate.
# -----------------------------------------------------------------------------

# Specification of the inflation predictor of a fit ("none" for fits made
# before the option existed)
.zero_mode <- function(fit) {
  z <- fit$zero_spatial
  if (is.null(z) || identical(z, FALSE)) "none" else if (isTRUE(z)) "field" else z
}

# Inflation predictor zeta at (gamma, rho0), with g at its prediction
.zisar_zeta <- function(fit, gamma = fit$params$gamma, rho0 = fit$params$rho0) {
  zeta <- as.numeric(fit$Z %*% gamma)
  switch(.zero_mode(fit),
         none  = zeta,
         field = zeta + fit$g_hat,
         reduced = zeta + fit$g_hat,
         sar   = as.numeric(solve(Diagonal(fit$n) - rho0 * as(fit$W, "CsparseMatrix"), zeta)))
}

# diag((I - rho W)^{-1}): dense for small n, sparse LU solves in blocks otherwise
.sar_diag <- function(W, rho, block = 2000) {
  n <- nrow(W)
  A <- Diagonal(n) - rho * W
  if (n <= 3000) return(diag(solve(as.matrix(A))))
  luA <- lu(A)
  out <- numeric(n)
  for (s in seq(1, n, by = block)) {
    j <- s:min(n, s + block - 1)
    E <- matrix(0, n, length(j))
    E[cbind(j, seq_along(j))] <- 1
    Sj <- solve(luA, E)
    out[j] <- Sj[cbind(j, seq_along(j))]
  }
  out
}

zisar_effects <- function(fit, B = 500) {
  n <- fit$n
  p <- fit$p
  W <- as(fit$W, "CsparseMatrix")
  rho <- fit$params$rho
  beta <- fit$params$beta
  gamma <- fit$params$gamma
  eta_hat <- fit$eta_hat
  I_n <- Diagonal(n)

  # Fitted values at the mode
  mu_hat <- exp(eta_hat)
  zsar   <- .zero_mode(fit) == "sar"
  rho0   <- fit$params$rho0
  pi_hat <- 1 / (1 + exp(-.zisar_zeta(fit)))
  d_vec  <- (1 - pi_hat) * mu_hat

  # Spatial multiplier S(rho) = (I - rho*W)^{-1}: diagonal and row sums
  s_ii <- .sar_diag(W, rho)
  s_1  <- as.numeric(solve(I_n - rho * W, rep(1, n)))
  # and of S0 for the inflation channel
  s0_ii <- if (zsar) .sar_diag(W, rho0) else 1
  s0_1  <- if (zsar) as.numeric(solve(I_n - rho0 * W, rep(1, n))) else 1

  # Position in Z of each X covariate entering both processes (NA otherwise)
  xn <- colnames(fit$X); zn <- colnames(fit$Z)
  in_z <- if (is.null(xn) || is.null(zn)) rep(NA_integer_, p) else match(xn, zn)
  in_z[seq_len(p) == 1] <- NA_integer_          # intercepts
  g_x <- function(g) ifelse(is.na(in_z), 0, g[in_z])

  zch <- mu_hat * pi_hat * (1 - pi_hat)         # d E[y] / d (z'gamma), negated
  effects <- list()
  for (r in seq_len(p)) {
    zero_i   <- -zch * s0_ii * g_x(gamma)[r]
    zero_t   <- -zch * s0_1 * g_x(gamma)[r]
    direct_i <- d_vec * s_ii * beta[r] + zero_i
    total_i  <- d_vec * s_1 * beta[r] + zero_t
    effects[[r]] <- list(
      covariate    = r,
      name         = if (is.null(xn)) NULL else xn[r],
      avg_direct   = mean(direct_i),
      avg_indirect = mean(total_i - direct_i),
      avg_total    = mean(total_i),
      avg_zero     = mean(zero_i),
      avg_zero_indirect = mean(zero_t - zero_i)
    )
  }

  # Simulation-based SEs
  if (B > 0 && all(is.finite(fit$vcov))) {
    sim_effects <- array(NA_real_, dim = c(B, p, 3))
    theta_draws <- MASS::mvrnorm(B, fit$coefficients, fit$vcov)
    u_hat <- as.numeric((I_n - rho * W) %*% eta_hat - fit$X %*% beta)

    idx_rho <- p + fit$q + 1
    rho_draws <- tanh(theta_draws[, idx_rho])
    rho0_draws <- if (zsar) tanh(theta_draws[, names(fit$coefficients) == "rho0_raw"])

    # diag(S) along the draws: exact at Chebyshev nodes spanning the draws,
    # Lagrange interpolation in between (S_ii is smooth in rho)
    diag_at <- function(draws) {
      if (n <= 3000) return(function(r0) .sar_diag(W, r0))
      a <- min(draws); b <- max(draws)
      nodes <- (a + b) / 2 + (b - a) / 2 * cos((2 * (1:7) - 1) * pi / 14)
      s_nodes <- sapply(nodes, function(r0) .sar_diag(W, r0))
      lagrange <- function(x) sapply(seq_along(nodes), function(k)
        prod((x - nodes[-k]) / (nodes[k] - nodes[-k])))
      function(r0) as.numeric(s_nodes %*% lagrange(r0))
    }
    s_ii_at <- diag_at(rho_draws)
    s0_ii_at <- if (zsar) diag_at(rho0_draws)

    for (b in seq_len(B)) {
      beta_b  <- theta_draws[b, 1:p]
      gamma_b <- theta_draws[b, p + seq_len(fit$q)]
      rho_b   <- rho_draws[b]
      rho0_b  <- if (zsar) rho0_draws[b]
      if (abs(rho_b) >= 1 || (zsar && abs(rho0_b) >= 1)) next

      A_b   <- I_n - rho_b * W
      eta_b <- as.numeric(solve(A_b, fit$X %*% beta_b + u_hat))
      pi_b  <- 1 / (1 + exp(-.zisar_zeta(fit, gamma_b, rho0_b)))
      d_b   <- (1 - pi_b) * exp(eta_b)
      dS_b  <- mean(d_b * s_ii_at(rho_b))
      dT_b  <- mean(d_b * as.numeric(solve(A_b, rep(1, n))))
      zc_b  <- exp(eta_b) * pi_b * (1 - pi_b)
      zS_b  <- if (zsar) mean(zc_b * s0_ii_at(rho0_b)) else mean(zc_b)
      zT_b  <- if (zsar) mean(zc_b * as.numeric(solve(I_n - rho0_b * W, rep(1, n)))) else mean(zc_b)

      sim_effects[b, , 1] <- beta_b * dS_b - zS_b * g_x(gamma_b)
      sim_effects[b, , 3] <- beta_b * dT_b - zT_b * g_x(gamma_b)
      sim_effects[b, , 2] <- beta_b * (dT_b - dS_b) - (zT_b - zS_b) * g_x(gamma_b)
    }

    for (r in seq_len(p)) {
      valid <- is.finite(sim_effects[, r, 1])
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
    cat(sprintf("Covariate %d%s:\n", eff$covariate,
                if (is.null(eff$name)) "" else paste0(" (", eff$name, ")")))
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
  pi_hat <- 1 / (1 + exp(-.zisar_zeta(fit)))

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
# zero_spatial - inflation predictor as in zisar_fit: "none", "field" (SAR
#           field with rho0, sigma0) or "sar" (SAR filter with rho0); "none"
#           leaves the random-number stream unchanged
# -----------------------------------------------------------------------------
rzisar <- function(n_side, beta, rho, gamma, kappa = 2, sigma = 0.5,
                   seed = NULL, zero_spatial = c("none", "field", "sar"),
                   rho0 = 0, sigma0 = 0) {
  zero_spatial <- match.arg(zero_spatial)
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
  zeta <- as.numeric(Z %*% gamma)
  if (zero_spatial == "field")
    zeta <- zeta + as.numeric(solve(I_n - rho0 * W_sp, rnorm(n, 0, sigma0)))
  if (zero_spatial == "sar")
    zeta <- as.numeric(solve(I_n - rho0 * W_sp, zeta))
  pi_i <- 1 / (1 + exp(-zeta))
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
                  kappa = kappa, sigma = sigma, rho0 = rho0, sigma0 = sigma0),
    eta_true = eta,
    mu_true  = mu,
    pi_true  = as.numeric(pi_i),
    s_true   = s_i
  )
}
