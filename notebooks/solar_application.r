# =============================================================================
#  ZISAR application: domestic solar PV adoption across English LSOAs
#
#  y_i : domestic PV installations of up to 10 kW commissioned under the
#        Feed-in Tariff after the tariff reform of 8 February 2016 and before
#        the scheme closed to new entrants (31 March 2019), in LSOA i
#        (n = 32,322 harmonised 2001/2011 LSOAs, see solar_data.r)
#  X   : count-process covariates (2011 Census, IoD 2019, VOA 2015, DESNZ
#        2010, location), including the covariates of Z
#  Z   : structural-zero covariates: housing form and age, density, tenure,
#        communal and vacant dwellings, and barriers measured before the
#        reform (EPC roof form and floor area, conservation areas, listed
#        buildings, leasehold houses)
#  W   : queen contiguity, row-standardised (k-NN in the sensitivity check)
#  The cumulative 2010-2019 count, the post-reform count of all capacities
#  and of up to 4 kW, and the specification without the pre-reform barriers
#  are robustness checks; a score-type test checks the restricted model for
#  spatial dependence in the structural regime; the general ZISAR (SAR field
#  in the inflation predictor), neighbour averages WZ in the inflation
#  predictor and a SAR filter on it are sensitivity checks.
#
#  Run solar_data.r first. Writes solar_results.txt and solar_results.rds,
#  and the figure ../tex/figures/solar_maps.png
# =============================================================================

source("zisar.r")
suppressPackageStartupMessages({
  library(pscl); library(spdep); library(sf); library(ggplot2); library(patchwork)
})
set.seed(2026)

solar <- readRDS("../data/solar/solar_lsoa.rds")
Ws    <- readRDS("../data/solar/W_lsoa.rds")
W     <- Ws$queen
n     <- nrow(solar)
y     <- solar$pv_dom_late

std <- function(x) (x - mean(x)) / sd(x)
cov_x <- with(solar, data.frame(
  households = std(log(households)),    # log number of households
  detached   = std(sh_detached),        # detached-house share
  social     = std(sh_social),          # social-rented share (landlord programmes)
  income_dep = std(income_dep),         # IoD income deprivation rate
  managerial = std(sh_nssec12),         # NS-SEC 1-2: managerial and professional
  age65      = std(sh_age65),           # share aged 65+
  children   = std(sh_children),        # households with dependent children
  offgas     = std(sh_offgas),          # households without gas central heating
  band_eh    = std(sh_band_eh),         # Council Tax bands E-H (VOA 2015)
  elec       = std(log(elec_kwh)),      # log domestic kWh per meter, 2010
  north      = std(north / 1000)))      # northing (km): solar resource gradient
cov_z <- with(solar, data.frame(
  flat       = std(sh_flat),            # share of household spaces in flats
  density    = std(log(density)),       # log residents per km^2
  private    = std(sh_private),         # private-rented share
  communal   = std(sh_communal),        # residents in communal establishments
  empty      = std(sh_empty),           # household spaces with no usual residents
  pre1919    = std(sh_pre1919)))        # dwellings built before 1919 (VOA 2015)
# Barriers to installation measured before the reform; units without an EPC
# dwelling or a house sale to measure them on (a few dozen) take the mean
imp <- function(v) ifelse(is.na(v), mean(v, na.rm = TRUE), v)
cov_b <- with(solar, data.frame(
  no_roof      = std(imp(epc_noroof)),     # EPC dwellings with a dwelling or premises above
  flat_roof    = std(imp(epc_flatroof)),   # EPC dwellings with a flat roof
  roof_room    = std(imp(epc_roofroom)),   # EPC dwellings with rooms in the roof
  floor_area   = std(imp(epc_floor)),      # log mean floor area of EPC dwellings
  epc_cover    = std(log(imp(epc_cover))), # EPC dwellings per household (log), coverage
  conservation = std(imp(ca_share)),       # EPC dwellings in a conservation area
  listed       = std(log1p(listed_rate)),  # listed buildings per 1,000 households (log1p)
  leasehold    = std(imp(lease_house))))   # leasehold share of houses sold 2008-2015
# Housing form also shifts the intensity among eligible LSOAs: the Z covariates
# enter X too (same names, so zisar_effects adds the inflation channel)
X <- cbind(intercept = 1, as.matrix(cov_x), as.matrix(cov_z), as.matrix(cov_b))
Z <- cbind(intercept = 1, as.matrix(cov_z), as.matrix(cov_b))
p <- ncol(X); q <- ncol(Z)
# specification without the pre-reform barriers (robustness)
X0 <- X[, setdiff(colnames(X), colnames(cov_b))]; Z0 <- Z[, setdiff(colnames(Z), colnames(cov_b))]
df <- data.frame(y = y, cov_x, cov_z, cov_b)
f_count <- paste(c(colnames(cov_x), colnames(cov_z), colnames(cov_b)), collapse = " + ")
f_zero  <- paste(c(colnames(cov_z), colnames(cov_b)), collapse = " + ")

lw <- nb2listw(Ws$nb_queen, style = "W")
timing <- list()

# -----------------------------------------------------------------------------
# 1. Benchmarks
# -----------------------------------------------------------------------------
fit_nb <- MASS::glm.nb(as.formula(paste("y ~", f_count)), data = df)

fit_zinb <- zeroinfl(as.formula(paste("y ~", f_count, "|", f_zero)),
                     data = df, dist = "negbin")

# SAR-NB: NB GLM with the spatial lag W y instrumented by [X, WX, W^2X]
WX  <- as.matrix(W %*% X[, -1]); W2X <- as.matrix(W %*% WX)
H   <- cbind(X, WX, W2X)
Wy  <- as.numeric(W %*% y)
df$Wy_hat <- as.numeric(H %*% solve(crossprod(H), crossprod(H, Wy)))
fit_sarnb <- MASS::glm.nb(as.formula(paste("y ~", f_count, "+ Wy_hat")), data = df)

# -----------------------------------------------------------------------------
# 2. ZISAR, its nested restrictions, the alternative W and the cumulative
#    and capacity-restricted outcomes, and the specification without the
#    pre-reform barriers (same template, sparse log-determinant). The nine
#    fits run as parallel processes and are cached in solar_fits.rds.
# -----------------------------------------------------------------------------
fits_file <- "../data/solar/solar_fits.rds"
fits <- if (file.exists(fits_file)) readRDS(fits_file) else list()
specs <- list(
  zisar   = list(y = y, W = W),
  # rho = 0: aspatial ZINB with Gaussian random effects (LR test of rho)
  zinbre  = list(y = y, W = W, fixed = list(rho_raw = 0)),
  # pi = 0: spatial NB with random effects, no zero-inflation
  sarnbre = list(y = y, W = W, fixed = list(gamma = c(-20, rep(0, q - 1)))),
  k6      = list(y = y, W = Ws$k6),
  k10     = list(y = y, W = Ws$k10),
  cum     = list(y = solar$pv_dom, W = W),
  # robustness: all capacities, up to 4 kW, no pre-reform barriers
  all     = list(y = solar$pv_dom_late_all, W = W),
  kw4     = list(y = solar$pv_dom_late_4, W = W),
  nobar   = list(y = y, W = W, X = X0, Z = Z0))
fit_spec <- function(nm, init = NULL) {
  s <- specs[[nm]]
  t0 <- Sys.time()
  f <- zisar_fit(s$y, if (is.null(s$X)) X else s$X, if (is.null(s$Z)) Z else s$Z,
                 s$W, fixed = s$fixed, init = init)
  f$elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  f$tmb_obj <- NULL                     # external pointers do not survive the fork
  f
}
# Only the specifications missing from the cache are fitted
todo <- setdiff(names(specs), names(fits))
if (length(todo)) {
  .zisar_compile()
  fits[todo] <- parallel::mclapply(todo, fit_spec, mc.cores = length(todo), mc.preschedule = FALSE)
  stopifnot(all(sapply(fits[todo], inherits, "zisar")))
  # A fit of the post-reform outcome that stops with a false-convergence flag
  # (a spurious optimum of the Laplace objective, with a huge gradient) is
  # restarted from the converged main estimates and, if it stops again, from
  # that stopping point, in both cases with kappa near the upper bound where
  # all the converged fits of this application end; the main fit itself is
  # restarted from its stopping point
  warm <- function(f0, f) {
    b <- f0$params$beta[match(colnames(f$X), colnames(f0$X))]
    g <- f0$params$gamma[match(colnames(f$Z), colnames(f0$Z))]
    list(beta = ifelse(is.na(b), 0, b), gamma = ifelse(is.na(g), 0, g),
         rho_raw = atanh(f0$params$rho), log_kappa = log(499),
         log_sigma = log(f0$params$sigma), eta = f0$eta_hat)
  }
  restart <- function(f, nm, from) {
    r <- fit_spec(nm, init = warm(from, f))
    r$elapsed <- r$elapsed + f$elapsed
    r$restarted <- TRUE
    r
  }
  if ("zisar" %in% todo && fits$zisar$convergence != 0)
    fits$zisar <- restart(fits$zisar, "zisar", fits$zisar)
  redo <- Filter(function(nm) fits[[nm]]$convergence != 0, setdiff(todo, c("zisar", "cum")))
  if (length(redo))
    fits[redo] <- parallel::mclapply(redo, function(nm) {
      f <- restart(fits[[nm]], nm, fits$zisar)
      if (f$convergence != 0) f <- restart(f, nm, f)
      f
    }, mc.cores = length(redo), mc.preschedule = FALSE)
  saveRDS(fits, fits_file)
}
fit_zisar   <- fits$zisar
fit_zinbre  <- fits$zinbre
fit_sarnbre <- fits$sarnbre
timing$zisar <- fit_zisar$elapsed

lr_rho <- 2 * (fit_zisar$loglik - fit_zinbre$loglik)
wald_rho <- fit_zisar$coef_table["rho_raw", ]

# -----------------------------------------------------------------------------
# 3. Effects decomposition
# -----------------------------------------------------------------------------
timing$effects <- system.time(eff <- zisar_effects(fit_zisar, B = 500))["elapsed"]

# -----------------------------------------------------------------------------
# 4. Fit of the count distribution: observed vs expected frequencies
#    (ZISAR expected frequencies integrate the latent field by simulation
#     from the fitted SAR prior, i.e. they are marginal, like the others)
# -----------------------------------------------------------------------------
K <- 0:5
obs_freq <- sapply(K, function(k) sum(y == k))

nb_pmf <- function(k, mu, size) dnbinom(k, size = size, mu = mu)
exp_nb <- sapply(K, function(k) sum(nb_pmf(k, fitted(fit_nb), fit_nb$theta)))
exp_zinb <- colSums(predict(fit_zinb, type = "prob")[, K + 1])
exp_sarnb <- sapply(K, function(k) sum(nb_pmf(k, fitted(fit_sarnb), fit_sarnb$theta)))

zisar_marginal_freq <- function(fit, K, R = 200) {
  A <- Diagonal(fit$n) - fit$params$rho * as(fit$W, "CsparseMatrix")
  luA <- lu(A)
  Xb <- as.numeric(fit$X %*% fit$params$beta)
  pi_i <- plogis(as.numeric(fit$Z %*% fit$params$gamma))
  kap <- fit$params$kappa
  acc <- numeric(length(K))
  for (r in seq_len(R)) {
    eta <- as.numeric(solve(luA, Xb + rnorm(fit$n, 0, fit$params$sigma)))
    for (j in seq_along(K))
      acc[j] <- acc[j] + sum((K[j] == 0) * pi_i +
                             (1 - pi_i) * dnbinom(K[j], size = kap, mu = exp(eta)))
  }
  acc / R
}
exp_zisar <- zisar_marginal_freq(fit_zisar, K)
freq_tab <- rbind(observed = obs_freq, NB = exp_nb, ZINB = exp_zinb,
                  `SAR-NB` = exp_sarnb, ZISAR = exp_zisar)
colnames(freq_tab) <- K

# Posterior probability that an observed zero is structural (ZISAR, at the mode)
pi_hat <- zisar_predict(fit_zisar, type = "zero")
mu_hat <- zisar_predict(fit_zisar, type = "count")
p0_nb  <- (fit_zisar$params$kappa / (fit_zisar$params$kappa + mu_hat))^fit_zisar$params$kappa
post_struct <- ifelse(y == 0, pi_hat / (pi_hat + (1 - pi_hat) * p0_nb), 0)

# -----------------------------------------------------------------------------
# 5. Residual spatial autocorrelation (Pearson residuals)
# -----------------------------------------------------------------------------
pearson_zisar <- function(fit) {
  m <- zisar_predict(fit, type = "count"); pz <- zisar_predict(fit, type = "zero")
  v <- (1 - pz) * m * (1 + m / fit$params$kappa + pz * m)
  (fit$y - (1 - pz) * m) / sqrt(v)
}
moran <- list(
  y      = moran.test(y, lw),
  NB     = moran.test(residuals(fit_nb, type = "pearson"), lw),
  ZINB   = moran.test(residuals(fit_zinb, type = "pearson"), lw),
  SARNB  = moran.test(residuals(fit_sarnb, type = "pearson"), lw),
  ZISAR  = moran.test(pearson_zisar(fit_zisar), lw))

# -----------------------------------------------------------------------------
# 6. Sensitivity to W
# -----------------------------------------------------------------------------
sens <- fits[c("k6", "k10")]

# -----------------------------------------------------------------------------
# 6b. Robustness: cumulative count over the whole scheme (2010-2019)
# -----------------------------------------------------------------------------
y_cum  <- solar$pv_dom
df_cum <- transform(df, y = y_cum)
rob <- list(
  nb    = MASS::glm.nb(as.formula(paste("y ~", f_count)), data = df_cum),
  zinb  = zeroinfl(as.formula(paste("y ~", f_count, "|", f_zero)),
                   data = df_cum, dist = "negbin"),
  zisar = fits$cum)

# With counts in the hundreds the optimiser stops with a false-convergence flag.
# Restart twice from the stopping point, the second time with kappa held at its
# value there; the random-number stream is restored so that the draws used
# below are unaffected.
if (is.null(fits$cum_restarts)) {
  seed <- .Random.seed
  start_from <- function(f) list(beta = f$params$beta, gamma = f$params$gamma,
    rho_raw = atanh(f$params$rho), log_kappa = log(f$params$kappa),
    log_sigma = log(f$params$sigma), eta = f$eta_hat)
  r1 <- zisar_fit(y_cum, X, Z, W, init = start_from(fits$cum))
  i2 <- start_from(r1); i2$log_kappa <- NULL
  r2 <- zisar_fit(y_cum, X, Z, W, init = i2,
                  fixed = list(log_kappa = log(r1$params$kappa)))
  fits$cum_restarts <- lapply(list(r1, r2), `[`,
    c("params", "coef_table", "loglik", "AIC", "convergence"))
  saveRDS(fits, fits_file)
  assign(".Random.seed", seed, envir = globalenv())
}

# -----------------------------------------------------------------------------
# 6c. Sensitivity: general ZISAR, logit(pi) = Z gamma + g with its own SAR
#     field g. On these data its Laplace approximation breaks down, so the
#     check records (i) the Laplace log-likelihood over a (rho0, sigma0) grid,
#     with the other parameters at the restricted estimates, and (ii) one
#     optimiser run from the ridge of that grid: where it stops, the AD and
#     finite-difference gradients there, and a scan of the objective in
#     log(sigma0) around the stopping point
# -----------------------------------------------------------------------------
if (is.null(fits$general)) {
  seed <- .Random.seed
  ob <- zisar_objective(fit_zisar, rho0 = 0.9, sigma0 = 0.8)
  L <- length(ob$par)
  grid <- expand.grid(rho0 = c(0, 0.5, 0.8, 0.9, 0.95, 0.98, 0.995),
                      sigma0 = c(0.02, 0.05, 0.1, 0.2, 0.35, 0.6, 1))
  grid$gain <- -unlist(parallel::mclapply(seq_len(nrow(grid)), function(i)
    ob$fn(replace(ob$par, c(L - 1, L), c(atanh(grid$rho0[i]), log(grid$sigma0[i])))),
    mc.cores = 14, mc.preschedule = FALSE)) - fit_zisar$loglik
  t0 <- Sys.time()
  opt <- nlminb(ob$par, ob$fn, ob$gr, lower = attr(ob, "lower"), upper = attr(ob, "upper"),
                control = list(eval.max = 1000, iter.max = 500, rel.tol = 1e-10))
  pnames <- c(paste0("beta:", colnames(X)), paste0("gamma:", colnames(Z)),
              "rho_raw", "log_kappa", "log_sigma", "rho0_raw", "log_sigma0")
  gr <- ob$gr(opt$par); j <- which.max(abs(gr))
  h <- 1e-4; e <- replace(numeric(L), j, h)
  scan_t <- seq(-2e-3, 2e-3, length.out = 11)
  fits$general <- list(
    grid = grid, opt = opt, par = setNames(opt$par, pnames),
    elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs")),
    grad_which = pnames[j], grad_ad = gr[j],
    grad_fd = (ob$fn(opt$par + e) - ob$fn(opt$par - e)) / (2 * h),
    scan_t = scan_t,
    scan_f = sapply(scan_t, function(t) ob$fn(replace(opt$par, L, opt$par[L] + t))))
  saveRDS(fits, fits_file)
  assign(".Random.seed", seed, envir = globalenv())
}
gen <- fits$general

# -----------------------------------------------------------------------------
# 6d. Sensitivity: spatial inflation equation at the LSOA level without new
#     random effects: (i) neighbour averages WZ of the zero covariates added
#     to Z; (ii) a SAR filter on the inflation predictor,
#     logit(pi) = (I - rho0 W)^{-1} Z gamma. Both start from the restricted
#     estimates. With the fourteen inflation covariates the WZ fit stops at a
#     spurious optimum (false convergence, gradient ~1e8) from this start, from
#     its own stopping point and from the main estimates with kappa near the
#     bound, at a log-likelihood that differs between starts: it is reported
#     as not estimable
# -----------------------------------------------------------------------------
WZ <- as.matrix(W %*% Z[, -1]); colnames(WZ) <- paste0("W_", colnames(Z)[-1])
if (is.null(fits[["sar"]])) {
  seed <- .Random.seed
  pr <- fit_zisar$params
  start <- list(beta = pr$beta, rho_raw = atanh(pr$rho), log_kappa = log(min(pr$kappa, 499)),
                log_sigma = log(pr$sigma), eta = fit_zisar$eta_hat)
  specs <- list(
    slx = list(Z = cbind(Z, WZ), zero_spatial = "none",
               init = c(start, list(gamma = c(pr$gamma, rep(0, ncol(WZ)))))),
    sar = list(Z = Z, zero_spatial = "sar",
               init = c(start, list(gamma = pr$gamma, rho0_raw = 0))))
  fits[names(specs)] <- parallel::mclapply(specs, function(s) {
    t0 <- Sys.time()
    f <- zisar_fit(y, X, s$Z, W, init = s$init, zero_spatial = s$zero_spatial)
    f$elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    f$tmb_obj <- NULL
    f
  }, mc.cores = 2, mc.preschedule = FALSE)
  stopifnot(all(sapply(fits[names(specs)], inherits, "zisar")))
  saveRDS(fits, fits_file)
  assign(".Random.seed", seed, envir = globalenv())
}
fit_slx <- fits[["slx"]]   # [[ ]]: fits$sar would match fits$sarnbre
fit_sar <- fits[["sar"]]

# Validation of the posterior structural-zero probability against adoption
# before the reform, which neither model uses: mean probability by number of
# pre-reform installations, and AUC for the zeros with none
pre <- solar$pv_dom - y
validate_struct <- function(fit) {
  pz <- zisar_predict(fit, type = "zero"); m <- zisar_predict(fit, type = "count")
  p0 <- (fit$params$kappa / (fit$params$kappa + m))^fit$params$kappa
  ps <- ifelse(y == 0, pz / (pz + (1 - pz) * p0), 0)
  a <- ps[y == 0 & pre == 0]; b <- ps[y == 0 & pre > 0]
  c(none = mean(a), `1-4` = mean(ps[y == 0 & pre %in% 1:4]), `5+` = mean(ps[y == 0 & pre >= 5]),
    AUC = mean(outer(a, b, ">") + 0.5 * outer(a, b, "==")), expected = sum(ps))
}
validation <- validate_struct(fit_zisar)

# Posterior predictive check on the clustering of zeros: Moran's I of
# 1{y = 0} in the data and in R draws from the fitted model, with the latent
# fields drawn from their SAR priors (England and Greater London, queen W).
# The random-number stream is restored afterwards.
london <- grepl("^E09", solar$lad19cd)
W_lon <- W[london, london]
W_lon <- Diagonal(x = 1 / pmax(rowSums(W_lon), 1e-12)) %*% W_lon
moran_cols <- function(M, Wm) {
  M <- scale(M, scale = FALSE)
  colSums(M * as.matrix(Wm %*% M)) / colSums(M^2)
}
zero_ppc <- function(fit, R = 200) {
  Wc <- as(fit$W, "CsparseMatrix"); pr <- fit$params; nn <- fit$n
  E <- as.matrix(solve(Diagonal(nn) - pr$rho * Wc,
                       as.numeric(fit$X %*% pr$beta) + matrix(rnorm(nn * R, 0, pr$sigma), nn, R)))
  field <- .zero_mode(fit) == "field"
  G <- if (field)
    as.matrix(solve(Diagonal(nn) - pr$rho0 * Wc, matrix(rnorm(nn * R, 0, pr$sigma0), nn, R))) else 0
  P <- plogis((if (field) as.numeric(fit$Z %*% pr$gamma) else .zisar_zeta(fit)) + G)
  Y0 <- 1 * (matrix(rbinom(nn * R, 1, P), nn, R) == 1 |
             matrix(rnbinom(nn * R, size = pr$kappa, mu = exp(E)), nn, R) == 0)
  list(zeros = colSums(Y0), zeros_lon = colSums(Y0[london, ]),
       moran = moran_cols(Y0, Wc), moran_lon = moran_cols(Y0[london, ], W_lon))
}
seed <- .Random.seed
ppc <- list(observed = list(zeros = sum(y == 0), zeros_lon = sum(y[london] == 0),
                            moran = moran_cols(matrix(1 * (y == 0)), W),
                            moran_lon = moran_cols(matrix(1 * (y[london] == 0)), W_lon)),
            zisar = zero_ppc(fit_zisar), slx = zero_ppc(fit_slx), sar = zero_ppc(fit_sar))
assign(".Random.seed", seed, envir = globalenv())

# -----------------------------------------------------------------------------
# 6e. Score-type test of spatial dependence in the structural regime. Under the
#     restricted model the score of the inflation predictor zeta_i is
#     r_i = P(s_i = 1 | y_i) - pi_i, and the score test of sigma0 = 0 against
#     a SAR field g is a quadratic form in r; the statistic is Moran's I of r
#     (England, Greater London, outside London, each with its own row-
#     standardised W). Null distribution: parametric bootstrap from the fitted
#     restricted model, with the parameters held at their estimates and the
#     inner Laplace mode re-solved on each draw. Main specification, all
#     capacities, up to 4 kW and without the pre-reform barriers.
# -----------------------------------------------------------------------------
W_out <- W[!london, !london]
W_out <- Diagonal(x = 1 / pmax(rowSums(W_out), 1e-12)) %*% W_out
struct_test <- function(fit, B = 200, seed = 9000) {
  .zisar_compile()
  pr <- fit$params; nn <- fit$n; Wc <- as(fit$W, "CsparseMatrix")
  pz <- plogis(as.numeric(fit$Z %*% pr$gamma)); Xb <- as.numeric(fit$X %*% pr$beta)
  data0 <- .zisar_data(fit$y, fit$X, fit$Z, fit$W, "sparse", 0L)
  eta_mode <- function(yy) {
    d <- data0; d$y <- yy
    par <- list(beta = pr$beta, gamma = pr$gamma, rho_raw = atanh(pr$rho), log_kappa = log(pr$kappa),
                log_sigma = log(pr$sigma), rho0_raw = 0, log_sigma0 = log(0.5), eta = Xb, g = numeric(0))
    map <- lapply(par[c("beta", "gamma", "rho_raw", "log_kappa", "log_sigma", "rho0_raw", "log_sigma0")],
                  function(v) factor(rep(NA, length(v))))
    ob <- MakeADFun(d, par, random = "eta", map = map, DLL = "zisar_tmb", silent = TRUE)
    ob$fn(ob$par); lp <- ob$env$last.par.best; as.numeric(lp[names(lp) == "eta"])
  }
  resid <- function(yy, eta) {
    p0 <- (pr$kappa / (pr$kappa + exp(eta)))^pr$kappa
    ifelse(yy == 0, pz / (pz + (1 - pz) * p0), 0) - pz
  }
  stat <- function(r) c(england = moran_cols(matrix(r), Wc), london = moran_cols(matrix(r[london]), W_lon),
                        outside = moran_cols(matrix(r[!london]), W_out))
  obs <- stat(resid(fit$y, fit$eta_hat))
  sim <- do.call(rbind, parallel::mclapply(seq_len(B), function(b) {
    set.seed(seed + b)
    eta <- as.numeric(solve(Diagonal(nn) - pr$rho * Wc, Xb + rnorm(nn, 0, pr$sigma)))
    yy <- ifelse(rbinom(nn, 1, pz) == 1, 0L, rnbinom(nn, size = pr$kappa, mu = exp(eta)))
    stat(resid(yy, eta_mode(yy)))
  }, mc.cores = 12, mc.preschedule = FALSE))
  list(obs = obs, null_mean = colMeans(sim), null_sd = apply(sim, 2, sd),
       z = (obs - colMeans(sim)) / apply(sim, 2, sd),
       p = sapply(names(obs), function(k) (1 + sum(sim[, k] >= obs[k])) / (B + 1)), B = B)
}
if (is.null(fits$struct)) {
  seed <- .Random.seed
  fits$struct <- lapply(list(zisar = fit_zisar, all = fits$all, kw4 = fits$kw4, nobar = fits$nobar),
                        struct_test)
  saveRDS(fits, fits_file)
  assign(".Random.seed", seed, envir = globalenv())
}
struct <- fits$struct

# -----------------------------------------------------------------------------
# 7. Results
# -----------------------------------------------------------------------------
saveRDS(list(fit_nb = fit_nb, fit_zinb = fit_zinb, fit_sarnb = fit_sarnb,
             zisar = fit_zisar[setdiff(names(fit_zisar), c("tmb_obj", "W", "X", "Z"))],
             zinbre = fit_zinbre[c("params", "coef_table", "loglik", "AIC", "BIC")],
             sarnbre = fit_sarnbre[c("params", "coef_table", "loglik", "AIC", "BIC")],
             sens = lapply(sens, `[`, c("params", "coef_table", "transformed", "loglik", "AIC")),
             eff = eff, freq_tab = freq_tab, moran = moran, timing = timing,
             pi_hat = pi_hat, mu_hat = mu_hat, post_struct = post_struct,
             rob = list(nb = rob$nb, zinb = rob$zinb,
                        zisar = rob$zisar[c("params", "coef_table", "transformed",
                                            "loglik", "AIC", "BIC", "convergence")]),
             general = gen,
             ext = lapply(list(slx = fit_slx, sar = fit_sar), `[`,
                          c("params", "coef_table", "transformed", "loglik", "AIC", "BIC",
                            "convergence", "max_grad", "elapsed")),
             robkw = lapply(fits[c("all", "kw4", "nobar")], `[`,
                            c("params", "coef_table", "transformed", "loglik", "AIC", "BIC",
                              "convergence", "max_grad", "elapsed")),
             validation = validation, ppc = ppc, struct = struct),
        "solar_results.rds")

ll <- function(f) as.numeric(logLik(f))
k_of <- function(f) attr(logLik(f), "df")
info <- rbind(
  NB      = c(ll(fit_nb),  k_of(fit_nb)),
  ZINB    = c(ll(fit_zinb), k_of(fit_zinb)),
  `SAR-NB` = c(ll(fit_sarnb), k_of(fit_sarnb)),
  `ZINB-RE` = c(fit_zinbre$loglik, length(fit_zinbre$coefficients)),
  `SAR-NB-RE` = c(fit_sarnbre$loglik, length(fit_sarnbre$coefficients)),
  ZISAR   = c(fit_zisar$loglik, length(fit_zisar$coefficients)))
info <- cbind(info, AIC = 2 * info[, 2] - 2 * info[, 1],
              BIC = log(n) * info[, 2] - 2 * info[, 1])
colnames(info)[1:2] <- c("logLik", "k")

sink("solar_results.txt")
cat("Domestic PV up to 10 kW installed under the FiT, 8 Feb 2016 - 31 Mar 2019, per LSOA, England\n")
cat("n =", n, " installations =", sum(y), " zeros =", sum(y == 0),
    sprintf("(%.2f%%)", 100 * mean(y == 0)), " ones =", sum(y == 1),
    " mean =", round(mean(y), 2), " var =", round(var(y), 1),
    " median =", median(y), " max =", max(y), "\n")
cat("W queen: mean neighbours", round(mean(card(lw$neighbours)), 2),
    " islands relinked:", length(Ws$islands), "\n\n")
cat("--- Descriptives (raw covariates) ---\n")
print(round(t(sapply(transform(solar, floor_m2 = exp(epc_floor))[, c("pv_dom_late", "pv_dom_late_all", "pv_dom_late_4",
  "pv_dom", "households", "sh_detached",
  "sh_social", "income_dep", "sh_nssec12", "sh_age65", "sh_children", "sh_offgas", "sh_band_eh", "elec_kwh", "north",
  "sh_flat", "density", "sh_private", "sh_communal", "sh_empty", "sh_pre1919",
  "epc_noroof", "epc_flatroof", "epc_roofroom", "floor_m2", "epc_cover", "ca_share", "listed_rate", "lease_house")],
  function(v) c(mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE), min = min(v, na.rm = TRUE),
                median = median(v, na.rm = TRUE), max = max(v, na.rm = TRUE), missing = sum(is.na(v))))), 3))
cat("\n--- Correlations of standardised covariates ---\n")
print(round(cor(cbind(cov_x, cov_z, cov_b)), 2))
cat("\n--- NB ---\n"); print(summary(fit_nb)$coefficients); cat("theta =", fit_nb$theta, "SE", fit_nb$SE.theta, "\n")
cat("\n--- ZINB ---\n"); print(summary(fit_zinb)$coefficients); cat("theta =", fit_zinb$theta, "\n")
cat("\n--- SAR-NB ---\n"); print(summary(fit_sarnb)$coefficients); cat("theta =", fit_sarnb$theta, "\n")
cat("\n--- ZISAR (convergence", fit_zisar$convergence, ", elapsed", round(timing$zisar, 1), "s) ---\n")
print(fit_zisar$coef_table); print(fit_zisar$transformed)
cat("\n--- ZINB-RE (rho = 0) ---\n"); print(fit_zinbre$coef_table); print(unlist(fit_zinbre$params[c("kappa", "sigma")]))
cat("\n--- SAR-NB-RE (pi = 0) ---\n"); print(fit_sarnbre$coef_table); print(fit_sarnbre$transformed)
cat("\n--- Information criteria ---\n"); print(round(info, 1))
cat(sprintf("\nLR test rho = 0: %.2f (p = %.3g); Wald z on atanh(rho): %.2f\n",
            lr_rho, pchisq(lr_rho, 1, lower.tail = FALSE), wald_rho$z))
cat("\n--- Effects on E[y] per 1 SD (avg over LSOAs; 'zero' = inflation channel, in direct) ---\n")
for (e in eff[-1]) cat(sprintf("%-11s direct %7.4f (%.4f) indirect %7.4f (%.4f) total %7.4f (%.4f) zero %7.4f\n",
  e$name, e$avg_direct, e$se_direct, e$avg_indirect, e$se_indirect,
  e$avg_total, e$se_total, e$avg_zero))
cat("elapsed", round(timing$effects, 1), "s\n")
cat("\n--- Observed vs expected frequencies ---\n"); print(round(freq_tab, 1))
cat("\n--- Structural zeros ---\n")
cat("mean pi_hat:", round(mean(pi_hat), 4), " sum pi_hat:", round(sum(pi_hat), 1),
    " expected structural zeros among observed zeros:", round(sum(post_struct), 1),
    " zeros with P(structural) > 0.5:", sum(post_struct > 0.5), "\n")
cat("\n--- Moran's I (queen) ---\n")
for (m in names(moran)) cat(sprintf("%-6s I = %.4f  (p = %.3g)\n", m,
                                    moran[[m]]$estimate[1], moran[[m]]$p.value))
cat("\n--- Sensitivity to W ---\n")
st <- sapply(c(list(queen = fit_zisar), sens), function(f)
  c(rho = f$params$rho, f$params$beta[-1], kappa = f$params$kappa,
    sigma = f$params$sigma, f$params$gamma[-1], logLik = f$loglik, AIC = f$AIC))
rownames(st) <- c("rho", colnames(X)[-1], "kappa", "sigma", colnames(Z)[-1], "logLik", "AIC")
print(round(st, 3))
cat("\n--- Robustness: cumulative 2010-2019 count ---\n")
cat("zeros =", sum(y_cum == 0), " mean =", round(mean(y_cum), 2), " var =", round(var(y_cum), 1), "\n")
cat("AIC: NB", round(AIC(rob$nb), 1), " ZINB", round(AIC(rob$zinb), 1),
    " ZISAR", round(rob$zisar$AIC, 1), " (convergence", rob$zisar$convergence, ")\n")
cat("expected zeros: NB", round(sum(dnbinom(0, size = rob$nb$theta, mu = fitted(rob$nb))), 1),
    " ZINB", round(sum(predict(rob$zinb, type = "prob")[, 1]), 1),
    " ZISAR", round(zisar_marginal_freq(rob$zisar, 0), 1), "\n")
print(rob$zisar$coef_table); print(rob$zisar$transformed)
print(round(summary(rob$zinb)$coefficients$zero, 3))
cat("\nRestarts from the stopping point (the second with kappa fixed):\n")
for (r in fits$cum_restarts) cat(sprintf(
  "  convergence %d  logLik %.1f  rho %.3f  kappa %.1f  sigma %.3f  gamma %s\n",
  r$convergence, r$loglik, r$params$rho, r$params$kappa, r$params$sigma,
  paste(sprintf("%.3f", r$params$gamma), collapse = " ")))
cat("\n--- Robustness: post-reform outcome by capacity, and without the pre-reform barriers ---\n")
for (nm in c("all", "kw4", "nobar")) {
  f <- fits[[nm]]
  cat(sprintf("%-6s installations %d  zeros %d  logLik %.2f  k %d  AIC %.1f  BIC %.1f  rho %.4f  sigma %.3f  kappa %.1f  convergence %d  max_grad %.2g\n",
              nm, sum(f$y), sum(f$y == 0), f$loglik, length(f$coefficients), f$AIC, f$BIC,
              f$params$rho, f$params$sigma, f$params$kappa, f$convergence, f$max_grad))
}
lr_bar <- 2 * (fit_zisar$loglik - fits$nobar$loglik)
cat(sprintf("LR pre-reform barriers (%d coefficients): %.1f (p = %.3g); AIC %+.1f, BIC %+.1f\n",
            2 * ncol(cov_b), lr_bar, pchisq(lr_bar, 2 * ncol(cov_b), lower.tail = FALSE),
            fit_zisar$AIC - fits$nobar$AIC, fit_zisar$BIC - fits$nobar$BIC))
cat("Coefficients, all capacities:\n"); print(round(fits$all$coef_table, 3))
cat("Coefficients, up to 4 kW:\n"); print(round(fits$kw4$coef_table, 3))
cat("Coefficients, without the pre-reform barriers:\n"); print(round(fits$nobar$coef_table, 3))
cat("\n--- Score-type test of spatial dependence in the structural regime (Moran's I of r, bootstrap) ---\n")
for (nm in names(struct)) for (k in names(struct[[nm]]$obs))
  cat(sprintf("%-6s %-8s observed %.4f  null mean %.4f  sd %.4f  z %.2f  p %.3f  (B = %d)\n", nm, k,
              struct[[nm]]$obs[k], struct[[nm]]$null_mean[k], struct[[nm]]$null_sd[k],
              struct[[nm]]$z[k], struct[[nm]]$p[k], struct[[nm]]$B))
cat("\n--- Sensitivity: general ZISAR (SAR field in logit(pi)) ---\n")
cat("Laplace log-likelihood gain over the restricted model, other parameters at the restricted estimates:\n")
print(round(xtabs(gain ~ sigma0 + rho0, gen$grid), 2))
cat(sprintf("\nOptimiser from rho0 = 0.9, sigma0 = 0.8: %s after %d iterations (%.0f s), logLik %.2f (gain %.2f)\n",
            gen$opt$message, gen$opt$iterations, gen$elapsed, -gen$opt$objective,
            -gen$opt$objective - fit_zisar$loglik))
cat(sprintf("rho %.4f  rho0 %.4f  sigma0 %.4f\n", tanh(gen$par["rho_raw"]),
            tanh(gen$par["rho0_raw"]), exp(gen$par["log_sigma0"])))
cat(sprintf("largest gradient component: %s, AD %.1f, central difference %.1f\n",
            gen$grad_which, gen$grad_ad, gen$grad_fd))
cat("objective along log(sigma0) around the stopping point:\n")
cat(sprintf("  %+.4f %.3f\n", gen$scan_t, gen$scan_f), sep = "")
cat("\n--- Sensitivity: spatial inflation equation at the LSOA level ---\n")
for (nm in c("zisar", "slx", "sar")) {
  f <- list(zisar = fit_zisar, slx = fit_slx, sar = fit_sar)[[nm]]
  cat(sprintf("%-6s logLik %.2f  k %d  AIC %.1f  BIC %.1f  rho %.4f  convergence %d  max_grad %.2g\n",
              nm, f$loglik, length(f$coefficients), f$AIC, f$BIC, f$params$rho,
              f$convergence, f$max_grad))
}
cat(sprintf("LR WZ in Z: %.2f on %d df (p = %.3g)%s\n", 2 * (fit_slx$loglik - fit_zisar$loglik),
            ncol(WZ), pchisq(2 * (fit_slx$loglik - fit_zisar$loglik), ncol(WZ), lower.tail = FALSE),
            if (fit_slx$convergence != 0) "  NOT VALID: the WZ fit did not converge" else ""))
cat(sprintf("SAR filter: rho0 %.4f (atanh %.4f, SE %.4f); LR %.2f on 1 df (p = %.3g)\n",
            fit_sar$params$rho0, fit_sar$coef_table["rho0_raw", "estimate"],
            fit_sar$coef_table["rho0_raw", "se"], 2 * (fit_sar$loglik - fit_zisar$loglik),
            pchisq(2 * (fit_sar$loglik - fit_zisar$loglik), 1, lower.tail = FALSE)))
cat("Inflation coefficients, WZ in Z:\n")
print(round(fit_slx$coef_table[grep("^gamma", rownames(fit_slx$coef_table)), ], 3))
cat("\nPosterior structural-zero probability by pre-reform installations, AUC\n")
print(round(validation, 3))
cat("\nZero clustering: observed vs draws from each fit (mean [2.5%, 97.5%])\n")
for (v in c("zeros", "zeros_lon", "moran", "moran_lon")) {
  cat(sprintf("%-10s observed %9.3f\n", v, ppc$observed[[v]]))
  for (m in c("zisar", "slx", "sar")) cat(sprintf("  %-6s %9.3f [%9.3f, %9.3f]\n", m,
    mean(ppc[[m]][[v]]), quantile(ppc[[m]][[v]], .025), quantile(ppc[[m]][[v]], .975)))
}
zg <- as.numeric(Z %*% fit_zisar$params$gamma)
cat(sprintf("cor(Z gamma, W Z gamma) %.3f; London zeros %d, expected structural among them %.1f\n",
            cor(zg, as.numeric(W %*% zg)), sum(y[london] == 0), sum(post_struct[london])))
sink()

# -----------------------------------------------------------------------------
# 8. Maps: adoption rate and structural-zero probability (England, London)
# -----------------------------------------------------------------------------
map <- readRDS("../data/solar/lsoa_map.rds")
idx <- match(map$unit, solar$unit)
map$rate <- 1000 * y[idx] / solar$households[idx]
map$pi   <- pi_hat[idx]
map$zero <- y[idx] == 0
london <- substr(solar$lad19cd[idx], 1, 3) == "E09"
theme_map <- theme_void(base_size = 8) +
  theme(legend.position = "bottom", legend.key.width = unit(0.6, "cm"),
        legend.key.height = unit(0.22, "cm"), legend.title.position = "top",
        legend.title = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5))
g1 <- ggplot(map) + geom_sf(aes(fill = pmin(rate, 10)), colour = NA) +
  scale_fill_viridis_c(option = "magma", name = "installations per 1,000 households",
                       breaks = c(0, 5, 10), labels = c("0", "5", "10+")) +
  ggtitle("(a) England") + theme_map
g2 <- ggplot(map[london, ]) + geom_sf(aes(fill = zero), colour = NA) +
  scale_fill_manual(values = c(`FALSE` = "grey85", `TRUE` = "#b2182b"),
                    labels = c("at least one", "none"), name = "installations") +
  ggtitle("(b) London: observed zeros") + theme_map
g3 <- ggplot(map[london, ]) + geom_sf(aes(fill = pi), colour = NA) +
  scale_fill_viridis_c(option = "viridis", trans = "sqrt", limits = c(0, 1),
                       breaks = c(0.01, 0.1, 0.3, 0.6, 1),
                       name = "structural-zero probability") +
  ggtitle(expression(paste("(c) London: ", hat(pi)[i]))) + theme_map
dir.create("../tex/figures", showWarnings = FALSE)
ggsave("../tex/figures/solar_maps.png", g1 + g2 + g3, width = 7.2, height = 3.6, dpi = 300)
