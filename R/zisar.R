#' Fit a zero-inflated spatial autoregressive (ZISAR) model
#'
#' Fits the ZISAR model of Salvini, Arbia and Nardelli (2026): a
#' zero-inflated negative binomial model whose log-mean follows a
#' simultaneous autoregressive (SAR) process with unit-level random effects.
#' By default the inflation probability depends on covariates only (the
#' restricted model); with `zero_spatial = TRUE` it has a spatial random
#' field of its own (the general model). `zisar()` takes a two-part formula
#' and a data frame; `zisar_fit()` is the matrix interface underneath it.
#'
#' @section Model:
#' For units \eqn{i = 1, \dots, n} linked by a spatial weights matrix
#' \eqn{W},
#' \deqn{\eta = \rho W \eta + X\beta + u, \qquad u \sim N(0, \sigma^2 I),}
#' so that \eqn{\eta = (I - \rho W)^{-1}(X\beta + u)}. Each count is zero with
#' probability \eqn{\pi_i} (a structural zero) and otherwise negative
#' binomial,
#' \deqn{y_i \sim \pi_i \, 1\{y_i = 0\} + (1 - \pi_i)\,
#'   \mathrm{NB}(\mu_i, \kappa), \qquad \mu_i = \exp(\eta_i + o_i),
#'   \qquad \mathrm{logit}(\pi_i) = z_i'\gamma,}
#' with \eqn{\mathrm{Var}(y_i \mid \text{not a structural zero}) = \mu_i +
#' \mu_i^2/\kappa} and \eqn{o_i} an optional offset (log exposure). The
#' offset scales unit \eqn{i} only: it does not pass through the spatial
#' multiplier.
#'
#' In the general model (`zero_spatial = TRUE`) the inflation predictor has
#' a SAR random field of its own,
#' \deqn{\mathrm{logit}(\pi_i) = z_i'\gamma + g_i, \qquad
#'   g = \rho_0 W g + e_0, \qquad e_0 \sim N(0, \sigma_0^2 I),}
#' independent of \eqn{u}, so that structural zeros can cluster in space
#' beyond what the covariates explain. The restricted model is the limit
#' \eqn{\sigma_0 \to 0}.
#'
#' @section Estimation:
#' The latent field \eqn{\eta} is integrated out with the Laplace
#' approximation in 'TMB' (together with \eqn{g} in the general model), and
#' the marginal likelihood is maximised with [stats::nlminb()] over
#' \eqn{(\beta, \gamma, \rho, \kappa, \sigma)}, plus \eqn{(\rho_0,
#' \sigma_0)} in the general model, with \eqn{\rho = \tanh(\rho_{raw})}
#' restricted to \eqn{(-1, 1)}. Standard errors come from the inverse Hessian
#' of the Laplace objective; those of \eqn{\rho}, \eqn{\kappa},
#' \eqn{\sigma}, \eqn{\rho_0} and \eqn{\sigma_0} from the delta method.
#'
#' The Laplace approximation can break down. For a zero count, the
#' conditional distribution of \eqn{\eta_i} can be bimodal (a structural zero
#' with a large mean, or a negative binomial zero with a small one), and
#' when \eqn{\sigma} is large the approximation produces spurious optima with
#' a positive log-likelihood and a huge gradient. This happens mostly with
#' small samples and misspecified models (an omitted covariate inflates
#' \eqn{\sigma}). Each fit is therefore checked: when the optimiser does not
#' converge, or the optimum has a gradient above 0.01 or a positive
#' log-likelihood, the fit is restarted from three other values of
#' \eqn{\kappa} and \eqn{\sigma} and from up to three random perturbations of
#' the best point (this uses the random number generator). If no usable
#' optimum is found, `zisar()` warns and `laplace_ok` is `FALSE`.
#'
#' The general model is more exposed to this breakdown: a zero count can
#' come from a large \eqn{g_i} (a structural zero) or from a small
#' \eqn{\eta_i} (a negative binomial zero), and the two latent fields
#' compete for the same zeros. Fit the restricted model first, and treat a
#' general fit with `laplace_ok = FALSE` as uninformative about
#' \eqn{(\rho_0, \sigma_0)}.
#'
#' The log-determinant \eqn{\log|I - \rho W|} is computed either from the
#' eigenvalues of \eqn{W}, found once in \eqn{O(n^3)} time and \eqn{O(n^2)}
#' memory (`logdet = "eigen"`), or from a sparse Cholesky factorisation of
#' \eqn{(I - \rho W)'(I - \rho W)} at each evaluation (`logdet = "sparse"`),
#' which scales to tens of thousands of units. `"auto"` uses the sparse
#' method when \eqn{n > 5000}.
#'
#' @section Nested models:
#' `fixed` holds parameters at given values, which gives the special cases
#' of the model: `fixed = list(rho = 0)` is the aspatial zero-inflated
#' negative binomial with random effects (ZINB-RE), and with an
#' intercept-only inflation part (`y ~ x | 1`), `fixed = list(gamma = -20)`
#' sets \eqn{\pi_i \approx 0}, the SAR negative binomial with random effects
#' (SAR-NB-RE). Likelihood-ratio tests between nested fits use
#' [logLik()][logLik.zisar]. The restricted model sits on the boundary of the
#' general one (\eqn{\sigma_0 = 0}), so the likelihood-ratio statistic
#' between them does not have the usual chi-squared distribution: its
#' p-value from \eqn{\chi^2_2} is conservative.
#'
#' @param formula A two-part formula `y ~ x1 + x2 | z1 + z2`: the terms
#'   before `|` enter the count process, those after it the zero-inflation
#'   process. Without `|`, both processes use the same terms. Use `| 1` for
#'   a constant inflation probability. An `offset()` term in the count part
#'   enters the log-mean with coefficient one.
#' @param data A data frame with one row per spatial unit, in the order of
#'   the rows of `W`. Missing values are not allowed: every unit enters the
#'   spatial multiplier, so none can be dropped.
#' @param W Spatial weights: an \eqn{n \times n} matrix (base or sparse
#'   'Matrix') or an `spdep` `listw` object, with zero diagonal and
#'   preferably row-standardised.
#' @param zero_spatial If `TRUE`, fit the general model, with a SAR random
#'   field in the inflation predictor (see Model).
#' @param fixed Optional named list of parameters held at fixed values:
#'   `rho`, `kappa`, `sigma`, and in the general model `rho0` and `sigma0`
#'   (scalars, natural scale), `beta` (length `ncol(X)`) or `gamma` (length
#'   `ncol(Z)`). In `beta` and `gamma`, `NA` entries are estimated.
#' @param start Optional named list of starting values, with the same names
#'   and scales as `fixed`. Defaults: `beta` from a Poisson regression,
#'   inflation intercept from the share of zeros, `rho = 0`, `kappa = 1`,
#'   `sigma = 0.5`, and `rho0 = 0`, `sigma0 = 0.5`.
#' @param logdet Method for the log-determinant: `"auto"`, `"eigen"` or
#'   `"sparse"`. See Estimation.
#' @param control List of control parameters for [stats::nlminb()],
#'   merged with the defaults `eval.max = 1000`, `iter.max = 500`,
#'   `rel.tol = 1e-10`.
#' @param silent If `FALSE`, print the 'TMB' tracing output.
#'
#' @return An object of class `"zisar"`, a list with, among others:
#'   \item{coefficients}{List with the `count` (\eqn{\beta}) and `zero`
#'     (\eqn{\gamma}) coefficients.}
#'   \item{rho, kappa, sigma}{Spatial autoregressive parameter, negative
#'     binomial shape and random-effect standard deviation.}
#'   \item{rho0, sigma0}{In the general model, the spatial autoregressive
#'     parameter and standard deviation of the inflation field; `NULL`
#'     otherwise.}
#'   \item{eta}{Mode of the latent field \eqn{\eta} given the data.}
#'   \item{g}{Mode of the inflation field \eqn{g} (zero in the restricted
#'     model).}
#'   \item{fitted.values}{Fitted means \eqn{E[y_i] = (1 - \pi_i)\mu_i}.}
#'   \item{loglik}{Laplace-approximated marginal log-likelihood.}
#'   \item{convergence, message}{Convergence code and message of
#'     [stats::nlminb()].}
#'   \item{max_grad, laplace_ok}{Largest gradient component at the optimum,
#'     and whether the optimum passed the checks described in Estimation.}
#'   \item{tmb}{The 'TMB' object, for advanced use.}
#'
#'   Use the methods [summary()][summary.zisar], [coef()][coef.zisar],
#'   [vcov()][vcov.zisar], [logLik()][logLik.zisar],
#'   [predict()][predict.zisar], [residuals()][residuals.zisar] and
#'   [zisar_effects()] rather than the components directly.
#'
#' @references
#' Salvini, N., Arbia, G. and Nardelli, V. (2026). ZISAR: A zero-inflated
#' spatial autoregressive model for count data with random effects.
#' Working paper, Università Cattolica del Sacro Cuore.
#'
#' Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H. and Bell, B. M.
#' (2016). TMB: Automatic differentiation and Laplace approximation.
#' *Journal of Statistical Software*, 70(5), 1--21.
#' \doi{10.18637/jss.v070.i05}
#'
#' @seealso [zisar_effects()] for direct, indirect and total effects,
#'   [rzisar()] to simulate from the model.
#' @export
#' @examples
#' set.seed(1)
#' sim <- rzisar(12, beta = c(1, 0.5), rho = 0.4, gamma = c(-1, 0.8))
#' fit <- zisar(y ~ x1 | z1, data = sim$data, W = sim$W)
#' summary(fit)
#'
#' # Aspatial ZINB with random effects, and a likelihood-ratio test of rho = 0
#' fit0 <- zisar(y ~ x1 | z1, data = sim$data, W = sim$W,
#'               fixed = list(rho = 0))
#' lr <- 2 * (logLik(fit) - logLik(fit0))
#' pchisq(as.numeric(lr), df = 1, lower.tail = FALSE)
#'
#' \donttest{
#' # General model: zeros clustered by a spatial field of their own
#' set.seed(2)
#' sim2 <- rzisar(15, beta = c(1, 0.5), rho = 0.4, gamma = c(-1, 0.8),
#'                rho0 = 0.6, sigma0 = 1)
#' fit2 <- zisar(y ~ x1 | z1, data = sim2$data, W = sim2$W,
#'               zero_spatial = TRUE)
#' summary(fit2)
#' }
zisar <- function(formula, data, W, zero_spatial = FALSE, fixed = NULL,
                  start = NULL, logdet = c("auto", "eigen", "sparse"),
                  control = list(), silent = TRUE) {
  cl <- match.call()
  if (missing(data)) data <- environment(formula)
  f <- .split_formula(formula)

  mf_count <- model.frame(f$count, data = data, na.action = stats::na.pass,
                          drop.unused.levels = TRUE)
  mf_zero <- model.frame(f$zero, data = data, na.action = stats::na.pass,
                         drop.unused.levels = TRUE)
  if (anyNA(mf_count) || anyNA(mf_zero))
    stop("Missing values in the model variables. Every unit in `W` enters ",
         "the spatial multiplier, so units cannot be dropped: impute the ",
         "missing values, or remove the units and rebuild `W`.", call. = FALSE)
  if (f$bar && !is.null(attr(terms(mf_zero), "offset")))
    warning("`offset()` in the zero-inflation part is ignored.", call. = FALSE)

  mt_count <- terms(mf_count)
  mt_zero <- terms(mf_zero)
  y <- model.response(mf_count)
  X <- model.matrix(mt_count, mf_count)
  Z <- model.matrix(mt_zero, mf_zero)
  offset <- model.offset(mf_count)

  fit <- zisar_fit(y, X, Z, W, offset = offset, zero_spatial = zero_spatial,
                   fixed = fixed, start = start, logdet = logdet,
                   control = control, silent = silent)
  fit$call <- cl
  fit$formula <- formula
  fit$terms <- list(count = mt_count, zero = mt_zero)
  fit$xlevels <- list(count = stats::.getXlevels(mt_count, mf_count),
                      zero = stats::.getXlevels(mt_zero, mf_zero))
  fit$contrasts <- list(count = attr(X, "contrasts"),
                        zero = attr(Z, "contrasts"))
  fit
}

#' @rdname zisar
#' @param y Vector of non-negative integer counts.
#' @param X Design matrix of the count process, \eqn{n \times p}; include a
#'   column of ones for an intercept.
#' @param Z Design matrix of the zero-inflation process, \eqn{n \times q}.
#'   `NULL` gives an intercept only.
#' @param offset Optional vector of length \eqn{n} added to the count
#'   log-mean (for example the log of an exposure).
#' @export
zisar_fit <- function(y, X, Z = NULL, W, offset = NULL, zero_spatial = FALSE,
                      fixed = NULL, start = NULL,
                      logdet = c("auto", "eigen", "sparse"),
                      control = list(), silent = TRUE) {
  n <- length(y)
  if (!is.numeric(y) || anyNA(y) || any(y < 0) || any(y != round(y)))
    stop("`y` must be a vector of non-negative integer counts without ",
         "missing values.", call. = FALSE)
  X <- .as_design(X, n, "X")
  if (is.null(Z)) Z <- matrix(1, n, 1L, dimnames = list(NULL, "(Intercept)"))
  Z <- .as_design(Z, n, "Z")
  p <- ncol(X)
  q <- ncol(Z)
  if (is.null(offset)) offset <- rep(0, n)
  if (length(offset) != n || any(!is.finite(offset)))
    stop("`offset` must be a finite vector with one value per observation.",
         call. = FALSE)
  W <- .as_weights(W, n)
  if (!is.logical(zero_spatial) || length(zero_spatial) != 1L ||
      is.na(zero_spatial))
    stop("`zero_spatial` must be TRUE or FALSE.", call. = FALSE)
  zs <- zero_spatial
  logdet <- match.arg(logdet)
  if (logdet == "auto") logdet <- if (n > 5000) "sparse" else "eigen"

  if (logdet == "eigen") {
    eig <- eigen(as.matrix(W), only.values = TRUE)$values
    eig_re <- Re(eig)
    eig_im <- Im(eig)
    I_n <- WpWt <- WtW <- as(Matrix::Matrix(0, 1, 1, sparse = TRUE),
                             "TsparseMatrix")
  } else {
    eig_re <- eig_im <- numeric(0)
    I_n <- as(Matrix::Diagonal(n), "TsparseMatrix")
    WpWt <- as(W + Matrix::t(W), "TsparseMatrix")
    WtW <- as(Matrix::crossprod(W), "TsparseMatrix")
  }

  data <- list(
    y = as.numeric(y),
    X = X,
    Z = Z,
    offset = as.numeric(offset),
    W = as(W, "TsparseMatrix"),
    eigenvalues_W_re = eig_re,
    eigenvalues_W_im = eig_im,
    logdet_sparse = as.integer(logdet == "sparse"),
    I_n = I_n,
    WpWt = WpWt,
    WtW = WtW,
    zero_spatial = as.integer(zs)
  )

  # Starting values: Poisson GLM for beta, share of zeros for the inflation
  # intercept
  init_beta <- tryCatch(
    suppressWarnings(stats::glm.fit(X, y, family = stats::poisson(),
                                    offset = offset)$coefficients),
    error = function(e) rep(0, p))
  init_beta[is.na(init_beta)] <- 0
  init_gamma <- rep(0, q)
  int_z <- .intercept_col(Z)
  if (!is.na(int_z))
    init_gamma[int_z] <- qlogis(max(0.01, min(0.99, mean(y == 0))))

  parameters <- list(
    beta       = unname(init_beta),
    gamma      = init_gamma,
    rho_raw    = 0,
    log_kappa  = 0,
    log_sigma  = log(0.5),
    rho0_raw   = 0,
    log_sigma0 = log(0.5),
    eta        = as.numeric(X %*% init_beta),
    g          = if (zs) rep(0, n) else numeric(0)
  )
  st <- .internal_fixed(start, p, q, zs, "start")
  for (b in names(st)) {
    v <- st[[b]]
    parameters[[b]][!is.na(v)] <- v[!is.na(v)]
  }
  if (!is.null(st$beta)) parameters$eta <- as.numeric(X %*% parameters$beta)
  fx <- .internal_fixed(fixed, p, q, zs)
  map <- list()
  for (b in names(fx)) {
    v <- fx[[b]]
    parameters[[b]][!is.na(v)] <- v[!is.na(v)]
    map[[b]] <- factor(ifelse(is.na(v), seq_along(v), NA))
  }
  # The inflation field exists only in the general model
  if (!zs) map[c("rho0_raw", "log_sigma0")] <- list(factor(NA), factor(NA))

  # Bounds: beta and gamma free, rho_raw (rho0_raw) in [-4, 4]
  # (|rho| < 0.9993), kappa in [0.01, 500], sigma (sigma0) in [0.001, 10]
  free <- .free_mask(fx, p, q, zs)
  lower <- c(rep(-Inf, p + q), -4, log(0.01), log(0.001),
             if (zs) c(-4, log(0.001)))[free]
  upper <- c(rep(Inf, p + q), 4, log(500), log(10),
             if (zs) c(4, log(10)))[free]
  if (!any(free)) stop("`fixed` leaves no parameter to estimate.", call. = FALSE)

  ctrl <- utils::modifyList(list(eval.max = 1000, iter.max = 500,
                                 rel.tol = 1e-10), control)
  make_obj <- function(par0)
    TMB::MakeADFun(data = data, parameters = par0,
                   random = c("eta", if (zs) "g"), map = map,
                   DLL = "zisar", silent = silent)
  best <- .optimise(make_obj(parameters), NULL, lower, upper, ctrl)

  # A failed fit is restarted from other values of kappa and sigma, then from
  # perturbations of the best point found (as in the code of the paper)
  if (!.fit_ok(best) || best$opt$convergence != 0) {
    alts <- list(c(log_kappa = 0, log_sigma = log(0.1)),
                 c(log_kappa = log(5), log_sigma = log(0.2)),
                 c(log_kappa = log(2), log_sigma = log(1)))
    for (a in alts) {
      par_a <- parameters
      for (b in setdiff(names(a), names(fx))) par_a[[b]] <- a[[b]]
      best <- .better_fit(best, .optimise(make_obj(par_a), NULL, lower,
                                          upper, ctrl))
    }
    for (attempt in 1:3) {
      if (is.null(best) || (.fit_ok(best) && best$opt$convergence == 0)) break
      x0 <- pmin(pmax(best$opt$par + rnorm(length(best$opt$par), 0, 0.1),
                      lower), upper)
      best <- .better_fit(best, .optimise(best$obj, x0, lower, upper, ctrl))
    }
  }
  if (is.null(best))
    stop("The optimiser failed from every starting point. Check the model ",
         "specification, or supply starting values with `start`.",
         call. = FALSE)
  obj <- best$obj
  opt <- best$opt
  if (!.fit_ok(best)) {
    warning("The Laplace approximation is unreliable at the optimum ",
            "(gradient ", signif(opt$max_grad, 3), ", log-likelihood ",
            signif(-opt$objective, 5), "). With many zeros and a large ",
            "random-effect variance the conditional distribution of eta can ",
            "be bimodal for zero counts. Check the specification, simplify ",
            "the model, or try other starting values with `start`.",
            call. = FALSE)
  } else if (opt$convergence != 0) {
    warning("The optimiser did not converge: ", opt$message, call. = FALSE)
  }

  # Hessian of the Laplace objective, shared with sdreport
  H <- tryCatch(suppressWarnings(optimHess(opt$par, obj$fn, obj$gr)),
                error = function(e) NULL)
  sdr <- tryCatch(suppressWarnings(
    TMB::sdreport(obj, par.fixed = opt$par, hessian.fixed = H)),
    error = function(e) NULL)

  par_names <- c(paste0("count_", colnames(X)), paste0("zero_", colnames(Z)),
                 "rho_raw", "log_kappa", "log_sigma",
                 if (zs) c("rho0_raw", "log_sigma0"))
  vc <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(vc) || any(!is.finite(vc)) || any(diag(vc) <= 0)) {
    warning("The Hessian is not positive definite: standard errors are ",
            "not available. The model may be over-parameterised, or the ",
            "optimum on a boundary.", call. = FALSE)
    vc <- matrix(NA_real_, sum(free), sum(free))
  }
  dimnames(vc) <- list(par_names[free], par_names[free])

  pl <- obj$env$parList(opt$par)
  par <- stats::setNames(c(pl$beta, pl$gamma, pl$rho_raw, pl$log_kappa,
                           pl$log_sigma,
                           if (zs) c(pl$rho0_raw, pl$log_sigma0)), par_names)
  beta <- stats::setNames(pl$beta, colnames(X))
  gamma <- stats::setNames(pl$gamma, colnames(Z))
  obj$fn(opt$par)
  lp <- obj$env$last.par
  eta <- as.numeric(lp[names(lp) == "eta"])
  g <- if (zs) as.numeric(lp[names(lp) == "g"]) else rep(0, n)

  # Delta-method SEs of the natural-scale parameters
  nat_names <- c("rho", "kappa", "sigma", if (zs) c("rho0", "sigma0"))
  rep_se <- stats::setNames(rep(NA_real_, length(nat_names)), nat_names)
  if (!is.null(sdr) && all(is.finite(vc))) {
    s <- summary(sdr, "report")
    rep_se[rownames(s)] <- s[, "Std. Error"]
  }
  rep_se[!free[p + q + seq_along(nat_names)]] <- NA_real_

  mu <- exp(eta + offset)
  pi <- as.numeric(plogis(as.numeric(Z %*% gamma) + g))

  structure(
    list(
      coefficients  = list(count = beta, zero = gamma),
      rho           = tanh(pl$rho_raw),
      kappa         = exp(pl$log_kappa),
      sigma         = exp(pl$log_sigma),
      zero_spatial  = zs,
      rho0          = if (zs) tanh(pl$rho0_raw),
      sigma0        = if (zs) exp(pl$log_sigma0),
      se_natural    = rep_se,
      par           = par,
      free          = stats::setNames(free, par_names),
      vcov_internal = vc,
      eta           = eta,
      g             = g,
      mu            = mu,
      pi            = pi,
      fitted.values = (1 - pi) * mu,
      loglik        = -opt$objective,
      df            = sum(free),
      n             = n,
      y             = as.numeric(y),
      X             = X,
      Z             = Z,
      offset        = as.numeric(offset),
      W             = W,
      logdet        = logdet,
      convergence   = opt$convergence,
      message       = opt$message,
      max_grad      = opt$max_grad,
      laplace_ok    = .fit_ok(best),
      iterations    = opt$iterations,
      tmb           = obj,
      sdreport      = sdr
    ),
    class = "zisar"
  )
}

# One nlminb run from x0 (default: the object's starting values). Records the
# largest gradient component away from the bounds; NULL if the run fails.
.optimise <- function(obj, x0, lower, upper, ctrl) {
  if (is.null(x0)) x0 <- obj$par
  # NA/NaN evaluations along the way are expected; the result is checked below
  opt <- tryCatch(suppressWarnings(nlminb(x0, obj$fn, obj$gr, lower = lower,
                                          upper = upper, control = ctrl)),
                  error = function(e) NULL)
  if (is.null(opt) || !is.finite(opt$objective)) return(NULL)
  g <- tryCatch(as.numeric(obj$gr(opt$par)), error = function(e) NA_real_)
  inside <- opt$par > lower + 1e-6 & opt$par < upper - 1e-6
  opt$max_grad <- if (anyNA(g)) Inf else max(abs(g[inside]), 0)
  list(obj = obj, opt = opt)
}

# A usable optimum: small gradient and a log-likelihood below zero (the
# Laplace approximation of a probability above one signals a breakdown)
.fit_ok <- function(fit) {
  !is.null(fit) && fit$opt$max_grad < 0.01 && fit$opt$objective > 0
}

.better_fit <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  if (.fit_ok(a) != .fit_ok(b)) return(if (.fit_ok(a)) a else b)
  if (b$opt$objective < a$opt$objective) b else a
}

# Split `y ~ count | zero` into its two parts; without `|`, the zero part
# reuses the count terms
.split_formula <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("`formula` must be a two-sided formula such as `y ~ x1 + x2 | z1`.",
         call. = FALSE)
  rhs <- formula[[3L]]
  env <- environment(formula)
  bar <- is.call(rhs) && identical(rhs[[1L]], as.name("|"))
  count_rhs <- if (bar) rhs[[2L]] else rhs
  zero_rhs <- if (bar) rhs[[3L]] else rhs
  list(count = stats::as.formula(call("~", formula[[2L]], count_rhs), env = env),
       zero = stats::as.formula(call("~", zero_rhs), env = env),
       bar = bar)
}

# Numeric design matrix with n rows and column names; an unnamed column of
# ones is named "(Intercept)"
.as_design <- function(M, n, what) {
  M <- as.matrix(M)
  if (!is.numeric(M) || nrow(M) != n || any(!is.finite(M)))
    stop(sprintf("`%s` must be a finite numeric matrix with %d rows.", what, n),
         call. = FALSE)
  cn <- colnames(M)
  if (is.null(cn)) cn <- rep("", ncol(M))
  empty <- is.na(cn) | cn == ""
  cn[empty] <- paste0(what, seq_len(ncol(M)))[empty]
  k <- .intercept_col(M)
  if (!is.na(k) && empty[k]) cn[k] <- "(Intercept)"
  colnames(M) <- cn
  if (anyDuplicated(colnames(M)))
    stop(sprintf("Column names of `%s` must be unique.", what), call. = FALSE)
  M
}

.intercept_col <- function(M) {
  k <- which(apply(M, 2, function(v) all(v == 1)))
  if (length(k)) k[1] else NA_integer_
}

# User-facing `fixed` (natural scale) -> internal parameter blocks
.internal_fixed <- function(fixed, p, q, zero_spatial = FALSE, arg = "fixed") {
  if (is.null(fixed) || !length(fixed)) return(list())
  if (!is.list(fixed) || is.null(names(fixed)) || any(names(fixed) == ""))
    stop(sprintf("`%s` must be a named list.", arg), call. = FALSE)
  bad <- setdiff(names(fixed), c("beta", "gamma", "rho", "kappa", "sigma",
                                 "rho0", "sigma0"))
  if (length(bad))
    stop(sprintf("Unknown parameter(s) in `%s`: ", arg),
         paste(bad, collapse = ", "),
         ". Use beta, gamma, rho, kappa, sigma, rho0 or sigma0.",
         call. = FALSE)
  if (!zero_spatial && any(c("rho0", "sigma0") %in% names(fixed)))
    stop(sprintf("`%s$rho0` and `%s$sigma0` need the general model ", arg,
                 arg), "(`zero_spatial = TRUE`).", call. = FALSE)
  out <- list()
  for (b in c("beta", "gamma")) {
    if (is.null(fixed[[b]])) next
    len <- if (b == "beta") p else q
    if (length(fixed[[b]]) != len)
      stop(sprintf("`%s$%s` must have length %d.", arg, b, len), call. = FALSE)
    out[[b]] <- as.numeric(fixed[[b]])
  }
  scal <- function(b, ok, tr, msg) {
    v <- fixed[[b]]
    if (length(v) != 1L || !is.finite(v) || !ok(v))
      stop(sprintf("`%s$%s` must be %s.", arg, b, msg), call. = FALSE)
    tr(v)
  }
  if (!is.null(fixed[["rho"]]))
    out$rho_raw <- scal("rho", function(v) abs(v) < 1, atanh,
                        "a number in (-1, 1)")
  if (!is.null(fixed[["kappa"]]))
    out$log_kappa <- scal("kappa", function(v) v > 0, log, "a positive number")
  if (!is.null(fixed[["sigma"]]))
    out$log_sigma <- scal("sigma", function(v) v > 0, log, "a positive number")
  if (!is.null(fixed[["rho0"]]))
    out$rho0_raw <- scal("rho0", function(v) abs(v) < 1, atanh,
                         "a number in (-1, 1)")
  if (!is.null(fixed[["sigma0"]]))
    out$log_sigma0 <- scal("sigma0", function(v) v > 0, log,
                           "a positive number")
  out
}

# Which of the p + q + 3 (+ 2 in the general model) fixed-effect parameters
# are estimated
.free_mask <- function(fx, p, q, zero_spatial = FALSE) {
  free_block <- function(b, len)
    if (is.null(fx[[b]])) rep(TRUE, len) else is.na(fx[[b]])
  c(free_block("beta", p), free_block("gamma", q), free_block("rho_raw", 1),
    free_block("log_kappa", 1), free_block("log_sigma", 1),
    if (zero_spatial) c(free_block("rho0_raw", 1),
                        free_block("log_sigma0", 1)))
}
