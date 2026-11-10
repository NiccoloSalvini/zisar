#' Fit a zero-inflated spatial autoregressive (ZISAR) model
#'
#' Fits the ZISAR model of Salvini, Arbia and Nardelli (2026): a
#' zero-inflated negative binomial model whose log-mean follows a
#' simultaneous autoregressive (SAR) process with unit-level random effects.
#' By default the inflation probability depends on covariates only (the
#' restricted model); `zero_spatial` adds spatial structure to it (the
#' general model). `zisar()` takes a two-part formula and a data frame;
#' `zisar_fit()` is the matrix interface underneath it.
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
#'   \qquad \mathrm{logit}(\pi_i) = \zeta_i,}
#' with \eqn{\mathrm{Var}(y_i \mid \text{not a structural zero}) = \mu_i +
#' \mu_i^2/\kappa} and \eqn{o_i} an optional offset (log exposure). The
#' offset scales unit \eqn{i} only: it does not pass through the spatial
#' multiplier.
#'
#' The inflation predictor \eqn{\zeta} is set by `zero_spatial`:
#' \describe{
#'   \item{`"none"`}{The restricted model, \eqn{\zeta = Z\gamma}.}
#'   \item{`"reduced"`}{A reduced-rank SAR field, \eqn{\zeta = Z\gamma +
#'     B\alpha} with \eqn{\alpha = \rho_0 W_K \alpha + e_0}, \eqn{e_0 \sim
#'     N(0, \sigma_0^2 I_K)}: \eqn{\alpha} holds the values of the field at
#'     the \eqn{K} knots of a lattice, \eqn{W_K} is the contiguity of the
#'     knots and \eqn{B} a basis of local functions built by
#'     [zisar_basis()], so that structural zeros can cluster in space beyond
#'     what the covariates explain. \eqn{\rho_0 \ge 0}. This is the
#'     recommended general model.}
#'   \item{`"field"`}{A SAR field with one element per unit, \eqn{\zeta =
#'     Z\gamma + g} with \eqn{g = \rho_0 W g + e_0}, \eqn{e_0 \sim N(0,
#'     \sigma_0^2 I)}. Its Laplace approximation fails (see Estimation).}
#'   \item{`"sar"`}{A SAR filter on the predictor, \eqn{\zeta = (I - \rho_0
#'     W)^{-1} Z\gamma}: the covariates of a unit's neighbours shift its
#'     inflation probability. No random effects are added.}
#' }
#' The random fields are independent of \eqn{u}. The restricted model is the
#' limit \eqn{\sigma_0 \to 0} of the fields, and \eqn{\rho_0 = 0} of the
#' filter.
#'
#' @section Estimation:
#' The latent field \eqn{\eta} (and \eqn{g} or \eqn{\alpha}) is integrated out
#' with the Laplace approximation in 'TMB', and the marginal likelihood is
#' maximised with [stats::nlminb()] over \eqn{(\beta, \gamma, \rho, \kappa,
#' \sigma)} and, in the general model, \eqn{\rho_0} and \eqn{\sigma_0}, with
#' \eqn{\rho = \tanh(\rho_{raw})} restricted to \eqn{(-1, 1)}. Standard
#' errors come from the inverse Hessian of the Laplace objective; those of
#' \eqn{\rho}, \eqn{\kappa}, \eqn{\sigma}, \eqn{\rho_0} and \eqn{\sigma_0}
#' from the delta method.
#'
#' With a random field in the inflation predictor (`"reduced"` or
#' `"field"`), the restricted model is fitted first, the Laplace likelihood of
#' the general model is evaluated on a grid of \eqn{(\rho_0, \sigma_0)} with
#' the other parameters at the restricted estimates, and the optimiser starts
#' from the best grid point. The grid and the gain over the restricted fit are
#' returned in `field_profile`; a gain that never rises above zero means the
#' data hold no field.
#'
#' The unit-level field (`"field"`) defeats the Laplace approximation: each
#' \eqn{g_i} is informed by a single outcome that is zero or not, and the
#' error of the approximation accumulates over the \eqn{n} elements of the
#' field. On data without a field the approximate likelihood rises with
#' \eqn{\sigma_0}, and the likelihood-ratio statistic against the restricted
#' model can reach thousands. Use `"reduced"`, whose knot values are each
#' informed by many units, and check any general fit on data simulated
#' without a field (`simulate()` of the restricted fit).
#'
#' Each fit is checked: when the optimiser does not converge, or the optimum
#' has a gradient above 0.01 or a positive log-likelihood, the fit is
#' restarted from three other values of \eqn{\kappa} and \eqn{\sigma} and from
#' up to three random perturbations of the best point (this uses the random
#' number generator; `restarts = FALSE` skips it). If no clean optimum is
#' found, `zisar()` warns and `laplace_ok` is `FALSE`. A positive
#' log-likelihood is always a spurious optimum of the approximation. A large
#' gradient is not specific: it marks spurious optima, but it also appears at
#' genuine optima when structural zeros occur among units of high intensity,
#' since the conditional distribution of \eqn{\eta_i} of a zero unit then has
#' two modes and the objective is not smooth. Fit the model to data simulated
#' from the fit to tell the two apart, and use a parametric bootstrap
#' (`simulate()` and refits started from the estimates) for standard errors.
#'
#' The log-determinant \eqn{\log|I - \rho W|} is computed either from the
#' eigenvalues of \eqn{W}, found once in \eqn{O(n^3)} time and \eqn{O(n^2)}
#' memory (`logdet = "eigen"`), or from a sparse Cholesky factorisation of
#' \eqn{(I - \rho W)'(I - \rho W)} at each evaluation (`logdet = "sparse"`),
#' which scales to tens of thousands of units. `"auto"` uses the sparse
#' method when \eqn{n > 5000}. The knot field uses the eigenvalues of
#' \eqn{W_K}.
#'
#' @section Nested models:
#' `fixed` holds parameters at given values, which gives the special cases
#' of the model: `fixed = list(rho = 0)` is the aspatial zero-inflated
#' negative binomial with random effects (ZINB-RE), and with an
#' intercept-only inflation part (`y ~ x | 1`), `fixed = list(gamma = -20)`
#' sets \eqn{\pi_i \approx 0}, the SAR negative binomial with random effects
#' (SAR-NB-RE). Likelihood-ratio tests between nested fits use
#' [logLik()][logLik.zisar]. The restricted model sits on the boundary of the
#' models with a field (\eqn{\sigma_0 = 0}), so the likelihood-ratio
#' statistic between them does not have the usual chi-squared distribution:
#' its p-value from \eqn{\chi^2_2} is conservative.
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
#' @param zero_spatial Inflation predictor: `"none"` (restricted model),
#'   `"reduced"`, `"field"` or `"sar"`; see Model. `FALSE` and `TRUE` stand for
#'   `"none"` and `"field"`.
#' @param basis For `zero_spatial = "reduced"`, the basis from
#'   [zisar_basis()].
#' @param fixed Optional named list of parameters held at fixed values:
#'   `rho`, `kappa`, `sigma`, `rho0` and `sigma0` (scalars, natural scale),
#'   `beta` (length `ncol(X)`) or `gamma` (length `ncol(Z)`). In `beta` and
#'   `gamma`, `NA` entries are estimated.
#' @param start Optional named list of starting values, with the same names
#'   and scales as `fixed`. Defaults: `beta` from a Poisson regression,
#'   inflation intercept from the share of zeros, `rho = 0`, `kappa = 1`,
#'   `sigma = 0.5`, `rho0 = 0` and, for the fields, `rho0` and `sigma0` from
#'   the grid described in Estimation. To restart from a previous fit, use
#'   `start = list(beta = coef(fit, "count"), gamma = coef(fit, "zero"),
#'   rho = fit$rho, kappa = fit$kappa, sigma = fit$sigma)`.
#' @param logdet Method for the log-determinant: `"auto"`, `"eigen"` or
#'   `"sparse"`. See Estimation.
#' @param restarts If `TRUE`, restart a fit that fails the checks described in
#'   Estimation. Set it to `FALSE` when a single fit is slow and a large
#'   gradient is expected.
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
#'   \item{zero_spatial}{The inflation predictor, as in the argument.}
#'   \item{rho0, sigma0}{Parameters of the inflation field or filter; `NULL`
#'     when absent.}
#'   \item{eta}{Mode of the latent field \eqn{\eta} given the data.}
#'   \item{g}{Mode of the inflation field (\eqn{g}, or \eqn{B\alpha} for the
#'     reduced-rank field); zero without a field.}
#'   \item{zeta}{The inflation predictor, \eqn{\mathrm{logit}(\pi)}.}
#'   \item{fitted.values}{Fitted means \eqn{E[y_i] = (1 - \pi_i)\mu_i}.}
#'   \item{loglik}{Laplace-approximated marginal log-likelihood.}
#'   \item{field_profile}{For the fields: the grid of \eqn{(\rho_0,
#'     \sigma_0)} with the gain in log-likelihood over the restricted fit,
#'     and the restricted log-likelihood.}
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
#'   [zisar_basis()] for the reduced-rank field, [rzisar()] to simulate from
#'   the model.
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
#' # General model: zeros clustered by a spatial field of their own, fitted
#' # with a reduced-rank field on knots 4 cells apart
#' set.seed(2)
#' sim2 <- rzisar(20, beta = c(1, 0.5), rho = 0.4, gamma = c(-1, 0.8),
#'                zero_spatial = "field", rho0 = 0.9, sigma0 = 0.5)
#' b <- zisar_basis(sim2$coords, h = 4)
#' fit2 <- zisar(y ~ x1 | z1, data = sim2$data, W = sim2$W,
#'               zero_spatial = "reduced", basis = b)
#' summary(fit2)
#' }
zisar <- function(formula, data, W,
                  zero_spatial = c("none", "reduced", "field", "sar"),
                  basis = NULL, fixed = NULL, start = NULL,
                  logdet = c("auto", "eigen", "sparse"), restarts = TRUE,
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
                   basis = basis, fixed = fixed, start = start,
                   logdet = logdet, restarts = restarts, control = control,
                   silent = silent)
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
zisar_fit <- function(y, X, Z = NULL, W, offset = NULL,
                      zero_spatial = c("none", "reduced", "field", "sar"),
                      basis = NULL, fixed = NULL, start = NULL,
                      logdet = c("auto", "eigen", "sparse"), restarts = TRUE,
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
  zm <- .zero_spatial_arg(zero_spatial)
  if (zm == "reduced") {
    if (!inherits(basis, "zisar_basis"))
      stop("`zero_spatial = \"reduced\"` needs a `basis` from zisar_basis().",
           call. = FALSE)
    if (nrow(basis$B) != n)
      stop(sprintf("`basis` was built for %d units, the data have %d.",
                   nrow(basis$B), n), call. = FALSE)
  } else if (!is.null(basis)) {
    warning("`basis` is used only with `zero_spatial = \"reduced\"`.",
            call. = FALSE)
    basis <- NULL
  }
  field <- zm %in% c("field", "reduced")
  logdet <- match.arg(logdet)
  if (logdet == "auto") logdet <- if (n > 5000) "sparse" else "eigen"
  fx <- .internal_fixed(fixed, p, q, zm)
  st <- .internal_fixed(start, p, q, zm, "start")

  # A field starts from the restricted fit and a grid of (rho0, sigma0)
  profile <- NULL
  if (field && is.null(st$rho0_raw) && is.null(st$log_sigma0) &&
      is.null(fx$rho0_raw) && is.null(fx$log_sigma0)) {
    keep <- c("beta", "gamma", "rho", "kappa", "sigma")
    r_fit <- zisar_fit(y, X, Z, W, offset = offset, zero_spatial = "none",
                       fixed = fixed[intersect(names(fixed), keep)],
                       start = start[intersect(names(start), keep)],
                       logdet = logdet, restarts = restarts,
                       control = control, silent = silent)
    st <- utils::modifyList(st, list(
      beta = unname(r_fit$coefficients$count),
      gamma = unname(r_fit$coefficients$zero),
      rho_raw = atanh(r_fit$rho), log_kappa = log(min(r_fit$kappa, 499)),
      log_sigma = log(r_fit$sigma)))
    profile <- list(restricted_loglik = r_fit$loglik, eta = r_fit$eta)
  }

  data <- list(
    y = as.numeric(y),
    X = X,
    Z = Z,
    offset = as.numeric(offset),
    W = as(W, "TsparseMatrix"),
    logdet_sparse = as.integer(logdet == "sparse"),
    zero_spatial = match(zm, c("none", "field", "sar", "reduced")) - 1L
  )
  empty <- as(Matrix::Matrix(0, 1, 1, sparse = TRUE), "TsparseMatrix")
  if (logdet == "eigen") {
    eig <- eigen(as.matrix(W), only.values = TRUE)$values
    data$eigenvalues_W_re <- Re(eig)
    data$eigenvalues_W_im <- Im(eig)
    data$I_n <- data$WpWt <- data$WtW <- empty
  } else {
    data$eigenvalues_W_re <- data$eigenvalues_W_im <- numeric(0)
    data$I_n <- as(Matrix::Diagonal(n), "TsparseMatrix")
    data$WpWt <- as(W + Matrix::t(W), "TsparseMatrix")
    data$WtW <- as(Matrix::crossprod(W), "TsparseMatrix")
  }
  data$B <- if (zm == "reduced") basis$B else empty
  data$Wk <- if (zm == "reduced") basis$Wk else empty
  data$eig_k_re <- if (zm == "reduced") basis$eig_re else numeric(0)
  data$eig_k_im <- if (zm == "reduced") basis$eig_im else numeric(0)

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
    g          = if (zm == "field") rep(0, n) else numeric(0),
    alpha      = if (zm == "reduced") rep(0, basis$K) else numeric(0)
  )
  for (b in names(st)) {
    v <- st[[b]]
    parameters[[b]][!is.na(v)] <- v[!is.na(v)]
  }
  if (!is.null(st$beta)) parameters$eta <- as.numeric(X %*% parameters$beta)
  if (!is.null(profile)) parameters$eta <- profile$eta
  map <- list()
  for (b in names(fx)) {
    v <- fx[[b]]
    parameters[[b]][!is.na(v)] <- v[!is.na(v)]
    map[[b]] <- factor(ifelse(is.na(v), seq_along(v), NA))
  }
  # rho0 only with a field or filter, sigma0 only with a field
  if (zm == "none") map$rho0_raw <- factor(NA)
  if (!field) map$log_sigma0 <- factor(NA)

  # Bounds: beta and gamma free, rho_raw (rho0_raw) in [-4, 4]
  # (|rho| < 0.9993), kappa in [0.01, 500], sigma (sigma0) in [0.001, 10];
  # the knot field has rho0 >= 0
  zextra <- .zero_extra(zm)
  free <- .free_mask(fx, p, q, zm)
  lower <- c(rep(-Inf, p + q), -4, log(0.01), log(0.001),
             c(if (zm == "reduced") 0 else -4, log(0.001))[seq_along(zextra)])
  upper <- c(rep(Inf, p + q), 4, log(500), log(10),
             c(4, log(10))[seq_along(zextra)])
  lower <- lower[free]
  upper <- upper[free]
  if (!any(free)) stop("`fixed` leaves no parameter to estimate.", call. = FALSE)

  ctrl <- utils::modifyList(list(eval.max = 1000, iter.max = 500,
                                 rel.tol = 1e-10), control)
  make_obj <- function(par0)
    TMB::MakeADFun(data = data, parameters = par0,
                   random = c("eta", if (zm == "field") "g",
                              if (zm == "reduced") "alpha"),
                   map = map, DLL = "zisar", silent = silent)

  obj0 <- make_obj(parameters)
  x0 <- NULL
  if (!is.null(profile)) {
    grid <- expand.grid(rho0 = c(0, 0.5, 0.8, 0.9, 0.95),
                        sigma0 = c(0.1, 0.3, 0.6, 1, 2))
    at <- function(r0, s0) {
      x <- obj0$par
      x[names(x) == "rho0_raw"] <- atanh(r0)
      x[names(x) == "log_sigma0"] <- log(s0)
      x
    }
    nll <- vapply(seq_len(nrow(grid)), function(i)
      tryCatch(suppressWarnings(obj0$fn(at(grid$rho0[i], grid$sigma0[i]))),
               error = function(e) NA_real_), numeric(1))
    grid$gain <- -nll - profile$restricted_loglik
    if (any(is.finite(grid$gain))) {
      k <- which.max(replace(grid$gain, !is.finite(grid$gain), -Inf))
      x0 <- pmin(pmax(at(grid$rho0[k], grid$sigma0[k]), lower), upper)
    }
    profile <- list(grid = grid, restricted_loglik = profile$restricted_loglik)
  }
  best <- .optimise(obj0, x0, lower, upper, ctrl)

  # A failed fit is restarted from other values of kappa and sigma, then from
  # perturbations of the best point found (as in the code of the paper)
  if (restarts && (!.fit_ok(best) || best$opt$convergence != 0)) {
    alts <- list(c(log_kappa = 0, log_sigma = log(0.1)),
                 c(log_kappa = log(5), log_sigma = log(0.2)),
                 c(log_kappa = log(2), log_sigma = log(1)))
    for (a in alts) {
      par_a <- parameters
      if (!is.null(x0)) {
        par_a$rho0_raw <- x0[names(x0) == "rho0_raw"]
        par_a$log_sigma0 <- x0[names(x0) == "log_sigma0"]
      }
      for (b in setdiff(names(a), names(fx))) par_a[[b]] <- a[[b]]
      best <- .better_fit(best, .optimise(make_obj(par_a), NULL, lower,
                                          upper, ctrl))
    }
    for (attempt in 1:3) {
      if (is.null(best) || (.fit_ok(best) && best$opt$convergence == 0)) break
      xp <- pmin(pmax(best$opt$par + rnorm(length(best$opt$par), 0, 0.1),
                      lower), upper)
      best <- .better_fit(best, .optimise(best$obj, xp, lower, upper, ctrl))
    }
  }
  if (is.null(best))
    stop("The optimiser failed from every starting point. Check the model ",
         "specification, or supply starting values with `start`.",
         call. = FALSE)
  obj <- best$obj
  opt <- best$opt
  if (opt$objective <= 0) {
    warning("The optimum has a positive log-likelihood (",
            signif(-opt$objective, 5), "), a spurious optimum of the Laplace ",
            "approximation. Check the specification, simplify the model, or ",
            "try other starting values with `start`; see ?zisar.",
            call. = FALSE)
  } else if (!.fit_ok(best)) {
    warning("The optimum has a large gradient (", signif(opt$max_grad, 3),
            "). It may be a spurious optimum of the Laplace approximation, ",
            "or a genuine one with structural zeros among units of high ",
            "intensity. Refit to data simulated from this fit to tell them ",
            "apart, and use a parametric bootstrap for standard errors; see ",
            "?zisar.", call. = FALSE)
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
                 "rho_raw", "log_kappa", "log_sigma", zextra)
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
                           c(pl$rho0_raw, pl$log_sigma0)[seq_along(zextra)]),
                         par_names)
  beta <- stats::setNames(pl$beta, colnames(X))
  gamma <- stats::setNames(pl$gamma, colnames(Z))
  obj$fn(opt$par)
  lp <- obj$env$last.par
  eta <- as.numeric(lp[names(lp) == "eta"])
  alpha <- if (zm == "reduced") as.numeric(lp[names(lp) == "alpha"])
  g <- switch(zm,
              field = as.numeric(lp[names(lp) == "g"]),
              reduced = as.numeric(basis$B %*% alpha),
              rep(0, n))
  rho0 <- if (zm != "none") tanh(pl$rho0_raw)
  zeta <- as.numeric(Z %*% gamma) + g
  if (zm == "sar")
    zeta <- as.numeric(Matrix::solve(Matrix::Diagonal(n) - rho0 * W, zeta))

  # Delta-method SEs of the natural-scale parameters
  nat_names <- c("rho", "kappa", "sigma",
                 c("rho0", "sigma0")[seq_along(zextra)])
  rep_se <- stats::setNames(rep(NA_real_, length(nat_names)), nat_names)
  if (!is.null(sdr) && all(is.finite(vc))) {
    s <- summary(sdr, "report")
    rep_se[rownames(s)] <- s[, "Std. Error"]
  }
  rep_se[!free[p + q + seq_along(nat_names)]] <- NA_real_

  mu <- exp(eta + offset)
  pi <- plogis(zeta)

  structure(
    list(
      coefficients  = list(count = beta, zero = gamma),
      rho           = tanh(pl$rho_raw),
      kappa         = exp(pl$log_kappa),
      sigma         = exp(pl$log_sigma),
      zero_spatial  = zm,
      rho0          = rho0,
      sigma0        = if (field) exp(pl$log_sigma0),
      se_natural    = rep_se,
      par           = par,
      free          = stats::setNames(free, par_names),
      vcov_internal = vc,
      eta           = eta,
      g             = g,
      alpha         = alpha,
      zeta          = zeta,
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
      basis         = basis,
      field_profile = profile,
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

# `zero_spatial` argument -> "none", "reduced", "field" or "sar"
.zero_spatial_arg <- function(zero_spatial) {
  choices <- c("none", "reduced", "field", "sar")
  if (is.logical(zero_spatial) && length(zero_spatial) == 1L &&
      !is.na(zero_spatial))
    return(if (zero_spatial) "field" else "none")
  if (!is.character(zero_spatial) || !length(zero_spatial) ||
      !zero_spatial[1] %in% choices)
    stop("`zero_spatial` must be one of \"none\", \"reduced\", \"field\" or ",
         "\"sar\".", call. = FALSE)
  zero_spatial[1]
}

# Extra internal parameters of each inflation predictor
.zero_extra <- function(zm) {
  switch(zm, none = character(0), sar = "rho0_raw",
         c("rho0_raw", "log_sigma0"))
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

# User-facing `fixed` (natural scale) -> internal parameter blocks; `zm` is
# the inflation predictor ("none", "reduced", "field" or "sar")
.internal_fixed <- function(fixed, p, q, zm = "none", arg = "fixed") {
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
  if (zm == "none" && any(c("rho0", "sigma0") %in% names(fixed)))
    stop(sprintf("`%s$rho0` and `%s$sigma0` need a spatial inflation ", arg,
                 arg), "predictor (`zero_spatial`).", call. = FALSE)
  if (zm == "sar" && "sigma0" %in% names(fixed))
    stop(sprintf("`%s$sigma0` needs a random inflation field ", arg),
         "(`zero_spatial = \"reduced\"` or `\"field\"`).", call. = FALSE)
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
  if (!is.null(fixed[["rho0"]])) {
    if (zm == "reduced")
      out$rho0_raw <- scal("rho0", function(v) v >= 0 && v < 1, atanh,
                           "a number in [0, 1) for the reduced-rank field")
    else
      out$rho0_raw <- scal("rho0", function(v) abs(v) < 1, atanh,
                           "a number in (-1, 1)")
  }
  if (!is.null(fixed[["sigma0"]]))
    out$log_sigma0 <- scal("sigma0", function(v) v > 0, log,
                           "a positive number")
  out
}

# Which of the p + q + 3 fixed-effect parameters (plus rho0 and sigma0 of the
# inflation predictor) are estimated
.free_mask <- function(fx, p, q, zm = "none") {
  free_block <- function(b, len)
    if (is.null(fx[[b]])) rep(TRUE, len) else is.na(fx[[b]])
  extra <- .zero_extra(zm)
  c(free_block("beta", p), free_block("gamma", q), free_block("rho_raw", 1),
    free_block("log_kappa", 1), free_block("log_sigma", 1),
    unlist(lapply(extra, free_block, len = 1)))
}
