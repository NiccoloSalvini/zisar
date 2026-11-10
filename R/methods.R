#' Methods for fitted ZISAR models
#'
#' Standard accessors for objects returned by [zisar()] and [zisar_fit()].
#'
#' `coef()` and `vcov()` cover the regression coefficients and \eqn{\rho}
#' (and \eqn{\rho_0} when the inflation predictor has a field or a filter);
#' the negative binomial shape \eqn{\kappa} and the standard deviations
#' \eqn{\sigma} and \eqn{\sigma_0} are reported by `summary()` and stored
#' in `object$kappa`, `object$sigma` and `object$sigma0`. The covariance of
#' \eqn{\rho} and \eqn{\rho_0} comes from the delta method; rows and
#' columns of parameters held fixed are `NA`.
#'
#' `logLik()` returns the Laplace-approximated marginal log-likelihood, with
#' degrees of freedom equal to the number of estimated parameters, so
#' [stats::AIC()] and [stats::BIC()] work on fitted models.
#'
#' @param object,x A fitted `"zisar"` object.
#' @param model Which coefficients: `"full"` (count, zero-inflation,
#'   \eqn{\rho} and, with a spatial inflation predictor, \eqn{\rho_0}),
#'   `"count"` or `"zero"`.
#' @param digits Number of significant digits to print.
#' @param ... Not used.
#'
#' @return `coef()` a named vector, `vcov()` a matrix, `logLik()` an object
#'   of class `"logLik"`, `nobs()` the number of spatial units, `fitted()`
#'   the fitted means \eqn{E[y_i] = (1 - \pi_i)\mu_i}, `summary()` an object
#'   of class `"summary.zisar"`.
#' @name zisar-methods
#' @examples
#' set.seed(1)
#' sim <- rzisar(10, beta = c(1, 0.5), rho = 0.4, gamma = c(-1, 0.8))
#' fit <- zisar(y ~ x1 | z1, data = sim$data, W = sim$W)
#' coef(fit)
#' sqrt(diag(vcov(fit)))
#' AIC(fit)
NULL

#' @rdname zisar-methods
#' @export
coef.zisar <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  cf <- object$coefficients
  switch(model,
    count = cf$count,
    zero = cf$zero,
    full = c(stats::setNames(cf$count, paste0("count_", names(cf$count))),
             stats::setNames(cf$zero, paste0("zero_", names(cf$zero))),
             rho = object$rho,
             if (.zmode(object) != "none") c(rho0 = object$rho0))
  )
}

#' @rdname zisar-methods
#' @export
vcov.zisar <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  p <- length(object$coefficients$count)
  q <- length(object$coefficients$zero)
  zs <- .zmode(object) != "none"
  keep <- c(seq_len(p + q + 1), if (zs) p + q + 4)
  V <- .vcov_full_internal(object)[keep, keep, drop = FALSE]
  # d rho / d rho_raw = 1 - rho^2; elementwise, so NA rows do not spread
  jac <- c(rep(1, p + q), 1 - object$rho^2, if (zs) 1 - object$rho0^2)
  V <- V * outer(jac, jac)
  nm <- names(coef(object, "full"))
  dimnames(V) <- list(nm, nm)
  idx <- switch(model, full = seq_along(keep), count = seq_len(p),
                zero = p + seq_len(q))
  V <- V[idx, idx, drop = FALSE]
  if (model != "full") dimnames(V) <- list(names(coef(object, model)),
                                           names(coef(object, model)))
  V
}

# Covariance of all internal parameters, NA where fixed
.vcov_full_internal <- function(object) {
  nm <- names(object$par)
  V <- matrix(NA_real_, length(nm), length(nm), dimnames = list(nm, nm))
  V[object$free, object$free] <- object$vcov_internal
  V
}

#' @rdname zisar-methods
#' @export
logLik.zisar <- function(object, ...) {
  structure(object$loglik, df = object$df, nobs = object$n, class = "logLik")
}

#' @rdname zisar-methods
#' @export
nobs.zisar <- function(object, ...) object$n

#' @rdname zisar-methods
#' @export
fitted.zisar <- function(object, ...) object$fitted.values

#' @rdname zisar-methods
#' @export
print.zisar <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!is.null(x$call)) {
    cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n", sep = "")
  }
  cat("\nCount model coefficients (SAR negative binomial, log link):\n")
  print(format(x$coefficients$count, digits = digits), print.gap = 2L,
        quote = FALSE)
  cat("\nZero-inflation model coefficients (binomial, logit link):\n")
  print(format(x$coefficients$zero, digits = digits), print.gap = 2L,
        quote = FALSE)
  cat(sprintf("\nrho = %s   kappa = %s   sigma = %s\n",
              format(x$rho, digits = digits), format(x$kappa, digits = digits),
              format(x$sigma, digits = digits)))
  zm <- .zmode(x)
  if (zm == "sar")
    cat(sprintf("Inflation SAR filter: rho0 = %s\n",
                format(x$rho0, digits = digits)))
  if (zm %in% c("field", "reduced"))
    cat(sprintf("Inflation field (%s): rho0 = %s   sigma0 = %s\n",
                if (zm == "reduced") paste(x$basis$K, "knots") else "unit level",
                format(x$rho0, digits = digits),
                format(x$sigma0, digits = digits)))
  cat(sprintf("Log-likelihood: %s on %d Df   (n = %d)\n\n",
              format(x$loglik, digits = digits + 2L), x$df, x$n))
  invisible(x)
}

#' @rdname zisar-methods
#' @export
summary.zisar <- function(object, ...) {
  se <- sqrt(diag(.vcov_full_internal(object)))
  p <- length(object$coefficients$count)
  q <- length(object$coefficients$zero)
  wald <- function(est, se) {
    z <- est / se
    cbind(Estimate = est, `Std. Error` = se, `z value` = z,
          `Pr(>|z|)` = 2 * pnorm(-abs(z)))
  }
  count <- wald(object$coefficients$count, se[seq_len(p)])
  zero <- wald(object$coefficients$zero, se[p + seq_len(q)])
  rownames(count) <- names(object$coefficients$count)
  rownames(zero) <- names(object$coefficients$zero)

  nat <- c(rho = object$rho, kappa = object$kappa, sigma = object$sigma,
           rho0 = object$rho0, sigma0 = object$sigma0)
  spatial <- wald(nat, object$se_natural[names(nat)])
  no_test <- intersect(c("kappa", "sigma", "sigma0"), names(nat))
  spatial[no_test, c("z value", "Pr(>|z|)")] <- NA
  rownames(spatial) <- names(nat)

  ll <- logLik(object)
  user_names <- c(rho_raw = "rho", log_kappa = "kappa", log_sigma = "sigma",
                  rho0_raw = "rho0", log_sigma0 = "sigma0")
  fixed <- names(object$free)[!object$free]
  fixed[fixed %in% names(user_names)] <- user_names[fixed[fixed %in%
                                                            names(user_names)]]
  structure(
    list(call = object$call, count = count, zero = zero, spatial = spatial,
         zero_spatial = .zmode(object), K = object$basis$K, fixed = fixed,
         loglik = object$loglik, df = object$df, n = object$n,
         zeros = mean(object$y == 0),
         aic = stats::AIC(ll), bic = stats::BIC(ll),
         logdet = object$logdet, convergence = object$convergence,
         laplace_ok = object$laplace_ok, max_grad = object$max_grad,
         message = object$message, iterations = object$iterations),
    class = "summary.zisar"
  )
}

#' @rdname zisar-methods
#' @param signif.stars Logical: show significance stars? See
#'   [stats::printCoefmat()].
#' @export
print.summary.zisar <- function(x, digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"),
                                ...) {
  if (!is.null(x$call)) {
    cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n", sep = "")
  }
  cat(sprintf("\n%d spatial units, %.1f%% zeros\n", x$n, 100 * x$zeros))
  cat("\nCount model coefficients (SAR negative binomial, log link):\n")
  stats::printCoefmat(x$count, digits = digits, signif.stars = signif.stars,
                      na.print = "NA")
  cat("\nZero-inflation model coefficients (binomial, logit link):\n")
  stats::printCoefmat(x$zero, digits = digits, signif.stars = signif.stars,
                      na.print = "NA")
  cat("\nSpatial and dispersion parameters",
      switch(x$zero_spatial,
             none = "",
             sar = " (rho0: SAR filter of the inflation predictor)",
             field = " (rho0, sigma0: unit-level inflation field)",
             reduced = sprintf(paste0(" (rho0, sigma0: reduced-rank inflation ",
                                      "field, %d knots)"), x$K)),
      ":\n", sep = "")
  stats::printCoefmat(x$spatial, digits = digits, signif.stars = signif.stars,
                      na.print = "", has.Pvalue = TRUE, P.values = TRUE)
  if (length(x$fixed))
    cat("Held fixed:", paste(x$fixed, collapse = ", "), "\n")
  cat(sprintf("\nLog-likelihood (Laplace): %s on %d Df\n",
              format(x$loglik, digits = digits + 2L), x$df))
  cat(sprintf("AIC: %s   BIC: %s\n", format(x$aic, digits = digits + 2L),
              format(x$bic, digits = digits + 2L)))
  cat(sprintf("Log-determinant: %s.  Optimiser: %s after %d iterations.\n\n",
              x$logdet,
              if (x$convergence == 0) "converged" else
                paste0("NOT converged (", x$message, ")"),
              x$iterations))
  if (!x$laplace_ok)
    cat("WARNING: the Laplace approximation is unreliable at this optimum",
        sprintf("(max gradient %s); see ?zisar.\n\n", signif(x$max_grad, 3)))
  invisible(x)
}

#' Predictions from a fitted ZISAR model
#'
#' Fitted or counterfactual means, inflation probabilities and latent field
#' of a ZISAR model.
#'
#' Without `newdata`, predictions are at the fitted values, with the latent
#' field \eqn{\eta} at its mode given the data.
#'
#' With `newdata`, predictions are a counterfactual scenario for the *same*
#' spatial units: `newdata` holds one row per unit, in the order of the
#' weights matrix, with some covariates changed. The latent field is
#' recomputed through the spatial multiplier,
#' \eqn{\eta^* = (I - \rho W)^{-1}(X^*\beta + \hat u)}, keeping the predicted
#' random effects \eqn{\hat u = (I - \rho W)\hat\eta - X\beta} at their
#' estimates, so a change in one unit spills over to its neighbours. An
#' inflation field \eqn{g} is kept at its mode; with the SAR filter
#' (`zero_spatial = "sar"`) the inflation predictor is recomputed through
#' \eqn{(I - \rho_0 W)^{-1}}, so changes in the inflation covariates spill
#' over too.
#' Predictions for new units would need a new weights matrix and are not
#' supported.
#'
#' @param object A fitted `"zisar"` object.
#' @param newdata Optional data frame of covariates for the same units, in
#'   the order of `W`; needs a model fitted with [zisar()].
#' @param type `"response"` for the mean \eqn{E[y_i] = (1 - \pi_i)\mu_i},
#'   `"count"` for the mean of the negative binomial component \eqn{\mu_i},
#'   `"zero"` for the probability of a structural zero \eqn{\pi_i}, or
#'   `"eta"` for the latent field \eqn{\eta_i} (without offset).
#' @param ... Not used.
#'
#' @return A numeric vector with one value per spatial unit.
#' @export
#' @examples
#' set.seed(1)
#' sim <- rzisar(10, beta = c(1, 0.5), rho = 0.4, gamma = c(-1, 0.8))
#' fit <- zisar(y ~ x1 | z1, data = sim$data, W = sim$W)
#' head(predict(fit, type = "zero"))
#'
#' # Raise x1 by one in unit 1 only: the change spills over to its neighbours
#' nd <- sim$data
#' nd$x1[1] <- nd$x1[1] + 1
#' delta <- predict(fit, newdata = nd) - predict(fit)
#' round(head(delta), 3)
predict.zisar <- function(object, newdata = NULL,
                          type = c("response", "count", "zero", "eta"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    eta <- object$eta
    mu <- object$mu
    pi <- object$pi
  } else {
    if (is.null(object$terms))
      stop("`newdata` needs a model fitted with zisar() (formula interface).",
           call. = FALSE)
    if (nrow(newdata) != object$n)
      stop(sprintf(paste0("`newdata` must have one row per spatial unit (%d), ",
                          "in the order of `W`."), object$n), call. = FALSE)
    des <- .new_design(object, newdata)
    A <- Matrix::Diagonal(object$n) - object$rho * object$W
    beta <- object$coefficients$count
    u_hat <- as.numeric(A %*% object$eta) - as.numeric(object$X %*% beta)
    eta <- as.numeric(Matrix::solve(A, as.numeric(des$X %*% beta) + u_hat))
    mu <- exp(eta + des$offset)
    pi <- plogis(.zeta(object, des$Z))
  }
  switch(type,
    response = (1 - pi) * mu,
    count = mu,
    zero = pi,
    eta = eta
  )
}

# Inflation predictor of a fit as a string (a logical in fits from
# version 0.1)
.zmode <- function(object) {
  z <- object$zero_spatial
  if (is.null(z) || isFALSE(z)) "none" else if (isTRUE(z)) "field" else z
}

# Inflation predictor zeta at design Z and (gamma, rho0): Z gamma plus the
# mode of the field, or the SAR filter of Z gamma
.zeta <- function(object, Z = object$Z, gamma = object$coefficients$zero,
                  rho0 = object$rho0) {
  lin <- as.numeric(Z %*% gamma)
  switch(.zmode(object),
         none = lin,
         sar = as.numeric(Matrix::solve(Matrix::Diagonal(object$n) -
                                          rho0 * object$W, lin)),
         lin + object$g)
}

# Design matrices and offset of newdata, with the fitted terms and levels
.new_design <- function(object, newdata) {
  tc <- stats::delete.response(object$terms$count)
  mfc <- model.frame(tc, newdata, xlev = object$xlevels$count,
                     na.action = stats::na.pass)
  mfz <- model.frame(object$terms$zero, newdata, xlev = object$xlevels$zero,
                     na.action = stats::na.pass)
  if (anyNA(mfc) || anyNA(mfz))
    stop("`newdata` has missing values in the model variables.", call. = FALSE)
  off <- model.offset(mfc)
  list(X = model.matrix(tc, mfc, contrasts.arg = object$contrasts$count),
       Z = model.matrix(object$terms$zero, mfz,
                        contrasts.arg = object$contrasts$zero),
       offset = if (is.null(off)) rep(0, object$n) else off)
}

#' Residuals of a fitted ZISAR model
#'
#' Response residuals \eqn{y_i - E[y_i]} or Pearson residuals, both
#' conditional on the estimated latent field. Pearson residuals divide by
#' the standard deviation of the zero-inflated negative binomial mixture,
#' \eqn{\mathrm{Var}(y_i) = (1 - \pi_i)\mu_i(1 + \mu_i/\kappa + \pi_i\mu_i)}.
#'
#' Because the residuals condition on the mode of \eqn{\eta}, which already
#' absorbs most of the spatial structure, tests such as Moran's I on them
#' have little power against remaining spatial dependence; the
#' likelihood-ratio test of \eqn{\rho = 0} is preferable (see [zisar()]).
#'
#' @param object A fitted `"zisar"` object.
#' @param type `"pearson"` or `"response"`.
#' @param ... Not used.
#' @return A numeric vector of residuals.
#' @export
residuals.zisar <- function(object, type = c("pearson", "response"), ...) {
  type <- match.arg(type)
  res <- object$y - object$fitted.values
  if (type == "response") return(res)
  mu <- object$mu
  pi <- object$pi
  res / sqrt((1 - pi) * mu * (1 + mu / object$kappa + pi * mu))
}

#' Simulate responses from a fitted ZISAR model
#'
#' Draws new count vectors from the fitted model, for example for a
#' parametric bootstrap. Each draw takes new random effects \eqn{u} (and a
#' new inflation field, \eqn{g} or \eqn{\alpha}, when the model has one),
#' new structural-zero indicators and new negative binomial counts, at the
#' estimated parameters and the observed covariates.
#'
#' To check a large gradient at the optimum (see [zisar()]), refit the model
#' to a few simulated responses: if the refits show the same gradient while
#' recovering the parameters, it belongs to the model and not to a spurious
#' optimum. For bootstrap standard errors, start each refit from the
#' estimates with `start`.
#'
#' @param object A fitted `"zisar"` object.
#' @param nsim Number of response vectors to simulate.
#' @param seed Optional seed, passed to [set.seed()].
#' @param ... Not used.
#' @return A data frame with `nsim` columns, one simulated response vector
#'   each.
#' @export
#' @examples
#' set.seed(1)
#' sim <- rzisar(10, beta = c(1, 0.5), rho = 0.4, gamma = c(-1, 0.8))
#' fit <- zisar(y ~ x1 | z1, data = sim$data, W = sim$W)
#' ys <- simulate(fit, nsim = 100, seed = 2)
#' # Share of zeros: observed against simulated
#' c(observed = mean(sim$data$y == 0), simulated = mean(as.matrix(ys) == 0))
simulate.zisar <- function(object, nsim = 1, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  zero <- list(mode = .zmode(object), rho0 = object$rho0,
               sigma0 = object$sigma0, basis = object$basis)
  sims <- replicate(nsim, .rzisar_core(
    object$X, object$Z, object$W, object$coefficients$count,
    object$coefficients$zero, object$rho, object$kappa, object$sigma,
    object$offset, zero = zero)$y, simplify = FALSE)
  out <- as.data.frame(sims, col.names = paste0("sim_", seq_len(nsim)))
  attr(out, "seed") <- seed
  out
}
