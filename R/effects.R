#' Direct, indirect and total effects of a ZISAR model
#'
#' Decomposes the effect of each covariate on the expected count
#' \eqn{E[y] = (1 - \pi)\mu} into a direct part (a unit's own covariate), an
#' indirect part (spillovers from the covariates of other units through the
#' spatial multiplier) and their sum, averaged over units as in LeSage and
#' Pace (2009).
#'
#' The matrix of partial derivatives of \eqn{E[y_i]} with respect to
#' \eqn{x_{jr}} is \eqn{\mathrm{diag}(d)\, S(\rho)\, \beta_r} with
#' \eqn{d_i = (1 - \pi_i)\mu_i} and \eqn{S(\rho) = (I - \rho W)^{-1}}. The
#' average direct effect is the mean of its diagonal, the average total
#' effect the mean of its row sums, and the indirect effect their
#' difference. Only the diagonal and the row sums of \eqn{S(\rho)} are
#' computed, so no dense \eqn{n \times n} matrix is formed for large \eqn{n}.
#'
#' A covariate that also enters the zero-inflation part (same name in both
#' parts of the formula) adds the inflation channel
#' \eqn{-\mu_i \pi_i (1 - \pi_i) \gamma_r}, which acts on the unit itself
#' only and so enters the direct and total effects. A covariate of the
#' inflation part alone has this direct effect only. With the SAR filter on
#' the inflation predictor (`zero_spatial = "sar"`) the inflation channel
#' passes through \eqn{S_0 = (I - \rho_0 W)^{-1}} as well, and has an
#' indirect part.
#'
#' Standard errors and intervals are simulation-based: parameters are drawn
#' from their estimated asymptotic normal distribution and the effects
#' recomputed for each draw, keeping the predicted random effects
#' \eqn{\hat u = (I - \rho W)\hat\eta - X\beta} (and an inflation field
#' at its mode) at their estimates. For \eqn{n > 3000}, the diagonal of
#' \eqn{S(\rho)} along the draws is interpolated between seven exact
#' evaluations.
#'
#' @param object A fitted `"zisar"` object.
#' @param B Number of simulation draws for standard errors and intervals;
#'   `0` gives point estimates only.
#' @param level Confidence level of the simulation intervals.
#'
#' @return A data frame of class `"zisar_effects"`, one row per covariate and
#'   effect (`direct`, `indirect`, `total`), with the point estimate and,
#'   when `B > 0`, the simulation standard error and interval bounds.
#'
#' @references
#' LeSage, J. and Pace, R. K. (2009). *Introduction to Spatial
#' Econometrics*. Chapman and Hall/CRC.
#' @export
#' @examples
#' set.seed(1)
#' sim <- rzisar(12, beta = c(1, 0.5), rho = 0.4, gamma = c(-1, 0.8))
#' fit <- zisar(y ~ x1 | z1, data = sim$data, W = sim$W)
#' zisar_effects(fit, B = 200)
zisar_effects <- function(object, B = 500, level = 0.95) {
  stopifnot(inherits(object, "zisar"), B >= 0, level > 0, level < 1)
  n <- object$n
  W <- object$W
  X <- object$X
  Z <- object$Z
  p <- ncol(X)
  q <- ncol(Z)
  I_n <- Matrix::Diagonal(n)
  zsar <- .zmode(object) == "sar"

  # Covariates (intercepts excluded) and their columns in X and Z
  terms <- union(.covariates(X), .covariates(Z))
  if (!length(terms)) stop("The model has no covariates.", call. = FALSE)
  bx <- match(terms, colnames(X))
  gz <- match(terms, colnames(Z))
  pick <- function(v, idx) ifelse(is.na(idx), 0, v[idx])

  # s0_ii: diagonal of S0 = (I - rho0 W)^{-1} with the SAR filter on the
  # inflation predictor, 1 otherwise
  avg <- function(beta, gamma, rho, rho0, eta, s_ii, s0_ii) {
    mu <- exp(eta + object$offset)
    pi <- plogis(.zeta(object, Z, gamma, rho0))
    d <- (1 - pi) * mu
    s_1 <- as.numeric(Matrix::solve(I_n - rho * W, rep(1, n)))
    s0_1 <- if (zsar) as.numeric(Matrix::solve(I_n - rho0 * W, rep(1, n)))
            else 1
    zc <- mu * pi * (1 - pi)
    b <- pick(beta, bx)
    g <- pick(gamma, gz)
    direct <- mean(d * s_ii) * b - mean(zc * s0_ii) * g
    total <- mean(d * s_1) * b - mean(zc * s0_1) * g
    cbind(direct = direct, indirect = total - direct, total = total)
  }

  beta <- object$coefficients$count
  gamma <- object$coefficients$zero
  rho <- object$rho
  rho0 <- object$rho0
  point <- avg(beta, gamma, rho, rho0, object$eta, .sar_diag(W, rho),
               if (zsar) .sar_diag(W, rho0) else 1)

  out <- data.frame(term = rep(terms, 3),
                    effect = rep(c("direct", "indirect", "total"),
                                 each = length(terms)),
                    estimate = as.vector(point))

  V <- object$vcov_internal
  if (B > 0 && !all(is.finite(V))) {
    warning("No standard errors: the covariance matrix is not available.",
            call. = FALSE)
    B <- 0
  }
  draws <- if (B > 0) tryCatch(
    matrix(MASS::mvrnorm(B, object$par[object$free], V), B),
    error = function(e) NULL)
  if (B > 0 && is.null(draws)) {
    warning("No standard errors: the covariance matrix is not positive ",
            "definite.", call. = FALSE)
    B <- 0
  }
  if (B > 0) {
    full <- matrix(object$par, B, length(object$par), byrow = TRUE)
    full[, object$free] <- draws
    rho_draws <- tanh(full[, p + q + 1])
    u_hat <- as.numeric((I_n - rho * W) %*% object$eta) -
      as.numeric(X %*% beta)
    s_ii_at <- .sar_diag_along(W, rho_draws, rho)
    if (zsar) {
      rho0_draws <- tanh(full[, p + q + 4])
      s0_ii_at <- .sar_diag_along(W, rho0_draws, rho0)
    }

    sims <- array(NA_real_, c(B, length(terms), 3))
    for (b in seq_len(B)) {
      beta_b <- full[b, seq_len(p)]
      eta_b <- as.numeric(Matrix::solve(I_n - rho_draws[b] * W,
                                        as.numeric(X %*% beta_b) + u_hat))
      sims[b, , ] <- avg(beta_b, full[b, p + seq_len(q)], rho_draws[b],
                         if (zsar) rho0_draws[b], eta_b,
                         s_ii_at(rho_draws[b]),
                         if (zsar) s0_ii_at(rho0_draws[b]) else 1)
    }
    # Draws with an overflowing mean are dropped
    valid <- apply(sims, 1, function(m) all(is.finite(m)))
    if (sum(valid) < min(B, 10)) {
      warning("No standard errors: too few draws gave finite effects.",
              call. = FALSE)
      B <- 0
    } else {
      sims <- sims[valid, , , drop = FALSE]
      alpha <- (1 - level) / 2
      out$se <- as.vector(apply(sims, c(2, 3), sd))
      out$lower <- as.vector(apply(sims, c(2, 3), quantile, alpha))
      out$upper <- as.vector(apply(sims, c(2, 3), quantile, 1 - alpha))
    }
  }
  out <- out[order(match(out$term, terms)), ]
  rownames(out) <- NULL
  structure(out, class = c("zisar_effects", "data.frame"), B = B,
            level = level)
}

#' @rdname zisar_effects
#' @param x A `"zisar_effects"` object.
#' @param digits Number of decimal places to print.
#' @param ... Not used.
#' @export
print.zisar_effects <- function(x, digits = 4, ...) {
  B <- attr(x, "B")
  cat("ZISAR average effects on E[y]")
  if (B > 0) cat(sprintf(" (simulation SEs, B = %d)", B))
  cat("\n\n")
  terms <- unique(x$term)
  tab <- sapply(c("direct", "indirect", "total"), function(e) {
    xe <- x[x$effect == e, ][match(terms, x$term[x$effect == e]), ]
    est <- formatC(xe$estimate, digits = digits, format = "f")
    if (B > 0) paste0(est, " (", formatC(xe$se, digits = digits,
                                         format = "f"), ")") else est
  })
  tab <- matrix(tab, nrow = length(terms),
                dimnames = list(terms, c("Direct", "Indirect", "Total")))
  print(tab, quote = FALSE, right = TRUE)
  invisible(x)
}

.covariates <- function(M) {
  k <- .intercept_col(M)
  if (is.na(k)) colnames(M) else colnames(M)[-k]
}

# diag((I - rho W)^{-1}): dense for small n, sparse LU solves in blocks
# otherwise
.sar_diag <- function(W, rho, block = 2000) {
  n <- nrow(W)
  A <- Matrix::Diagonal(n) - rho * W
  if (n <= 3000) return(diag(solve(as.matrix(A))))
  luA <- Matrix::lu(A)
  out <- numeric(n)
  for (s in seq(1, n, by = block)) {
    j <- s:min(n, s + block - 1)
    E <- matrix(0, n, length(j))
    E[cbind(j, seq_along(j))] <- 1
    Sj <- Matrix::solve(luA, E)
    out[j] <- Sj[cbind(j, seq_along(j))]
  }
  out
}

# Function giving diag(S(rho)) along simulation draws of rho: exact for small
# n; for large n, exact at Chebyshev nodes spanning the draws and Lagrange
# interpolation in between (S_ii is smooth in rho)
.sar_diag_along <- function(W, rho_draws, rho) {
  if (nrow(W) <= 3000) return(function(r0) .sar_diag(W, r0))
  a <- min(rho_draws)
  b <- max(rho_draws)
  if (b - a < 1e-10) {
    s0 <- .sar_diag(W, rho)
    return(function(r0) s0)
  }
  nodes <- (a + b) / 2 + (b - a) / 2 * cos((2 * (1:7) - 1) * pi / 14)
  s_nodes <- sapply(nodes, function(r0) .sar_diag(W, r0))
  lagrange <- function(x) sapply(seq_along(nodes), function(k)
    prod((x - nodes[-k]) / (nodes[k] - nodes[-k])))
  function(r0) as.numeric(s_nodes %*% lagrange(r0))
}
