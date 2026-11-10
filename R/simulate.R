#' Simulate data from the ZISAR model
#'
#' Draws one data set from the ZISAR model, on a regular lattice with queen
#' contiguity (the Monte Carlo design of the ZISAR paper) or on a weights
#' matrix supplied by the user. Covariates are independent standard normal;
#' the first element of `beta` and of `gamma` is the intercept.
#' `zero_spatial` adds spatial structure to the inflation predictor: a SAR
#' field with one element per unit (`"field"`, the general model of the
#' paper) or a SAR filter on \eqn{Z\gamma} (`"sar"`); see [zisar()].
#'
#' @param n_side Side of the square lattice, so \eqn{n} = `n_side^2` units.
#'   Ignored when `W` is supplied.
#' @param beta Count-process coefficients, intercept first.
#' @param rho Spatial autoregressive parameter, in \eqn{(-1, 1)}.
#' @param gamma Zero-inflation coefficients, intercept first.
#' @param kappa Negative binomial shape.
#' @param sigma Standard deviation of the random effects.
#' @param zero_spatial Inflation predictor: `"none"` (restricted model),
#'   `"field"` or `"sar"`.
#' @param rho0,sigma0 Spatial autoregressive parameter and standard deviation
#'   of the inflation field (`"field"`), or parameter of the filter
#'   (`"sar"`, which takes no `sigma0`).
#' @param W Optional row-standardised weights matrix (or `spdep` `listw`)
#'   to use instead of the lattice.
#'
#' @return A list with
#'   \item{data}{Data frame with the counts `y`, the count covariates
#'     `x1, x2, ...` and the inflation covariates `z1, z2, ...`, ready for
#'     `zisar(y ~ x1 + ... | z1 + ..., data = data, W = W)`.}
#'   \item{W}{The weights matrix (sparse).}
#'   \item{coords}{Lattice coordinates (`row`, `col`), for example for
#'     [zisar_basis()]; `NULL` with a user-supplied `W`.}
#'   \item{truth}{The parameters used.}
#'   \item{eta, g, zeta, mu, pi, s}{Latent field, inflation field (zero
#'     without one), inflation predictor, negative binomial means, inflation
#'     probabilities and structural-zero indicators.}
#' @export
#' @examples
#' set.seed(42)
#' sim <- rzisar(15, beta = c(1, 0.5, -0.3), rho = 0.5, gamma = c(-1, 0.7))
#' head(sim$data)
#' mean(sim$data$y == 0)
#'
#' # Structural zeros clustered by a field of their own
#' sim2 <- rzisar(15, beta = c(1, 0.5), rho = 0.5, gamma = c(-1, 0.7),
#'                zero_spatial = "field", rho0 = 0.8, sigma0 = 1)
#' cor(sim2$g, as.numeric(sim2$W %*% sim2$g))
rzisar <- function(n_side, beta, rho, gamma, kappa = 2, sigma = 0.5,
                   zero_spatial = c("none", "field", "sar"), rho0 = 0,
                   sigma0 = 0, W = NULL) {
  zero_spatial <- match.arg(zero_spatial)
  stopifnot(length(beta) >= 1, length(gamma) >= 1, abs(rho) < 1,
            kappa > 0, sigma >= 0, abs(rho0) < 1, sigma0 >= 0)
  if (zero_spatial == "none" && (rho0 != 0 || sigma0 != 0))
    stop("`rho0` and `sigma0` need `zero_spatial = \"field\"` or `\"sar\"`.",
         call. = FALSE)
  if (zero_spatial == "sar" && sigma0 != 0)
    stop("The SAR filter (`zero_spatial = \"sar\"`) has no `sigma0`.",
         call. = FALSE)
  if (is.null(W)) {
    W <- lattice_weights(n_side)
    coords <- expand.grid(row = seq_len(n_side), col = seq_len(n_side))
  } else {
    W <- .as_weights(W, if (inherits(W, "listw")) length(W$neighbours) else
      nrow(W))
    coords <- NULL
  }
  n <- nrow(W)
  p <- length(beta)
  q <- length(gamma)

  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  Z <- cbind(1, matrix(rnorm(n * (q - 1)), n, q - 1))
  colnames(X) <- c("(Intercept)", if (p > 1) paste0("x", seq_len(p - 1)))
  colnames(Z) <- c("(Intercept)", if (q > 1) paste0("z", seq_len(q - 1)))

  draw <- .rzisar_core(X, Z, W, beta, gamma, rho, kappa, sigma,
                       zero = list(mode = zero_spatial, rho0 = rho0,
                                   sigma0 = sigma0))
  data <- data.frame(y = draw$y, X[, -1, drop = FALSE], Z[, -1, drop = FALSE],
                     check.names = FALSE)

  list(data = data, W = W, coords = coords,
       truth = list(beta = beta, rho = rho, gamma = gamma, kappa = kappa,
                    sigma = sigma, zero_spatial = zero_spatial, rho0 = rho0,
                    sigma0 = sigma0),
       eta = draw$eta, g = draw$g, zeta = draw$zeta, mu = draw$mu,
       pi = draw$pi, s = draw$s)
}

# One draw of (eta, g, s, y) given design matrices and parameters. `zero`
# describes the inflation predictor: mode ("none", "field", "sar" or
# "reduced"), rho0, sigma0 and, for "reduced", the basis. Random numbers are
# drawn in the order u, inflation field, zero indicators, counts.
.rzisar_core <- function(X, Z, W, beta, gamma, rho, kappa, sigma,
                         offset = 0, zero = list(mode = "none")) {
  n <- nrow(X)
  I_n <- Matrix::Diagonal(n)
  u <- rnorm(n, 0, sigma)
  eta <- as.numeric(Matrix::solve(I_n - rho * W, as.numeric(X %*% beta) + u))
  mu <- exp(eta + offset)
  g <- switch(zero$mode,
    field = as.numeric(Matrix::solve(I_n - zero$rho0 * W,
                                     rnorm(n, 0, zero$sigma0))),
    reduced = {
      b <- zero$basis
      alpha <- Matrix::solve(Matrix::Diagonal(b$K) - zero$rho0 * b$Wk,
                             rnorm(b$K, 0, zero$sigma0))
      as.numeric(b$B %*% alpha)
    },
    rep(0, n))
  zeta <- as.numeric(Z %*% gamma) + g
  if (zero$mode == "sar")
    zeta <- as.numeric(Matrix::solve(I_n - zero$rho0 * W, zeta))
  pi <- plogis(zeta)
  s <- rbinom(n, 1, pi)
  y_count <- rnbinom(n, size = kappa, mu = mu)
  list(y = ifelse(s == 1, 0L, y_count), eta = eta, g = g, zeta = zeta,
       mu = mu, pi = pi, s = s)
}
