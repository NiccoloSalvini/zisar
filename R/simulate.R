#' Simulate data from the ZISAR model
#'
#' Draws one data set from the restricted ZISAR model, on a regular lattice
#' with queen contiguity (the Monte Carlo design of the ZISAR paper) or on a
#' weights matrix supplied by the user. Covariates are independent standard
#' normal; the first element of `beta` and of `gamma` is the intercept.
#'
#' @param n_side Side of the square lattice, so \eqn{n} = `n_side^2` units.
#'   Ignored when `W` is supplied.
#' @param beta Count-process coefficients, intercept first.
#' @param rho Spatial autoregressive parameter, in \eqn{(-1, 1)}.
#' @param gamma Zero-inflation coefficients, intercept first.
#' @param kappa Negative binomial shape.
#' @param sigma Standard deviation of the random effects.
#' @param W Optional row-standardised weights matrix (or `spdep` `listw`)
#'   to use instead of the lattice.
#'
#' @return A list with
#'   \item{data}{Data frame with the counts `y`, the count covariates
#'     `x1, x2, ...` and the inflation covariates `z1, z2, ...`, ready for
#'     `zisar(y ~ x1 + ... | z1 + ..., data = data, W = W)`.}
#'   \item{W}{The weights matrix (sparse).}
#'   \item{coords}{Lattice coordinates (`row`, `col`), or `NULL` with a
#'     user-supplied `W`.}
#'   \item{truth}{The parameters used.}
#'   \item{eta, mu, pi, s}{Latent field, negative binomial means,
#'     inflation probabilities and structural-zero indicators.}
#' @export
#' @examples
#' set.seed(42)
#' sim <- rzisar(15, beta = c(1, 0.5, -0.3), rho = 0.5, gamma = c(-1, 0.7))
#' head(sim$data)
#' mean(sim$data$y == 0)
rzisar <- function(n_side, beta, rho, gamma, kappa = 2, sigma = 0.5,
                   W = NULL) {
  stopifnot(length(beta) >= 1, length(gamma) >= 1, abs(rho) < 1,
            kappa > 0, sigma >= 0)
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

  draw <- .rzisar_core(X, Z, W, beta, gamma, rho, kappa, sigma)
  data <- data.frame(y = draw$y, X[, -1, drop = FALSE], Z[, -1, drop = FALSE],
                     check.names = FALSE)

  list(data = data, W = W, coords = coords,
       truth = list(beta = beta, rho = rho, gamma = gamma, kappa = kappa,
                    sigma = sigma),
       eta = draw$eta, mu = draw$mu, pi = draw$pi, s = draw$s)
}

# One draw of (eta, s, y) given design matrices and parameters
.rzisar_core <- function(X, Z, W, beta, gamma, rho, kappa, sigma,
                         offset = 0) {
  n <- nrow(X)
  u <- rnorm(n, 0, sigma)
  A <- Matrix::Diagonal(n) - rho * W
  eta <- as.numeric(Matrix::solve(A, as.numeric(X %*% beta) + u))
  mu <- exp(eta + offset)
  pi <- as.numeric(plogis(Z %*% gamma))
  s <- rbinom(n, 1, pi)
  y_count <- rnbinom(n, size = kappa, mu = mu)
  list(y = ifelse(s == 1, 0L, y_count), eta = eta, mu = mu, pi = pi, s = s)
}
