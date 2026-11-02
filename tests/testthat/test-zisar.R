test_that("zisar() recovers the parameters of a simulated data set", {
  set.seed(123)
  sim <- rzisar(25, beta = c(1, 0.5, -0.3), rho = 0.5, gamma = c(-1, 0.8),
                kappa = 2, sigma = 0.4)
  fit <- zisar(y ~ x1 + x2 | z1, data = sim$data, W = sim$W)
  expect_s3_class(fit, "zisar")
  expect_equal(fit$convergence, 0)
  expect_equal(unname(fit$coefficients$count), c(1, 0.5, -0.3),
               tolerance = 0.2)
  expect_equal(unname(fit$coefficients$zero), c(-1, 0.8), tolerance = 0.4)
  expect_equal(fit$rho, 0.5, tolerance = 0.3)
  expect_true(all(is.finite(fit$vcov_internal)))
})

test_that("formula and matrix interfaces agree", {
  sim <- sim_small()
  fit <- fit_small()
  X <- cbind(1, as.matrix(sim$data[, c("x1", "x2")]))
  fit_m <- zisar_fit(sim$data$y, X, cbind(1, sim$data$z1), sim$W)
  expect_equal(fit_m$loglik, fit$loglik, tolerance = 1e-8)
  expect_equal(unname(fit_m$coefficients$count),
               unname(fit$coefficients$count), tolerance = 1e-6)
  expect_named(fit_m$coefficients$count, c("(Intercept)", "x1", "x2"))
})

test_that("eigenvalue and sparse log-determinants give the same fit", {
  sim <- sim_small()
  fit_s <- zisar(y ~ x1 + x2 | z1, data = sim$data, W = sim$W,
                 logdet = "sparse")
  expect_equal(fit_s$logdet, "sparse")
  expect_equal(fit_s$loglik, fit_small()$loglik, tolerance = 1e-6)
  expect_equal(fit_s$rho, fit_small()$rho, tolerance = 1e-4)
})

test_that("a formula without `|` uses the same terms in both parts", {
  f <- zisar:::.split_formula(y ~ x1 + log(x2))
  expect_equal(f$count, y ~ x1 + log(x2), ignore_attr = TRUE)
  expect_equal(f$zero, ~ x1 + log(x2), ignore_attr = TRUE)
  expect_false(f$bar)
  f <- zisar:::.split_formula(y ~ x1 + x2 | z1)
  expect_equal(f$zero, ~ z1, ignore_attr = TRUE)

  sim <- sim_small()
  fit <- zisar(y ~ x1 + x2 + z1, data = sim$data, W = sim$W)
  expect_named(fit$coefficients$zero, c("(Intercept)", "x1", "x2", "z1"))
  fit1 <- zisar(y ~ x1 + x2 | 1, data = sim$data, W = sim$W)
  expect_named(fit1$coefficients$zero, "(Intercept)")
})

test_that("`fixed` gives the nested models", {
  f <- fits_rho0()
  expect_true(f$nested$laplace_ok)
  expect_equal(f$nested$rho, 0)
  expect_equal(f$nested$df, f$full$df - 1)
  expect_true(is.na(f$nested$se_natural[["rho"]]))
  expect_lte(f$nested$loglik, f$full$loglik + 1e-6)

  sim <- sim_small()
  # Only the second gamma held fixed
  fitg <- zisar(y ~ x1 + x2 | z1, data = sim$data, W = sim$W,
                fixed = list(gamma = c(NA, 0)))
  expect_equal(unname(fitg$coefficients$zero[2]), 0)
  expect_equal(fitg$df, fit_small()$df - 1)

  fitk <- zisar(y ~ x1 + x2 | z1, data = sim$data, W = sim$W,
                fixed = list(kappa = 3, sigma = 0.5))
  expect_equal(c(fitk$kappa, fitk$sigma), c(3, 0.5))
})

test_that("rho = 0 and sigma near 0 reproduce the aspatial ZINB of pscl", {
  skip_if_not_installed("pscl")
  sim <- sim_small()
  fit <- zisar(y ~ x1 + x2 | z1, data = sim$data, W = sim$W,
               fixed = list(rho = 0, sigma = 0.001))
  ref <- pscl::zeroinfl(y ~ x1 + x2 | z1, data = sim$data, dist = "negbin")
  expect_equal(fit$loglik, as.numeric(logLik(ref)), tolerance = 1e-4)
  expect_equal(unname(c(fit$coefficients$count, fit$coefficients$zero)),
               unname(coef(ref)), tolerance = 1e-3)
  expect_equal(fit$kappa, ref$theta, tolerance = 1e-3)
})

test_that("an offset enters the count log-mean with coefficient one", {
  set.seed(5)
  sim <- rzisar(15, beta = c(0, 0.5), rho = 0.3, gamma = c(-1.5))
  expo <- exp(rnorm(nrow(sim$data), 1, 0.5))
  # Rescaled counts: same structure with mean multiplied by the exposure
  d <- sim$data
  d$y <- rnbinom(nrow(d), size = 2, mu = expo * sim$mu) * (1 - sim$s)
  d$expo <- expo
  fit <- zisar(y ~ x1 + offset(log(expo)) | 1, data = d, W = sim$W)
  expect_equal(fit$offset, log(expo))
  expect_equal(unname(fit$coefficients$count[2]), 0.5, tolerance = 0.3)
  expect_equal(fitted(fit), (1 - fit$pi) * exp(fit$eta + log(expo)))
})

test_that("invalid input is rejected with informative errors", {
  sim <- sim_small()
  d <- sim$data
  d$x1[3] <- NA
  expect_error(zisar(y ~ x1 | z1, data = d, W = sim$W), "Missing values")
  d <- sim$data
  d$y[1] <- -1
  expect_error(zisar(y ~ x1 | z1, data = d, W = sim$W), "non-negative")
  d$y[1] <- 1.5
  expect_error(zisar(y ~ x1 | z1, data = d, W = sim$W), "integer")
  expect_error(zisar(y ~ x1 | z1, data = sim$data, W = sim$W[-1, -1]),
               "observations")
  expect_error(zisar(~ x1, data = sim$data, W = sim$W), "two-sided")
  expect_error(zisar(y ~ x1 | z1, data = sim$data, W = sim$W,
                     fixed = list(lambda = 1)), "Unknown parameter")
  expect_error(zisar(y ~ x1 | z1, data = sim$data, W = sim$W,
                     fixed = list(rho = 1)), "\\(-1, 1\\)")
  expect_error(zisar(y ~ x1 | z1, data = sim$data, W = sim$W,
                     fixed = list(beta = 1)), "length 2")
  expect_warning(
    zisar(y ~ x1 + x2 | z1 + offset(x1), data = sim$data, W = sim$W),
    "ignored"
  )
})

test_that("an spdep listw object can be passed as W", {
  skip_if_not_installed("spdep")
  sim <- sim_small()
  lw <- spdep::mat2listw(as.matrix(sim$W), style = "W")
  fit <- zisar(y ~ x1 + x2 | z1, data = sim$data, W = lw)
  expect_equal(fit$loglik, fit_small()$loglik, tolerance = 1e-8)
})

test_that("a breakdown of the Laplace approximation is reported", {
  # Small sample, omitted covariate: sigma inflates and the approximation
  # produces spurious optima from every start
  sim <- sim_small(n_side = 12)
  expect_warning(
    expect_warning(
      fit <- zisar(y ~ x1 | z1, data = sim$data, W = sim$W,
                   fixed = list(rho = 0)),
      "Hessian"
    ),
    "unreliable"
  )
  expect_false(fit$laplace_ok)
  expect_output(print(summary(fit)), "unreliable")
})
