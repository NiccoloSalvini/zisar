test_that("coef() and vcov() are consistent", {
  fit <- fit_small()
  cf <- coef(fit)
  expect_named(cf, c("count_(Intercept)", "count_x1", "count_x2",
                     "zero_(Intercept)", "zero_z1", "rho"))
  expect_equal(coef(fit, "count"), fit$coefficients$count)
  expect_equal(coef(fit, "zero"), fit$coefficients$zero)

  V <- vcov(fit)
  expect_equal(dim(V), c(6, 6))
  expect_equal(rownames(V), names(cf))
  expect_true(isSymmetric(V))
  expect_true(all(diag(V) > 0))
  expect_equal(rownames(vcov(fit, "zero")), names(coef(fit, "zero")))
  # rho by the delta method agrees with TMB's ADREPORT
  expect_equal(sqrt(V["rho", "rho"]), fit$se_natural[["rho"]],
               tolerance = 1e-4)
})

test_that("vcov() of a fit with rho fixed has NA for rho only", {
  V <- vcov(fits_rho0()$nested)
  expect_true(all(is.na(V["rho", ])))
  expect_true(all(is.finite(V[-6, -6])))
})

test_that("logLik(), AIC() and nobs() work", {
  fit <- fit_small()
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_equal(attr(ll, "df"), 8)
  expect_equal(nobs(fit), 225)
  expect_equal(AIC(fit), -2 * as.numeric(ll) + 2 * 8)
  expect_equal(BIC(fit), -2 * as.numeric(ll) + log(225) * 8)
})

test_that("predict() returns the fitted quantities", {
  fit <- fit_small()
  expect_equal(predict(fit), fitted(fit))
  expect_equal(predict(fit, type = "count"), fit$mu)
  expect_equal(predict(fit, type = "zero"), fit$pi)
  expect_equal(predict(fit, type = "eta"), fit$eta)
  expect_true(all(predict(fit, type = "zero") > 0 &
                    predict(fit, type = "zero") < 1))
})

test_that("predict() with newdata gives counterfactual scenarios", {
  sim <- sim_small()
  fit <- fit_small()
  # Unchanged covariates reproduce the fitted values
  expect_equal(predict(fit, newdata = sim$data), fitted(fit),
               tolerance = 1e-8)
  # A change in unit 1 moves unit 1 and its neighbours, not distant units
  nd <- sim$data
  nd$x1[1] <- nd$x1[1] + 1
  delta <- predict(fit, newdata = nd) - fitted(fit)
  nbrs <- which(sim$W[1, ] > 0)
  expect_gt(delta[1], 0)
  expect_true(all(delta[nbrs] > 0))
  expect_lt(abs(delta[225]), abs(delta[nbrs[1]]))
  # Inflation covariates act on pi only
  nd <- sim$data
  nd$z1 <- nd$z1 + 1
  expect_equal(predict(fit, newdata = nd, type = "count"), fit$mu,
               tolerance = 1e-8)
  expect_true(all(predict(fit, newdata = nd, type = "zero") > fit$pi))

  expect_error(predict(fit, newdata = sim$data[-1, ]), "one row per")
  X <- cbind(1, as.matrix(sim$data[, c("x1", "x2")]))
  fit_m <- zisar_fit(sim$data$y, X, NULL, sim$W)
  expect_error(predict(fit_m, newdata = sim$data), "formula interface")
})

test_that("residuals() of both types", {
  fit <- fit_small()
  expect_equal(residuals(fit, "response"), fit$y - fitted(fit))
  r <- residuals(fit)
  expect_length(r, 225)
  expect_true(all(is.finite(r)))
  expect_lt(abs(mean(r)), 0.5)
})

test_that("print() and summary() show the model", {
  fit <- fit_small()
  expect_output(print(fit), "SAR negative binomial")
  s <- summary(fit)
  expect_s3_class(s, "summary.zisar")
  expect_equal(dim(s$count), c(3, 4))
  expect_equal(rownames(s$spatial), c("rho", "kappa", "sigma"))
  expect_output(print(s), "Spatial and dispersion parameters")
  expect_output(print(s), "converged")

  expect_output(print(summary(fits_rho0()$nested)), "Held fixed: rho")
})

test_that("simulate() draws response vectors from the fit", {
  fit <- fit_small()
  s1 <- simulate(fit, nsim = 3, seed = 10)
  expect_s3_class(s1, "data.frame")
  expect_equal(dim(s1), c(225, 3))
  expect_true(all(as.matrix(s1) >= 0))
  expect_equal(simulate(fit, nsim = 3, seed = 10), s1)
})
