test_that("the unit-level field is reported", {
  f <- fits_general()
  fit <- f$general
  expect_equal(fit$zero_spatial, "field")
  expect_true(is.list(fit$field_profile))
  expect_equal(nrow(fit$field_profile$grid), 25)
  expect_equal(fit$field_profile$restricted_loglik, f$restricted$loglik,
               tolerance = 1e-6)
  expect_true(fit$laplace_ok)
  expect_equal(fit$df, f$restricted$df + 2)
  expect_true(abs(fit$rho0) < 1)
  expect_gt(fit$sigma0, 0)
  expect_length(fit$g, nrow(f$sim$data))
  expect_equal(fit$pi, plogis(as.numeric(fit$Z %*% fit$coefficients$zero) +
                                fit$g))
  expect_gt(cor(fit$g, f$sim$g), 0)
  expect_named(fit$se_natural, c("rho", "kappa", "sigma", "rho0", "sigma0"))
  expect_null(f$restricted$rho0)
  expect_equal(f$restricted$g, rep(0, f$restricted$n))
})

test_that("coef(), vcov() and summary() include rho0 and sigma0", {
  fit <- fits_general()$general
  cf <- coef(fit)
  expect_equal(tail(names(cf), 2), c("rho", "rho0"))
  V <- vcov(fit)
  expect_equal(dimnames(V), list(names(cf), names(cf)))
  expect_equal(sqrt(V["rho0", "rho0"]), fit$se_natural[["rho0"]],
               tolerance = 1e-6)
  s <- summary(fit)
  expect_equal(rownames(s$spatial), c("rho", "kappa", "sigma", "rho0",
                                      "sigma0"))
  expect_true(is.na(s$spatial["sigma0", "z value"]))
  expect_output(print(s), "unit-level inflation field")
  expect_output(print(fit), "rho0")
})

test_that("the restricted model is the general one with sigma0 near 0", {
  f <- fits_general()
  nested <- zisar(y ~ x1 | z1, data = f$sim$data, W = f$sim$W,
                  zero_spatial = TRUE,
                  fixed = list(rho0 = 0, sigma0 = 0.001))
  expect_equal(nested$zero_spatial, "field")
  expect_equal(nested$loglik, f$restricted$loglik, tolerance = 1e-4)
  expect_equal(nested$rho, f$restricted$rho, tolerance = 1e-3)
  expect_equal(nested$df, f$restricted$df)
  expect_setequal(summary(nested)$fixed, c("rho0", "sigma0"))
  expect_gte(f$general$loglik, f$restricted$loglik - 1e-6)
})

test_that("`fixed` matches parameter names exactly", {
  fx <- zisar:::.internal_fixed(list(rho0 = 0.5, sigma0 = 2), 2, 2, "field")
  expect_named(fx, c("rho0_raw", "log_sigma0"))
  fx <- zisar:::.internal_fixed(list(rho = 0.5), 2, 2, "field")
  expect_named(fx, "rho_raw")
})

test_that("rho0 and sigma0 need a spatial inflation predictor", {
  sim <- sim_small()
  expect_error(zisar(y ~ x1 | z1, data = sim$data, W = sim$W,
                     fixed = list(sigma0 = 0.5)), "zero_spatial")
  expect_error(zisar(y ~ x1 | z1, data = sim$data, W = sim$W,
                     zero_spatial = "sar", fixed = list(sigma0 = 0.5)),
               "random inflation field")
  expect_error(zisar(y ~ x1 | z1, data = sim$data, W = sim$W,
                     zero_spatial = NA), "must be one of")
  expect_error(zisar(y ~ x1 | z1, data = sim$data, W = sim$W,
                     zero_spatial = "reduced"), "zisar_basis")
  b <- zisar_basis(expand.grid(1:5, 1:5), h = 2)
  expect_error(zisar(y ~ x1 | z1, data = sim$data, W = sim$W,
                     zero_spatial = "reduced", basis = b), "25 units")
  expect_error(zisar:::.internal_fixed(list(rho0 = -0.2), 2, 2, "reduced"),
               "\\[0, 1\\)")
})

test_that("predict(), effects and simulate() use the inflation field", {
  f <- fits_general()
  fit <- f$general
  expect_equal(predict(fit, newdata = f$sim$data), fitted(fit),
               tolerance = 1e-6)
  eff <- zisar_effects(fit, B = 0)
  expect_equal(eff$estimate[eff$term == "z1" & eff$effect == "direct"],
               -mean(fit$mu * fit$pi * (1 - fit$pi)) *
                 fit$coefficients$zero[["z1"]])
  ys <- simulate(fit, nsim = 2, seed = 1)
  expect_equal(dim(ys), c(fit$n, 2))
})

test_that("the reduced-rank field recovers the clustering of the zeros", {
  f <- fits_reduced()
  fit <- f$fit
  expect_equal(fit$zero_spatial, "reduced")
  expect_length(fit$alpha, f$basis$K)
  expect_equal(fit$g, as.numeric(f$basis$B %*% fit$alpha))
  expect_equal(fit$pi, plogis(as.numeric(fit$Z %*% fit$coefficients$zero) +
                                fit$g))
  expect_gte(fit$rho0, 0)
  expect_gt(fit$sigma0, 0.1)
  expect_gt(cor(fit$g, f$sim$g), 0.5)
  expect_gt(max(fit$field_profile$grid$gain), 0)
  expect_named(fit$se_natural, c("rho", "kappa", "sigma", "rho0", "sigma0"))
  expect_output(print(summary(fit)), "reduced-rank inflation field")
  expect_output(print(fit), "knots")
  ys <- simulate(fit, nsim = 2, seed = 1)
  expect_equal(dim(ys), c(fit$n, 2))
  expect_equal(predict(fit, newdata = f$sim$data), fitted(fit),
               tolerance = 1e-6)
})

test_that("the SAR filter on the inflation predictor is estimated", {
  f <- fits_sar()
  fit <- f$fit
  expect_equal(fit$zero_spatial, "sar")
  expect_true(fit$laplace_ok)
  expect_null(fit$sigma0)
  expect_null(fit$field_profile)
  expect_equal(fit$g, rep(0, fit$n))
  expect_equal(fit$rho0, 0.6, tolerance = 0.4)
  A <- diag(fit$n) - fit$rho0 * as.matrix(fit$W)
  zeta <- as.numeric(solve(A, fit$Z %*% fit$coefficients$zero))
  expect_equal(fit$zeta, zeta, tolerance = 1e-8)
  expect_equal(fit$pi, plogis(zeta))
  expect_equal(fit$df, 8)
  expect_named(fit$se_natural, c("rho", "kappa", "sigma", "rho0"))
  expect_equal(tail(names(coef(fit)), 2), c("rho", "rho0"))
  expect_equal(rownames(summary(fit)$spatial),
               c("rho", "kappa", "sigma", "rho0"))
  expect_output(print(summary(fit)), "SAR filter")

  # The filter spreads a change in z1 to the neighbours
  nd <- f$sim$data
  nd$z1[1] <- nd$z1[1] + 1
  dpi <- predict(fit, newdata = nd, type = "zero") - predict(fit, type = "zero")
  expect_gt(sum(abs(dpi[-1]) > 1e-8), 0)

  # Inflation channel of z1: direct from diag(S0), total from S0 1
  eff <- zisar_effects(fit, B = 0)
  S0 <- solve(A)
  zc <- fit$mu * fit$pi * (1 - fit$pi)
  g1 <- fit$coefficients$zero[["z1"]]
  est <- function(e) eff$estimate[eff$term == "z1" & eff$effect == e]
  expect_equal(est("direct"), -mean(zc * diag(S0)) * g1)
  expect_equal(est("total"), -mean(zc * rowSums(S0)) * g1)
  expect_lt(abs(est("indirect")), abs(est("total")))
  expect_true(all(is.finite(zisar_effects(fit, B = 20)$se)))

  # rho0 = 0 is the restricted model
  r0 <- zisar(y ~ x1 | z1, data = f$sim$data, W = f$sim$W,
              zero_spatial = "sar", fixed = list(rho0 = 0))
  r <- zisar(y ~ x1 | z1, data = f$sim$data, W = f$sim$W)
  expect_equal(r0$loglik, r$loglik, tolerance = 1e-5)
})
