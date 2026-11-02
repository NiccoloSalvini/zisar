test_that("the general model reports the inflation field", {
  f <- fits_general()
  fit <- f$general
  expect_true(fit$zero_spatial)
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
  expect_output(print(s), "inflation field")
  expect_output(print(fit), "rho0")
})

test_that("the restricted model is the general one with sigma0 near 0", {
  f <- fits_general()
  nested <- zisar(y ~ x1 | z1, data = f$sim$data, W = f$sim$W,
                  zero_spatial = TRUE,
                  fixed = list(rho0 = 0, sigma0 = 0.001))
  expect_equal(nested$loglik, f$restricted$loglik, tolerance = 1e-4)
  expect_equal(nested$rho, f$restricted$rho, tolerance = 1e-3)
  expect_equal(nested$df, f$restricted$df)
  expect_setequal(summary(nested)$fixed, c("rho0", "sigma0"))
  expect_gte(f$general$loglik, f$restricted$loglik - 1e-6)
})

test_that("`fixed` matches parameter names exactly", {
  fx <- zisar:::.internal_fixed(list(rho0 = 0.5, sigma0 = 2), 2, 2, TRUE)
  expect_named(fx, c("rho0_raw", "log_sigma0"))
  fx <- zisar:::.internal_fixed(list(rho = 0.5), 2, 2, TRUE)
  expect_named(fx, "rho_raw")
})

test_that("rho0 and sigma0 need the general model", {
  sim <- sim_small()
  expect_error(zisar(y ~ x1 | z1, data = sim$data, W = sim$W,
                     fixed = list(sigma0 = 0.5)), "zero_spatial = TRUE")
  expect_error(zisar(y ~ x1 | z1, data = sim$data, W = sim$W,
                     zero_spatial = NA), "TRUE or FALSE")
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
