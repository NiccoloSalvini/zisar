test_that("effects match the dense matrix of partial derivatives", {
  fit <- fit_small()
  eff <- zisar_effects(fit, B = 0)
  expect_s3_class(eff, "zisar_effects")
  expect_equal(unique(eff$term), c("x1", "x2", "z1"))

  # Brute force: dE[y_i]/dx_jr = d_i S_ij beta_r, plus the inflation channel
  S <- solve(diag(fit$n) - fit$rho * as.matrix(fit$W))
  d <- (1 - fit$pi) * fit$mu
  M <- d * S
  for (r in c("x1", "x2")) {
    b <- fit$coefficients$count[[r]]
    e <- eff[eff$term == r, ]
    expect_equal(e$estimate[e$effect == "direct"], mean(diag(M)) * b)
    expect_equal(e$estimate[e$effect == "total"], mean(rowSums(M)) * b)
  }
  zch <- mean(fit$mu * fit$pi * (1 - fit$pi))
  ez <- eff[eff$term == "z1", ]
  expect_equal(ez$estimate[ez$effect == "direct"],
               -zch * fit$coefficients$zero[["z1"]])
  expect_equal(ez$estimate[ez$effect == "indirect"], 0)
})

test_that("direct + indirect = total, and no indirect effect when rho = 0", {
  eff <- zisar_effects(fit_small(), B = 0)
  w <- reshape(eff, idvar = "term", timevar = "effect", direction = "wide")
  expect_equal(w$estimate.direct + w$estimate.indirect, w$estimate.total)

  eff0 <- zisar_effects(fits_rho0()$nested, B = 20)
  expect_equal(eff0$estimate[eff0$effect == "indirect"], c(0, 0, 0))
})

test_that("a covariate in both parts adds the inflation channel", {
  sim <- sim_small()
  fit <- zisar(y ~ x1 + x2 | x1, data = sim$data, W = sim$W)
  eff <- zisar_effects(fit, B = 0)
  expect_equal(unique(eff$term), c("x1", "x2"))
  S <- solve(diag(fit$n) - fit$rho * as.matrix(fit$W))
  d <- (1 - fit$pi) * fit$mu
  zch <- mean(fit$mu * fit$pi * (1 - fit$pi))
  direct <- mean(d * diag(S)) * fit$coefficients$count[["x1"]] -
    zch * fit$coefficients$zero[["x1"]]
  expect_equal(eff$estimate[eff$term == "x1" & eff$effect == "direct"],
               direct)
})

test_that("simulation standard errors and intervals", {
  set.seed(3)
  eff <- zisar_effects(fit_small(), B = 100, level = 0.9)
  expect_true(all(c("se", "lower", "upper") %in% names(eff)))
  expect_true(all(eff$se[eff$term != "z1" | eff$effect != "indirect"] > 0))
  expect_true(all(eff$lower <= eff$upper))
  expect_equal(attr(eff, "B"), 100)
  expect_output(print(eff), "B = 100")
})
