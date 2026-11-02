test_that("rzisar() returns data ready for zisar()", {
  set.seed(1)
  sim <- rzisar(10, beta = c(1, 0.5, -0.3), rho = 0.4, gamma = c(-1, 0.8))
  expect_named(sim$data, c("y", "x1", "x2", "z1"))
  expect_equal(nrow(sim$data), 100)
  expect_equal(dim(sim$W), c(100, 100))
  expect_equal(nrow(sim$coords), 100)
  expect_true(all(sim$data$y >= 0))
  expect_true(all(sim$data$y[sim$s == 1] == 0))
  expect_equal(sim$truth$rho, 0.4)
})

test_that("rzisar() reproduces the latent field it reports", {
  set.seed(2)
  sim <- rzisar(8, beta = c(0.5, 1), rho = 0.6, gamma = -2, sigma = 0)
  X <- cbind(1, sim$data$x1)
  A <- diag(64) - 0.6 * as.matrix(sim$W)
  expect_equal(sim$eta, as.numeric(solve(A, X %*% c(0.5, 1))))
  expect_equal(sim$mu, exp(sim$eta))
  expect_equal(sim$pi, rep(plogis(-2), 64))
})

test_that("rzisar() accepts a user weights matrix", {
  set.seed(3)
  W <- lattice_weights(5, 6, type = "rook")
  sim <- rzisar(beta = 1, rho = 0.2, gamma = -1, W = W)
  expect_equal(nrow(sim$data), 30)
  expect_null(sim$coords)
  expect_named(sim$data, "y")
})

test_that("rzisar() draws an inflation field in the general model", {
  set.seed(4)
  a <- rzisar(8, beta = c(1, 0.5), rho = 0.3, gamma = c(-1, 0.5))
  set.seed(4)
  b <- rzisar(8, beta = c(1, 0.5), rho = 0.3, gamma = c(-1, 0.5),
              rho0 = 0.5, sigma0 = 0)
  expect_identical(a$data, b$data)
  expect_equal(a$g, rep(0, 64))

  set.seed(5)
  sim <- rzisar(8, beta = 1, rho = 0.3, gamma = -1, rho0 = 0.5, sigma0 = 1)
  expect_gt(sd(sim$g), 0)
  expect_equal(sim$pi, plogis(-1 + sim$g))
  expect_equal(sim$truth$sigma0, 1)
})
