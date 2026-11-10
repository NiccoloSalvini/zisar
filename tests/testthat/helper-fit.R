# A small simulated data set and its fit, shared by the tests
sim_small <- function(seed = 1, n_side = 15) {
  set.seed(seed)
  rzisar(n_side, beta = c(1, 0.5, -0.3), rho = 0.4, gamma = c(-1, 0.8))
}

fit_small <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      sim <- sim_small()
      cache <<- zisar(y ~ x1 + x2 | z1, data = sim$data, W = sim$W)
    }
    cache
  }
})

# Data without spatial dependence, where the nested rho = 0 model is well
# specified: list(sim, full fit, fit with rho held at 0)
fits_rho0 <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(11)
      sim <- rzisar(15, beta = c(1, 0.5, -0.3), rho = 0, gamma = c(-1, 0.8))
      cache <<- list(
        sim = sim,
        full = zisar(y ~ x1 + x2 | z1, data = sim$data, W = sim$W),
        nested = zisar(y ~ x1 + x2 | z1, data = sim$data, W = sim$W,
                       fixed = list(rho = 0))
      )
    }
    cache
  }
})

# Data from the general model (unit-level inflation field with rho0 = 0.6,
# sigma0 = 1): list(sim, restricted fit, fit with the unit-level field)
fits_general <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(2)
      sim <- rzisar(15, beta = c(1, 0.5), rho = 0.4, gamma = c(-1, 0.8),
                    zero_spatial = "field", rho0 = 0.6, sigma0 = 1)
      cache <<- list(
        sim = sim,
        restricted = zisar(y ~ x1 | z1, data = sim$data, W = sim$W),
        general = zisar(y ~ x1 | z1, data = sim$data, W = sim$W,
                        zero_spatial = "field")
      )
    }
    cache
  }
})

# Data with a strong inflation field on a 20 x 20 lattice, fitted with the
# reduced-rank field on knots 4 cells apart: list(sim, basis, fit)
fits_reduced <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(3)
      sim <- rzisar(20, beta = c(1, 0.5), rho = 0.4, gamma = c(-0.5, 0.8),
                    zero_spatial = "field", rho0 = 0.9, sigma0 = 0.6)
      b <- zisar_basis(sim$coords, h = 4)
      cache <<- list(
        sim = sim, basis = b,
        fit = zisar(y ~ x1 | z1, data = sim$data, W = sim$W,
                    zero_spatial = "reduced", basis = b)
      )
    }
    cache
  }
})

# Data with a SAR filter on the inflation predictor (rho0 = 0.6):
# list(sim, fit)
fits_sar <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(7)
      sim <- rzisar(15, beta = c(1, 0.5), rho = 0.4, gamma = c(-0.5, 1),
                    zero_spatial = "sar", rho0 = 0.6)
      cache <<- list(
        sim = sim,
        fit = zisar(y ~ x1 | z1, data = sim$data, W = sim$W,
                    zero_spatial = "sar")
      )
    }
    cache
  }
})
