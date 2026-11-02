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

# Data from the general model (inflation field with rho0 = 0.6, sigma0 = 1):
# list(sim, restricted fit, general fit)
fits_general <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(2)
      sim <- rzisar(15, beta = c(1, 0.5), rho = 0.4, gamma = c(-1, 0.8),
                    rho0 = 0.6, sigma0 = 1)
      cache <<- list(
        sim = sim,
        restricted = zisar(y ~ x1 | z1, data = sim$data, W = sim$W),
        general = zisar(y ~ x1 | z1, data = sim$data, W = sim$W,
                        zero_spatial = TRUE)
      )
    }
    cache
  }
})
