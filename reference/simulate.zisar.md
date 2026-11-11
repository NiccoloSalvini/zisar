# Simulate responses from a fitted ZISAR model

Draws new count vectors from the fitted model, for example for a
parametric bootstrap. Each draw takes new random effects \\u\\ (and a
new inflation field, \\g\\ or \\\alpha\\, when the model has one), new
structural-zero indicators and new negative binomial counts, at the
estimated parameters and the observed covariates.

## Usage

``` r
# S3 method for class 'zisar'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  A fitted `"zisar"` object.

- nsim:

  Number of response vectors to simulate.

- seed:

  Optional seed, passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

- ...:

  Not used.

## Value

A data frame with `nsim` columns, one simulated response vector each.

## Details

To check a large gradient at the optimum (see
[`zisar()`](https://www.niccolosalvini.me/packages/zisar/reference/zisar.md)),
refit the model to a few simulated responses: if the refits show the
same gradient while recovering the parameters, it belongs to the model
and not to a spurious optimum. For bootstrap standard errors, start each
refit from the estimates with `start`.

## Examples

``` r
set.seed(1)
sim <- rzisar(10, beta = c(1, 0.5), rho = 0.4, gamma = c(-1, 0.8))
fit <- zisar(y ~ x1 | z1, data = sim$data, W = sim$W)
ys <- simulate(fit, nsim = 100, seed = 2)
# Share of zeros: observed against simulated
c(observed = mean(sim$data$y == 0), simulated = mean(as.matrix(ys) == 0))
#>  observed simulated 
#>    0.3500    0.3528 
```
