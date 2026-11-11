# Predictions from a fitted ZISAR model

Fitted or counterfactual means, inflation probabilities and latent field
of a ZISAR model.

## Usage

``` r
# S3 method for class 'zisar'
predict(
  object,
  newdata = NULL,
  type = c("response", "count", "zero", "eta"),
  ...
)
```

## Arguments

- object:

  A fitted `"zisar"` object.

- newdata:

  Optional data frame of covariates for the same units, in the order of
  `W`; needs a model fitted with
  [`zisar()`](https://www.niccolosalvini.me/packages/zisar/reference/zisar.md).

- type:

  `"response"` for the mean \\E\[y_i\] = (1 - \pi_i)\mu_i\\, `"count"`
  for the mean of the negative binomial component \\\mu_i\\, `"zero"`
  for the probability of a structural zero \\\pi_i\\, or `"eta"` for the
  latent field \\\eta_i\\ (without offset).

- ...:

  Not used.

## Value

A numeric vector with one value per spatial unit.

## Details

Without `newdata`, predictions are at the fitted values, with the latent
field \\\eta\\ at its mode given the data.

With `newdata`, predictions are a counterfactual scenario for the *same*
spatial units: `newdata` holds one row per unit, in the order of the
weights matrix, with some covariates changed. The latent field is
recomputed through the spatial multiplier, \\\eta^\* = (I - \rho
W)^{-1}(X^\*\beta + \hat u)\\, keeping the predicted random effects
\\\hat u = (I - \rho W)\hat\eta - X\beta\\ at their estimates, so a
change in one unit spills over to its neighbours. An inflation field
\\g\\ is kept at its mode; with the SAR filter (`zero_spatial = "sar"`)
the inflation predictor is recomputed through \\(I - \rho_0 W)^{-1}\\,
so changes in the inflation covariates spill over too. Predictions for
new units would need a new weights matrix and are not supported.

## Examples

``` r
set.seed(1)
sim <- rzisar(10, beta = c(1, 0.5), rho = 0.4, gamma = c(-1, 0.8))
fit <- zisar(y ~ x1 | z1, data = sim$data, W = sim$W)
head(predict(fit, type = "zero"))
#> [1] 0.1924763 0.2642031 0.1660733 0.2783742 0.1892080 0.5107059

# Raise x1 by one in unit 1 only: the change spills over to its neighbours
nd <- sim$data
nd$x1[1] <- nd$x1[1] + 1
delta <- predict(fit, newdata = nd) - predict(fit)
round(head(delta), 3)
#> [1] 2.413 0.087 0.008 0.001 0.000 0.000
```
