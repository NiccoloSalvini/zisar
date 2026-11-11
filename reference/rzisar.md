# Simulate data from the ZISAR model

Draws one data set from the ZISAR model, on a regular lattice with queen
contiguity (the Monte Carlo design of the ZISAR paper) or on a weights
matrix supplied by the user. Covariates are independent standard normal;
the first element of `beta` and of `gamma` is the intercept.
`zero_spatial` adds spatial structure to the inflation predictor: a SAR
field with one element per unit (`"field"`, the general model of the
paper) or a SAR filter on \\Z\gamma\\ (`"sar"`); see
[`zisar()`](https://www.niccolosalvini.me/packages/zisar/reference/zisar.md).

## Usage

``` r
rzisar(
  n_side,
  beta,
  rho,
  gamma,
  kappa = 2,
  sigma = 0.5,
  zero_spatial = c("none", "field", "sar"),
  rho0 = 0,
  sigma0 = 0,
  W = NULL
)
```

## Arguments

- n_side:

  Side of the square lattice, so \\n\\ = `n_side^2` units. Ignored when
  `W` is supplied.

- beta:

  Count-process coefficients, intercept first.

- rho:

  Spatial autoregressive parameter, in \\(-1, 1)\\.

- gamma:

  Zero-inflation coefficients, intercept first.

- kappa:

  Negative binomial shape.

- sigma:

  Standard deviation of the random effects.

- zero_spatial:

  Inflation predictor: `"none"` (restricted model), `"field"` or
  `"sar"`.

- rho0, sigma0:

  Spatial autoregressive parameter and standard deviation of the
  inflation field (`"field"`), or parameter of the filter (`"sar"`,
  which takes no `sigma0`).

- W:

  Optional row-standardised weights matrix (or `spdep` `listw`) to use
  instead of the lattice.

## Value

A list with

- data:

  Data frame with the counts `y`, the count covariates `x1, x2, ...` and
  the inflation covariates `z1, z2, ...`, ready for
  `zisar(y ~ x1 + ... | z1 + ..., data = data, W = W)`.

- W:

  The weights matrix (sparse).

- coords:

  Lattice coordinates (`row`, `col`), for example for
  [`zisar_basis()`](https://www.niccolosalvini.me/packages/zisar/reference/zisar_basis.md);
  `NULL` with a user-supplied `W`.

- truth:

  The parameters used.

- eta, g, zeta, mu, pi, s:

  Latent field, inflation field (zero without one), inflation predictor,
  negative binomial means, inflation probabilities and structural-zero
  indicators.

## Examples

``` r
set.seed(42)
sim <- rzisar(15, beta = c(1, 0.5, -0.3), rho = 0.5, gamma = c(-1, 0.7))
head(sim$data)
#>   y         x1         x2         z1
#> 1 6  1.3709584 -0.1740182 -1.3038212
#> 2 0 -0.5646982  0.5156677 -0.2509145
#> 3 0  0.3631284 -0.2343653  0.1710074
#> 4 0  0.6328626 -0.6585034 -0.4034675
#> 5 7  0.4042683  1.2502366  0.1046594
#> 6 4 -0.1061245 -0.2717637 -0.3188808
mean(sim$data$y == 0)
#> [1] 0.4622222

# Structural zeros clustered by a field of their own
sim2 <- rzisar(15, beta = c(1, 0.5), rho = 0.5, gamma = c(-1, 0.7),
               zero_spatial = "field", rho0 = 0.8, sigma0 = 1)
cor(sim2$g, as.numeric(sim2$W %*% sim2$g))
#> [1] 0.7777463
```
