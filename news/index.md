# Changelog

## zisar 0.2.0

- `zero_spatial` now selects the inflation predictor: `"none"`
  (restricted model), `"reduced"`, `"field"` or `"sar"`. `TRUE` and
  `FALSE` still work and mean `"field"` and `"none"`.
- New reduced-rank inflation field (`zero_spatial = "reduced"`): a SAR
  field on the knots of a lattice, spread over the units by a bisquare
  basis built by the new
  [`zisar_basis()`](https://www.niccolosalvini.me/packages/zisar/reference/zisar_basis.md).
  It is the recommended general model: the Laplace approximation of the
  unit-level field (`"field"`) fails, with likelihood-ratio statistics
  in the thousands on data without a field.
- New SAR filter on the inflation predictor (`zero_spatial = "sar"`),
  `logit(pi) = (I - rho0 W)^{-1} Z gamma`, with spillovers of the
  inflation covariates in
  [`zisar_effects()`](https://www.niccolosalvini.me/packages/zisar/reference/zisar_effects.md)
  and [`predict()`](https://rdrr.io/r/stats/predict.html).
- Fits with a random inflation field start from the restricted fit and a
  grid of `(rho0, sigma0)`, returned in `field_profile`.
- The warning on a large gradient no longer calls the optimum spurious:
  with structural zeros among units of high intensity a genuine optimum
  has one too. Check it on data from
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) and use a
  parametric bootstrap for standard errors (see
  [`vignette("zisar")`](https://www.niccolosalvini.me/packages/zisar/articles/zisar.md)).
  A positive log-likelihood is still reported as spurious.
- New argument `restarts = FALSE` skips the restarts of a fit that fails
  the checks, for slow fits where a large gradient is expected.
- [`rzisar()`](https://www.niccolosalvini.me/packages/zisar/reference/rzisar.md)
  gains `zero_spatial`; `rho0` and `sigma0` now need it.
- The README and
  [`vignette("zisar")`](https://www.niccolosalvini.me/packages/zisar/articles/zisar.md)
  explain the model with two figures: the two kinds of zeros and the
  reach of the spatial multiplier (`data-raw/figures.R`).

## zisar 0.1.0

First release.

- [`zisar()`](https://www.niccolosalvini.me/packages/zisar/reference/zisar.md)
  fits the ZISAR model from a two-part formula
  `y ~ count terms | zero-inflation terms` and a spatial weights matrix
  (base matrix, sparse `Matrix` or `spdep` `listw`);
  [`zisar_fit()`](https://www.niccolosalvini.me/packages/zisar/reference/zisar.md)
  is the matrix interface.
- Estimation by Laplace approximation in ‘TMB’, with the log-determinant
  from the eigenvalues of `W` or, for large `n`, from a sparse Cholesky
  factor.
- `fixed` holds parameters at given values, for the nested ZINB-RE
  (`rho = 0`) and SAR-NB-RE models.
- `zero_spatial = TRUE` fits the general model, with a SAR random field
  in the zero-inflation predictor (parameters `rho0` and `sigma0`);
  [`rzisar()`](https://www.niccolosalvini.me/packages/zisar/reference/rzisar.md)
  simulates from it with `sigma0 > 0`.
- [`offset()`](https://rdrr.io/r/stats/offset.html) terms in the count
  part, for exposures.
- Methods: [`summary()`](https://rdrr.io/r/base/summary.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) (hence
  [`AIC()`](https://rdrr.io/r/stats/AIC.html),
  [`BIC()`](https://rdrr.io/r/stats/AIC.html)),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html) (including
  counterfactual scenarios on the same units) and
  [`simulate()`](https://rdrr.io/r/stats/simulate.html).
- [`zisar_effects()`](https://www.niccolosalvini.me/packages/zisar/reference/zisar_effects.md):
  average direct, indirect and total effects with simulation standard
  errors, including the zero-inflation channel.
- [`rzisar()`](https://www.niccolosalvini.me/packages/zisar/reference/rzisar.md)
  simulates from the model;
  [`lattice_weights()`](https://www.niccolosalvini.me/packages/zisar/reference/lattice_weights.md)
  builds lattice contiguity weights.
