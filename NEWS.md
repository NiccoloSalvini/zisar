# zisar 0.1.0

First release.

* `zisar()` fits the ZISAR model from a two-part formula
  `y ~ count terms | zero-inflation terms` and a spatial weights matrix
  (base matrix, sparse `Matrix` or `spdep` `listw`); `zisar_fit()` is the
  matrix interface.
* Estimation by Laplace approximation in 'TMB', with the log-determinant from
  the eigenvalues of `W` or, for large `n`, from a sparse Cholesky factor.
* `fixed` holds parameters at given values, for the nested ZINB-RE
  (`rho = 0`) and SAR-NB-RE models.
* `zero_spatial = TRUE` fits the general model, with a SAR random field in
  the zero-inflation predictor (parameters `rho0` and `sigma0`); `rzisar()`
  simulates from it with `sigma0 > 0`.
* `offset()` terms in the count part, for exposures.
* Methods: `summary()`, `coef()`, `vcov()`, `logLik()` (hence `AIC()`,
  `BIC()`), `nobs()`, `fitted()`, `residuals()`, `predict()` (including
  counterfactual scenarios on the same units) and `simulate()`.
* `zisar_effects()`: average direct, indirect and total effects with
  simulation standard errors, including the zero-inflation channel.
* `rzisar()` simulates from the model; `lattice_weights()` builds lattice
  contiguity weights.
