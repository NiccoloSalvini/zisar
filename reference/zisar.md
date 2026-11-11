# Fit a zero-inflated spatial autoregressive (ZISAR) model

Fits the ZISAR model of Salvini, Arbia and Nardelli (2026): a
zero-inflated negative binomial model whose log-mean follows a
simultaneous autoregressive (SAR) process with unit-level random
effects. By default the inflation probability depends on covariates only
(the restricted model); `zero_spatial` adds spatial structure to it (the
general model). `zisar()` takes a two-part formula and a data frame;
`zisar_fit()` is the matrix interface underneath it.

## Usage

``` r
zisar(
  formula,
  data,
  W,
  zero_spatial = c("none", "reduced", "field", "sar"),
  basis = NULL,
  fixed = NULL,
  start = NULL,
  logdet = c("auto", "eigen", "sparse"),
  restarts = TRUE,
  control = list(),
  silent = TRUE
)

zisar_fit(
  y,
  X,
  Z = NULL,
  W,
  offset = NULL,
  zero_spatial = c("none", "reduced", "field", "sar"),
  basis = NULL,
  fixed = NULL,
  start = NULL,
  logdet = c("auto", "eigen", "sparse"),
  restarts = TRUE,
  control = list(),
  silent = TRUE
)
```

## Arguments

- formula:

  A two-part formula `y ~ x1 + x2 | z1 + z2`: the terms before `|` enter
  the count process, those after it the zero-inflation process. Without
  `|`, both processes use the same terms. Use `| 1` for a constant
  inflation probability. An
  [`offset()`](https://rdrr.io/r/stats/offset.html) term in the count
  part enters the log-mean with coefficient one.

- data:

  A data frame with one row per spatial unit, in the order of the rows
  of `W`. Missing values are not allowed: every unit enters the spatial
  multiplier, so none can be dropped.

- W:

  Spatial weights: an \\n \times n\\ matrix (base or sparse 'Matrix') or
  an `spdep` `listw` object, with zero diagonal and preferably
  row-standardised.

- zero_spatial:

  Inflation predictor: `"none"` (restricted model), `"reduced"`,
  `"field"` or `"sar"`; see Model. `FALSE` and `TRUE` stand for `"none"`
  and `"field"`.

- basis:

  For `zero_spatial = "reduced"`, the basis from
  [`zisar_basis()`](https://www.niccolosalvini.me/packages/zisar/reference/zisar_basis.md).

- fixed:

  Optional named list of parameters held at fixed values: `rho`,
  `kappa`, `sigma`, `rho0` and `sigma0` (scalars, natural scale), `beta`
  (length `ncol(X)`) or `gamma` (length `ncol(Z)`). In `beta` and
  `gamma`, `NA` entries are estimated.

- start:

  Optional named list of starting values, with the same names and scales
  as `fixed`. Defaults: `beta` from a Poisson regression, inflation
  intercept from the share of zeros, `rho = 0`, `kappa = 1`,
  `sigma = 0.5`, `rho0 = 0` and, for the fields, `rho0` and `sigma0`
  from the grid described in Estimation. To restart from a previous fit,
  use
  `start = list(beta = coef(fit, "count"), gamma = coef(fit, "zero"), rho = fit$rho, kappa = fit$kappa, sigma = fit$sigma)`.

- logdet:

  Method for the log-determinant: `"auto"`, `"eigen"` or `"sparse"`. See
  Estimation.

- restarts:

  If `TRUE`, restart a fit that fails the checks described in
  Estimation. Set it to `FALSE` when a single fit is slow and a large
  gradient is expected.

- control:

  List of control parameters for
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html), merged with
  the defaults `eval.max = 1000`, `iter.max = 500`, `rel.tol = 1e-10`.

- silent:

  If `FALSE`, print the 'TMB' tracing output.

- y:

  Vector of non-negative integer counts.

- X:

  Design matrix of the count process, \\n \times p\\; include a column
  of ones for an intercept.

- Z:

  Design matrix of the zero-inflation process, \\n \times q\\. `NULL`
  gives an intercept only.

- offset:

  Optional vector of length \\n\\ added to the count log-mean (for
  example the log of an exposure).

## Value

An object of class `"zisar"`, a list with, among others:

- coefficients:

  List with the `count` (\\\beta\\) and `zero` (\\\gamma\\)
  coefficients.

- rho, kappa, sigma:

  Spatial autoregressive parameter, negative binomial shape and
  random-effect standard deviation.

- zero_spatial:

  The inflation predictor, as in the argument.

- rho0, sigma0:

  Parameters of the inflation field or filter; `NULL` when absent.

- eta:

  Mode of the latent field \\\eta\\ given the data.

- g:

  Mode of the inflation field (\\g\\, or \\B\alpha\\ for the
  reduced-rank field); zero without a field.

- zeta:

  The inflation predictor, \\\mathrm{logit}(\pi)\\.

- fitted.values:

  Fitted means \\E\[y_i\] = (1 - \pi_i)\mu_i\\.

- loglik:

  Laplace-approximated marginal log-likelihood.

- field_profile:

  For the fields: the grid of \\(\rho_0, \sigma_0)\\ with the gain in
  log-likelihood over the restricted fit, and the restricted
  log-likelihood.

- convergence, message:

  Convergence code and message of
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- max_grad, laplace_ok:

  Largest gradient component at the optimum, and whether the optimum
  passed the checks described in Estimation.

- tmb:

  The 'TMB' object, for advanced use.

Use the methods
[summary()](https://www.niccolosalvini.me/packages/zisar/reference/zisar-methods.md),
[coef()](https://www.niccolosalvini.me/packages/zisar/reference/zisar-methods.md),
[vcov()](https://www.niccolosalvini.me/packages/zisar/reference/zisar-methods.md),
[logLik()](https://www.niccolosalvini.me/packages/zisar/reference/zisar-methods.md),
[predict()](https://www.niccolosalvini.me/packages/zisar/reference/predict.zisar.md),
[residuals()](https://www.niccolosalvini.me/packages/zisar/reference/residuals.zisar.md)
and
[`zisar_effects()`](https://www.niccolosalvini.me/packages/zisar/reference/zisar_effects.md)
rather than the components directly.

## Model

For units \\i = 1, \dots, n\\ linked by a spatial weights matrix \\W\\,
\$\$\eta = \rho W \eta + X\beta + u, \qquad u \sim N(0, \sigma^2 I),\$\$
so that \\\eta = (I - \rho W)^{-1}(X\beta + u)\\. Each count is zero
with probability \\\pi_i\\ (a structural zero) and otherwise negative
binomial, \$\$y_i \sim \pi_i \\ 1\\y_i = 0\\ + (1 - \pi_i)\\
\mathrm{NB}(\mu_i, \kappa), \qquad \mu_i = \exp(\eta_i + o_i), \qquad
\mathrm{logit}(\pi_i) = \zeta_i,\$\$ with \\\mathrm{Var}(y_i \mid
\text{not a structural zero}) = \mu_i + \mu_i^2/\kappa\\ and \\o_i\\ an
optional offset (log exposure). The offset scales unit \\i\\ only: it
does not pass through the spatial multiplier.

The inflation predictor \\\zeta\\ is set by `zero_spatial`:

- `"none"`:

  The restricted model, \\\zeta = Z\gamma\\.

- `"reduced"`:

  A reduced-rank SAR field, \\\zeta = Z\gamma + B\alpha\\ with \\\alpha
  = \rho_0 W_K \alpha + e_0\\, \\e_0 \sim N(0, \sigma_0^2 I_K)\\:
  \\\alpha\\ holds the values of the field at the \\K\\ knots of a
  lattice, \\W_K\\ is the contiguity of the knots and \\B\\ a basis of
  local functions built by
  [`zisar_basis()`](https://www.niccolosalvini.me/packages/zisar/reference/zisar_basis.md),
  so that structural zeros can cluster in space beyond what the
  covariates explain. \\\rho_0 \ge 0\\. This is the recommended general
  model.

- `"field"`:

  A SAR field with one element per unit, \\\zeta = Z\gamma + g\\ with
  \\g = \rho_0 W g + e_0\\, \\e_0 \sim N(0, \sigma_0^2 I)\\. Its Laplace
  approximation fails (see Estimation).

- `"sar"`:

  A SAR filter on the predictor, \\\zeta = (I - \rho_0 W)^{-1}
  Z\gamma\\: the covariates of a unit's neighbours shift its inflation
  probability. No random effects are added.

The random fields are independent of \\u\\. The restricted model is the
limit \\\sigma_0 \to 0\\ of the fields, and \\\rho_0 = 0\\ of the
filter.

## Estimation

The latent field \\\eta\\ (and \\g\\ or \\\alpha\\) is integrated out
with the Laplace approximation in 'TMB', and the marginal likelihood is
maximised with [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html)
over \\(\beta, \gamma, \rho, \kappa, \sigma)\\ and, in the general
model, \\\rho_0\\ and \\\sigma_0\\, with \\\rho = \tanh(\rho\_{raw})\\
restricted to \\(-1, 1)\\. Standard errors come from the inverse Hessian
of the Laplace objective; those of \\\rho\\, \\\kappa\\, \\\sigma\\,
\\\rho_0\\ and \\\sigma_0\\ from the delta method.

With a random field in the inflation predictor (`"reduced"` or
`"field"`), the restricted model is fitted first, the Laplace likelihood
of the general model is evaluated on a grid of \\(\rho_0, \sigma_0)\\
with the other parameters at the restricted estimates, and the optimiser
starts from the best grid point. The grid and the gain over the
restricted fit are returned in `field_profile`; a gain that never rises
above zero means the data hold no field.

The unit-level field (`"field"`) defeats the Laplace approximation: each
\\g_i\\ is informed by a single outcome that is zero or not, and the
error of the approximation accumulates over the \\n\\ elements of the
field. On data without a field the approximate likelihood rises with
\\\sigma_0\\, and the likelihood-ratio statistic against the restricted
model can reach thousands. Use `"reduced"`, whose knot values are each
informed by many units, and check any general fit on data simulated
without a field ([`simulate()`](https://rdrr.io/r/stats/simulate.html)
of the restricted fit).

Each fit is checked: when the optimiser does not converge, or the
optimum has a gradient above 0.01 or a positive log-likelihood, the fit
is restarted from three other values of \\\kappa\\ and \\\sigma\\ and
from up to three random perturbations of the best point (this uses the
random number generator; `restarts = FALSE` skips it). If no clean
optimum is found, `zisar()` warns and `laplace_ok` is `FALSE`. A
positive log-likelihood is always a spurious optimum of the
approximation. A large gradient is not specific: it marks spurious
optima, but it also appears at genuine optima when structural zeros
occur among units of high intensity, since the conditional distribution
of \\\eta_i\\ of a zero unit then has two modes and the objective is not
smooth. Fit the model to data simulated from the fit to tell the two
apart, and use a parametric bootstrap
([`simulate()`](https://rdrr.io/r/stats/simulate.html) and refits
started from the estimates) for standard errors.

The log-determinant \\\log\|I - \rho W\|\\ is computed either from the
eigenvalues of \\W\\, found once in \\O(n^3)\\ time and \\O(n^2)\\
memory (`logdet = "eigen"`), or from a sparse Cholesky factorisation of
\\(I - \rho W)'(I - \rho W)\\ at each evaluation (`logdet = "sparse"`),
which scales to tens of thousands of units. `"auto"` uses the sparse
method when \\n \> 5000\\. The knot field uses the eigenvalues of
\\W_K\\.

## Nested models

`fixed` holds parameters at given values, which gives the special cases
of the model: `fixed = list(rho = 0)` is the aspatial zero-inflated
negative binomial with random effects (ZINB-RE), and with an
intercept-only inflation part (`y ~ x | 1`), `fixed = list(gamma = -20)`
sets \\\pi_i \approx 0\\, the SAR negative binomial with random effects
(SAR-NB-RE). Likelihood-ratio tests between nested fits use
[logLik()](https://www.niccolosalvini.me/packages/zisar/reference/zisar-methods.md).
The restricted model sits on the boundary of the models with a field
(\\\sigma_0 = 0\\), so the likelihood-ratio statistic between them does
not have the usual chi-squared distribution: its p-value from
\\\chi^2_2\\ is conservative.

## References

Salvini, N., Arbia, G. and Nardelli, V. (2026). ZISAR: A zero-inflated
spatial autoregressive model for count data with random effects. Working
paper, Università Cattolica del Sacro Cuore.

Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H. and Bell, B. M.
(2016). TMB: Automatic differentiation and Laplace approximation.
*Journal of Statistical Software*, 70(5), 1–21.
[doi:10.18637/jss.v070.i05](https://doi.org/10.18637/jss.v070.i05)

## See also

[`zisar_effects()`](https://www.niccolosalvini.me/packages/zisar/reference/zisar_effects.md)
for direct, indirect and total effects,
[`zisar_basis()`](https://www.niccolosalvini.me/packages/zisar/reference/zisar_basis.md)
for the reduced-rank field,
[`rzisar()`](https://www.niccolosalvini.me/packages/zisar/reference/rzisar.md)
to simulate from the model.

## Examples

``` r
set.seed(1)
sim <- rzisar(12, beta = c(1, 0.5), rho = 0.4, gamma = c(-1, 0.8))
fit <- zisar(y ~ x1 | z1, data = sim$data, W = sim$W)
summary(fit)
#> 
#> Call:
#> zisar(formula = y ~ x1 | z1, data = sim$data, W = sim$W)
#> 
#> 144 spatial units, 38.9% zeros
#> 
#> Count model coefficients (SAR negative binomial, log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   1.6222     0.8087   2.006   0.0449 *  
#> x1            0.6690     0.1141   5.865 4.49e-09 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Zero-inflation model coefficients (binomial, logit link):
#>             Estimate Std. Error z value Pr(>|z|)   
#> (Intercept)  -0.8054     0.2642  -3.048  0.00230 **
#> z1            0.6516     0.2243   2.905  0.00367 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Spatial and dispersion parameters:
#>       Estimate Std. Error z value Pr(>|z|)
#> rho     0.1056     0.4121   0.256    0.798
#> kappa   1.4609     0.8026                 
#> sigma   0.2034     0.6278                 
#> 
#> Log-likelihood (Laplace): -344.303 on 7 Df
#> AIC: 702.606   BIC: 723.394
#> Log-determinant: eigen.  Optimiser: converged after 24 iterations.
#> 

# Aspatial ZINB with random effects, and a likelihood-ratio test of rho = 0
fit0 <- zisar(y ~ x1 | z1, data = sim$data, W = sim$W,
              fixed = list(rho = 0))
lr <- 2 * (logLik(fit) - logLik(fit0))
pchisq(as.numeric(lr), df = 1, lower.tail = FALSE)
#> [1] 0.8091071

# \donttest{
# General model: zeros clustered by a spatial field of their own, fitted
# with a reduced-rank field on knots 4 cells apart
set.seed(2)
sim2 <- rzisar(20, beta = c(1, 0.5), rho = 0.4, gamma = c(-1, 0.8),
               zero_spatial = "field", rho0 = 0.9, sigma0 = 0.5)
b <- zisar_basis(sim2$coords, h = 4)
fit2 <- zisar(y ~ x1 | z1, data = sim2$data, W = sim2$W,
              zero_spatial = "reduced", basis = b)
summary(fit2)
#> 
#> Call:
#> zisar(formula = y ~ x1 | z1, data = sim2$data, W = sim2$W, zero_spatial = "reduced", 
#>     basis = b)
#> 
#> 400 spatial units, 40.8% zeros
#> 
#> Count model coefficients (SAR negative binomial, log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  1.21972    0.48548   2.512    0.012 *  
#> x1           0.59762    0.06994   8.545   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Zero-inflation model coefficients (binomial, logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -0.9816     0.2559  -3.836 0.000125 ***
#> z1            0.8217     0.1836   4.474 7.66e-06 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Spatial and dispersion parameters (rho0, sigma0: reduced-rank inflation field, 36 knots):
#>        Estimate Std. Error z value Pr(>|z|)
#> rho      0.3119     0.2558   1.219    0.223
#> kappa    1.1170     0.1783                 
#> sigma    0.0010     0.2721                 
#> rho0     0.0000     1.0220   0.000    1.000
#> sigma0   0.7970     0.6941                 
#> 
#> Log-likelihood (Laplace): -968.571 on 9 Df
#> AIC: 1955.14   BIC: 1991.07
#> Log-determinant: eigen.  Optimiser: converged after 23 iterations.
#> 
# }
```
