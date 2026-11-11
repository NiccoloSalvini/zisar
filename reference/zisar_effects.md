# Direct, indirect and total effects of a ZISAR model

Decomposes the effect of each covariate on the expected count \\E\[y\] =
(1 - \pi)\mu\\ into a direct part (a unit's own covariate), an indirect
part (spillovers from the covariates of other units through the spatial
multiplier) and their sum, averaged over units as in LeSage and Pace
(2009).

## Usage

``` r
zisar_effects(object, B = 500, level = 0.95)

# S3 method for class 'zisar_effects'
print(x, digits = 4, ...)
```

## Arguments

- object:

  A fitted `"zisar"` object.

- B:

  Number of simulation draws for standard errors and intervals; `0`
  gives point estimates only.

- level:

  Confidence level of the simulation intervals.

- x:

  A `"zisar_effects"` object.

- digits:

  Number of decimal places to print.

- ...:

  Not used.

## Value

A data frame of class `"zisar_effects"`, one row per covariate and
effect (`direct`, `indirect`, `total`), with the point estimate and,
when `B > 0`, the simulation standard error and interval bounds.

## Details

The matrix of partial derivatives of \\E\[y_i\]\\ with respect to
\\x\_{jr}\\ is \\\mathrm{diag}(d)\\ S(\rho)\\ \beta_r\\ with \\d_i =
(1 - \pi_i)\mu_i\\ and \\S(\rho) = (I - \rho W)^{-1}\\. The average
direct effect is the mean of its diagonal, the average total effect the
mean of its row sums, and the indirect effect their difference. Only the
diagonal and the row sums of \\S(\rho)\\ are computed, so no dense \\n
\times n\\ matrix is formed for large \\n\\.

A covariate that also enters the zero-inflation part (same name in both
parts of the formula) adds the inflation channel \\-\mu_i \pi_i (1 -
\pi_i) \gamma_r\\, which acts on the unit itself only and so enters the
direct and total effects. A covariate of the inflation part alone has
this direct effect only. With the SAR filter on the inflation predictor
(`zero_spatial = "sar"`) the inflation channel passes through \\S_0 =
(I - \rho_0 W)^{-1}\\ as well, and has an indirect part.

Standard errors and intervals are simulation-based: parameters are drawn
from their estimated asymptotic normal distribution and the effects
recomputed for each draw, keeping the predicted random effects \\\hat u
= (I - \rho W)\hat\eta - X\beta\\ (and an inflation field at its mode)
at their estimates. For \\n \> 3000\\, the diagonal of \\S(\rho)\\ along
the draws is interpolated between seven exact evaluations.

## References

LeSage, J. and Pace, R. K. (2009). *Introduction to Spatial
Econometrics*. Chapman and Hall/CRC.

## Examples

``` r
set.seed(1)
sim <- rzisar(12, beta = c(1, 0.5), rho = 0.4, gamma = c(-1, 0.8))
fit <- zisar(y ~ x1 | z1, data = sim$data, W = sim$W)
zisar_effects(fit, B = 200)
#> ZISAR average effects on E[y] (simulation SEs, B = 200)
#> 
#>              Direct        Indirect            Total
#> x1  3.4788 (1.4814) 0.4043 (1.6099)  3.8832 (1.7326)
#> z1 -0.9652 (0.4093) 0.0000 (0.0000) -0.9652 (0.4093)
```
