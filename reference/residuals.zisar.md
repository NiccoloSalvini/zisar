# Residuals of a fitted ZISAR model

Response residuals \\y_i - E\[y_i\]\\ or Pearson residuals, both
conditional on the estimated latent field. Pearson residuals divide by
the standard deviation of the zero-inflated negative binomial mixture,
\\\mathrm{Var}(y_i) = (1 - \pi_i)\mu_i(1 + \mu_i/\kappa + \pi_i\mu_i)\\.

## Usage

``` r
# S3 method for class 'zisar'
residuals(object, type = c("pearson", "response"), ...)
```

## Arguments

- object:

  A fitted `"zisar"` object.

- type:

  `"pearson"` or `"response"`.

- ...:

  Not used.

## Value

A numeric vector of residuals.

## Details

Because the residuals condition on the mode of \\\eta\\, which already
absorbs most of the spatial structure, tests such as Moran's I on them
have little power against remaining spatial dependence; the
likelihood-ratio test of \\\rho = 0\\ is preferable (see
[`zisar()`](https://www.niccolosalvini.me/packages/zisar/reference/zisar.md)).
