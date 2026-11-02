#' @keywords internal
#' @details
#' The main function is [zisar()], which fits the model from a two-part
#' formula `y ~ count terms | zero-inflation terms` and a spatial weights
#' matrix. [zisar_fit()] is the matrix interface underneath it. Fitted models
#' have the usual methods ([summary()][summary.zisar], [coef()][coef.zisar],
#' [vcov()][vcov.zisar], [logLik()][logLik.zisar], [predict()][predict.zisar],
#' [residuals()][residuals.zisar]); [zisar_effects()] decomposes covariate
#' effects into direct, indirect and total, and [rzisar()] simulates data
#' from the model.
#'
#' @references
#' Salvini, N., Arbia, G. and Nardelli, V. (2026). ZISAR: A zero-inflated
#' spatial autoregressive model for count data with random effects.
#' Working paper, Università Cattolica del Sacro Cuore.
"_PACKAGE"

## usethis namespace: start
#' @useDynLib zisar, .registration = TRUE
#' @importFrom methods as
#' @importFrom stats coef fitted logLik model.frame model.matrix model.offset
#'   model.response nlminb nobs optimHess plogis pnorm predict qlogis quantile
#'   rbinom residuals rnbinom rnorm sd simulate terms vcov
## usethis namespace: end
NULL
