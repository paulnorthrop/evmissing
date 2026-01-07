#' evmissing: Extreme Value Analyses with Missing Data
#'
#' @description Performs likelihood-based extreme value inferences with
#' adjustment for the presence of missing values. A Generalised Extreme Value
#' (GEV) distribution is fitted to block maxima using maximum likelihood
#' estimation, with the GEV location and scale parameters reflecting the
#' numbers of non-missing raw values in each block. A Bayesian version is also
#' provided. For the purposes of comparison, there are options to make no
#' adjustment for missing values or to discard any block maximum for which
#' greater than a percentage of the underlying raw values are missing.
#'
#' The `evmissing` package was created to accompany Simpson, E. S. and
#' Northrop, P. J. (2025) Accounting for missing data when modelling block
#' maxima. \doi{10.48550/arXiv.2512.15429}
#'
#' @details The main functions are
#'
#' * [`gev_mle`]: maximum likelihood inference for block maxima based on a GEV
#'   distribution, with [`S3 methods`][evmissing_methods] including `confint`.
#' * [`gev_bayes`]: Bayesian inference for block maxima based on a GEV
#'   distribution.
#'
#' For objects returned by `gev_mle`, inferences about return levels are
#' performed by [`gev_return`], with with [`S3 methods`][return_level_methods]
#' including `confint`.
#'
#' The function [`gev_influence`] quantifies the influence that individual
#' extreme (small or large) block maxima have on the maximum likelihood
#' estimators of GEV parameters.
#'
#' The following example datasets are provided.
#'
#' * [`BloomsburyOzoneMaxima`]: Annual maxima ozone levels at Bloomsbury,
#'   London, UK, 1992-2024.
#' * [`PlymouthOzoneMaxima`]: Annual maxima ozone levels at Plymouth, Devon,
#'   UK, 1998-2024.
#' * [`BrestSurgeMaxima`]: Annual maxima surge heights at Brest, France,
#'   1846-2007.
#'
#' @docType package
#' @aliases evmissing-package
#' @import revdbayes
#' @importFrom stats nobs vcov coef logLik confint
#' @importFrom graphics plot
"_PACKAGE"
