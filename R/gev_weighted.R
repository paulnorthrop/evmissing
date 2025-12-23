#' Weighted GEV ML Inference with Adjustment for Missing Data
#'
#' Fits a GEV distribution to block maxima using maximum likelihood estimation,
#' with the option to make an adjustment for the numbers of non-missing raw
#' values in each block using one of the two weighting schemes proposed in
#' McVittie and Murphy (2025).
#'
#' @param data Either
#'
#'   * a numeric vector containing a time series of raw data,
#'   * an object returned from [`block_maxima`], a list with components
#'     `maxima`, `notNA` and `n`,
#'   * a data frame or named list containing the same information (variables
#'     `maxima`, `notNA` and `n`) as an object returned from [`block_maxima`],
#'     such as the data frame [`BrestSurgeMaxima`].
#'
#' @param block_length A numeric scalar. Used calculate the maxima of disjoint
#'   blocks of `block_length` contiguous values in the vector `data`.
#'   If `length(data)` is not an integer multiple of `block_length` then
#'   the values at the end of `data` that do not constitute a complete block
#'   of length `block_length` are discarded, without warning.
#' @param block A numeric vector with the same length as `data`. The value of
#'   `block[i]` indicates the block into which `data[i]` falls. For example,
#'   `block` could provide the year in which observation `i` was observed.
#' @param scheme A numeric scalar. Pass `scheme = 1` for the first weighting
#'  scheme, or `scheme = 2` for the second weighting scheme, in McVittie and
#'  Murphy (2025). Any value other than `1` or `2` result in an unweighted fit,
#'  that is, all weight are set to 1. See **Details**.
#'
#' @param init Either a character scalar, one of `"quartiles"` or `"moments"`,
#'   or a numeric vector of length 3 giving initial estimates of the GEV
#'   location, scale and shape parameters: \eqn{\mu}, \eqn{\sigma} and
#'   \eqn{\xi}. If `init = "quartiles"` then initial estimates of \eqn{\mu} and
#'   \eqn{\sigma} are based on sample quartiles of block maxima, ignoring the
#'   underlying numbers of non-missing raw data, and a value of 0 for
#'   \eqn{\xi}. If `init = "moments"` then instead we use the sample mean and
#'   variance of these maxima and an initial value of 0.1 for \eqn{\xi}.
#' @param ... Further arguments to be passed to [`stats::optim`].
#' @details Suppose that a full (no missing values) block will contain \eqn{n}
#'   raw values. Let \eqn{n_i} be the number of non-missing values, and
#'   \eqn{m_i} the observed block maximum, in block \eqn{i}.
#'   The contribution of block maximum \eqn{m_i} to the GEV log-likelihood is
#'   weighted (multiplied) by the weight \eqn{w_i}. The set of weights depends
#'   on the weighting scheme chosen by `scheme` as follows.
#'
#'  * If `scheme = 1` then \eqn{w_i = n_i / n}, for \eqn{i = 1, ..., n}.
#'  * If `scheme = 2` then \eqn{w_i = \hat{F}(m_i) ^ {n - n_i}}, for
#'    \eqn{i = 1, ..., n}, where \eqn{\hat{F}} is the empirical distribution
#'    function of \eqn{m_1, ..., m_n}.
#'
#'  For any other value of `scheme` all weights are set to 1, that is, an
#'  unweighted fit is performed.
#'
#' For each weighting scheme, the larger the number \eqn{n - n_i} of missing
#' values the smaller the weight.
#' See McVittie and Murphy (2025) for further details.
#'
#' @return A list, returned from [`stats::optim`] (the MLEs are in the
#' component `par`), with the additional components:
#'
#' * `loglik`: value of the maximised log-likelihood.
#' * `vcov`: estimated variance-covariance matrix of the parameters.
#' * `se`: estimated standard errors of the parameters.
#' * `maxima`: the vector of block maxima used to fit the model.
#' * `notNA`: the number of non-missing raw values on which the maxima in
#'   `maxima` are based.
#' * `n`: the maximal block length, that is, the largest number of values that
#'     could have been observed in each of these blocks.
#' * `weights`: the weights used in the fitting.
#'
#' The call to `gev_mle` is provided in the attribute `"call"`.
#'
#' The class of the returned object is
#' `c("evmissing", "weighted_mle", "list")`.
#'
#' Objects inheriting from class `"evmissing"` have `coef`, `logLik`, `nobs`,
#' `summary`, `vcov` and `confint` methods.  See [`evmissing_methods`].
#' @references McVittie, J. H. and Murphy, O. A. (2025) Weighted Parameter
#'   Estimators of the Generalized Extreme Value Distribution in the Presence
#'   of Missing Observations. \doi{10.48550/arXiv.2506.15964}
#' @examples
#' ## Simulate raw data from an exponential distribution
#'
#' set.seed(13032025)
#' blocks <- 50
#' block_length <- 365
#' sdata <- sim_data(blocks = blocks, block_length = block_length)
#'
#' # sdata$data_full have no missing values
#' # sdata$data_miss have had missing values created artificially
#'
#' ## Fits to full data: fit0, fit 1 and fit2 should give the same results
#'
#' # Fit a GEV distribution to block maxima from the full data
#' fit0 <- gev_mle(sdata$data_full, block_length = sdata$block_length)
#' summary(fit0)
#'
#' # Fit to the full data using weighting scheme 1
#' fit1 <- gev_weighted(sdata$data_full, scheme = 1,
#'                      block_length = sdata$block_length)
#' summary(fit1)
#'
#' # Fit to the full data using weighting scheme 2
#' fit2 <- gev_weighted(sdata$data_full, scheme = 2,
#'                      block_length = sdata$block_length)
#' summary(fit2)
#'
#' ## Fits to the data with missings
#'
#' #  Make adjustment for the numbers of non-missing values per block
#' fit0 <- gev_mle(sdata$data_miss, block_length = sdata$block_length)
#' summary(fit0)
#'
#' # Make adjustment using weighting scheme 1
#' fit1 <- gev_weighted(sdata$data_miss, scheme = 1,
#'                      block_length = sdata$block_length)
#' summary(fit1)
#'
#' # Make adjustment using weighting scheme 2
#' fit2 <- gev_weighted(sdata$data_miss, scheme = 2,
#'                      block_length = sdata$block_length)
#' summary(fit2)
#' @export
gev_weighted <- function(data, scheme = 1, block_length, block,
                         init = "quartiles", ...) {

  # If data was created by block_maxima() or is a data frame that contains the
  # correct information then use it. Otherwise, use block_maxima() to calculate
  # the block maxima, the numbers of non-missing values in the blocks and the
  # largest possible number of non-missing values in each block
  if (inherits(data, "block_maxima") && inherits(data, "evmissing")) {
    maxima_notNA <- data
  } else if (is.data.frame(data)) {
    if (all(is.element(c("maxima", "notNA", "n"), colnames(data)))) {
      maxima_notNA <- as.list(data)
    } else {
      stop("Data frame ''data'' does not contain the required variables.")
    }
  } else if (is.list(data)) {
    if (all(is.element(c("maxima", "notNA", "n"), names(data)))) {
      maxima_notNA <- as.list(data)
    } else {
      stop("List ''data'' does not contain the required named components.")
    }
  } else {
    maxima_notNA <- block_maxima(data, block_length, block)
  }
  # If there are maxima = NA, notNA = 0 entries in the data then remove them
  no_data <- which(maxima_notNA$notNA == 0)
  if (length(no_data) > 0) {
    maxima_notNA <- lapply(maxima_notNA, function(x) x[-no_data])
  }
  # maxima_notNA is a list with 3 components
  #   maxima: the block maxima, the response
  #    notNA: the numbers of non-missing values in each block, the covariate
  #        n: the largest possible number of non-missing values in each block

  # If init is a has not been supplied then calculate initial estimates of mu
  # and sigma for assuming that xi = 0
  if (is.character(init)) {
    init_method <- match.arg(init, c("quartiles", "moments"))
    init <- gev_init(maxima_notNA, init_method = init_method)
  } else {
    names(init) <- c("mu", "sigma", "xi")
  }
  # Minimise the negated log-likelihood with respect to (mu, sigma, xi)
  optim_args <- list(...)
  # The default method is "Nelder-Mead". If optim_args$method is non-NULL then
  # the user has specified the method via ...
  # The optim method "L-BFGS-B" doesn't like returned values of Inf or NA
  # If this method is used then set the out-of-bounds big value big_val to
  # 10 ^ 10. Otherwise, use big_val = Inf.
  big_val <- Inf
  if (!is.null(optim_args$method)) {
    if (optim_args$method == "L-BFGS-B") {
      big_val <- 10 ^ 10
    }
  }
  # The weights depend only on the observed maxima and/or the number of
  # missing values in each block so calculate them once now rather than
  # repeatedly in the weighted negated log-likelihood functions
  # Set the weights vector
  # Scheme = 1: proportions of missing values per block
  # Scheme = 2: ecdf evaluated at the block maximum to the power of the number
  #             of missing values
  # Extract the block maxima, the numbers of non-missing values in the blocks
  # and the the largest possible number of non-missing values in each block
  maxima <- maxima_notNA$maxima
  n_i <- maxima_notNA$notNA
  n <- maxima_notNA$n
  if (scheme == 1) {
    weights <- n_i / n
  } else if (scheme == 2) {
    Fhat <- stats::ecdf(data)
    weights <- Fhat(maxima) ^ (n - n_i)
  } else {
    weights <- 1
  }
  fit <- try(stats::optim(par = init, fn = weighted_negated_gev_loglik,
                          maxima = maxima, weights = weights, hessian = TRUE,
                          ..., big_val = big_val), silent = TRUE)
  # If there is an error then return NA values for the estimates etc.
  if (inherits(fit, "try-error")) {
    optim_error <- attr(fit, "condition")
    fit <- list()
    fit$optim_error <- optim_error
    fit$par <- rep(NA, 3)
    par_names <- c("mu", "sigma", "xi")
    names(fit$par) <- par_names
    fit$vcov <- matrix(NA, 3, 3)
    dimnames(fit$vcov) <- list(par_names, par_names)
    fit$se <- fit$par
    fit$loglik <- NA
    warning("MLEs not found: NA values returned.")
  } else {
    var_cov <- qr(fit$hessian)
    if (var_cov$rank != ncol(var_cov$qr)) {
      fit$vcov <- matrix(NA, 3, 3)
      dimnames(fit$vcov) <- list(names(fit$par), names(fit$par))
      warning("The observed information matrix is singular.")
    } else {
      fit$vcov <- solve(fit$hessian)
    }
    fit$se <- sqrt(diag(fit$vcov))
    if (any(is.na(fit$se)) || any(fit$se <= 0)) {
      fit$se <- rep(NA, 3)
      names(fit$se) <- names(fit$par)
      warning("The observed information matrix is singular.")
    }
    fit$loglik <- -fit$value
  }
  fit$maxima <- maxima_notNA$maxima
  fit$notNA <- maxima_notNA$notNA
  fit$n <- maxima_notNA$n
  fit$scheme <- scheme
  fit$weights <- weights
  attr(fit, "call") <- match.call(expand.dots = TRUE)
  class(fit) <- c("evmissing", "weighted_mle", class(fit))
  return(fit)
}
