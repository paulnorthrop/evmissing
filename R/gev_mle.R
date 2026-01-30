#' GEV ML Inference with Adjustment for Missing Data
#'
#' Fits a GEV distribution to block maxima using maximum likelihood estimation,
#' with the option to make an adjustment for the numbers of non-missing raw
#' values in each block. The GEV location and scale parameters are adjusted to
#' reflect the proportion of raw values that are missing.
#'
#' @param data Either
#'
#'   * a numeric vector containing a time series of raw data,
#'   * an object returned from [`block_maxima`], a list with components
#'     `maxima`, `notNA` and `n`,
#'   * a data frame or named list containing the same information, that is, the
#'     variables `maxima`, `notNA` and `n`, as an object returned from
#'     [`block_maxima`], such as the data frame [`BrestSurgeMaxima`].
#'
#' @param block_length A numeric scalar. Used calculate the maxima of disjoint
#'   blocks of `block_length` contiguous values in the vector `data`.
#'   If `length(data)` is not an integer multiple of `block_length` then
#'   the values at the end of `data` that do not constitute a complete block
#'   of length `block_length` are discarded, without warning.
#' @param block A numeric vector with the same length as `data`. The value of
#'   `block[i]` indicates the block into which `data[i]` falls. For example,
#'   `block` could provide the year in which observation `i` was observed.
#' @param adjust A logical scalar or a numeric scalar in `[0, 100]`.
#'
#'  * If `adjust = TRUE` then the adjustment, described in **Details**, for the
#'    numbers of non-missing values underlying each block maximum is performed.
#'  * If `adjust = FALSE` then no adjustment is made, that is, the block maxima
#'    are treated as if the underlying raw data have no missing values.
#'
#' @param discard A numeric scalar. Any block maximum for which greater than
#'   `discard` percent of the underlying raw values were missing is discarded.
#'   Whether or not an adjustment for missingness is made for the block maxima
#'   that remain is determined by `adjust`.
#' @param init Either a character scalar, one of `"quartiles"` or `"moments"`,
#'   or a numeric vector of length 3 giving initial estimates of the GEV
#'   location, scale and shape parameters: \eqn{\mu}, \eqn{\sigma} and
#'   \eqn{\xi}. If `init = "quartiles"` then initial estimates of \eqn{\mu} and
#'   \eqn{\sigma} are based on sample quartiles of block maxima, ignoring the
#'   underlying numbers of non-missing raw data, and a value of 0 for
#'   \eqn{\xi}. If `init = "moments"` then instead we use the sample mean and
#'   variance of these maxima and an initial value of 0.1 for \eqn{\xi}.
#' @param ... Further arguments to be passed to [`stats::optim`].
#' @details If `data` is numeric vector then exactly one of the arguments
#'   `block_length` or `block` must be supplied. The parameters are fitted
#'   using maximum likelihood estimation.
#'
#' The adjustment for the numbers of non-missing values underlying the block
#' maxima is based on the strong assumption that missing values occur
#' completely at random. We suppose that a block maximum \eqn{M_n} based on
#' a full (no missing raw values) block of length \eqn{n} has a
#' \eqn{\text{GEV}(\mu, \sigma, \xi)} distribution, with distribution function
#' \eqn{G(x)}. Let \eqn{n_i} be the number of non-missing values in block \eqn{i}
#' and let \eqn{M_{n_i}} denote the block maximum of such a block. We suppose
#' that \eqn{M_{n_i}} has a \eqn{\text{GEV}(\mu(n_i), \sigma(n_i), \xi)}
#' distribution, where \deqn{\mu(n_i) = \mu + \sigma [(n_i/n)^\xi -1] / \xi,}
#' \deqn{\sigma(n_i) = \sigma (n_i/n)^\xi.}
#'
#' These expressions are based on inferring the parameters of an approximate
#' GEV distribution for \eqn{M_{n_i}} from its approximate distribution function
#' \eqn{[G(x)]^{n_i/n}}.
#'
#' A likelihood is constructed as the product of contributions from the maxima
#' from distinct blocks, under the assumption that these maxima are
#' independent. Let \eqn{\theta = (\mu, \sigma, \xi)} and let
#' \eqn{\ell_F(\underline{z}; \theta)} denote the usual, unadjusted, GEV
#' log-likelihood for the full-data case where there are no missing values.
#' It can be shown that our adjusted log-likelihood
#' \eqn{\ell(\theta, \underline{z})} is given by
#'
#' \deqn{\ell(\theta, \underline{z}) = \ell_F(\underline{z}; \theta) -
#'       \sum_{i=1}^n p_i \log G(z_i; \theta)}
#'
#' where \eqn{p_i = 1 - n_i / n} is the proportion of missing values in block
#' \eqn{i}.
#'
#' The negated log-likelihood is minimised using a call to
#' [`stats::optim`] with `hessian = TRUE`. If [`stats::optim`] throws an error
#' then a warning is produced and the returned object has `NA` values for
#' the components `par`, `loglik`, `vcov` and `se` and an extra component
#' `optim_error` containing the error message. If the estimated observed
#' information matrix is singular then a warning is produced and the returned
#' object has `NA` values for the components `vcov` and `se`.
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
#' * `adjust,discard` : the values of these input arguments.
#'
#' The call to `gev_mle` is provided in the attribute `"call"`.
#'
#' The class of the returned object is `c("evmissing", "mle", "list")`.
#'
#' Objects inheriting from class `"evmissing"` have `coef`, `logLik`, `nobs`,
#' `summary`, `vcov` and `confint` methods.  See [`evmissing_methods`].
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
#' # Fit a GEV distribution to block maxima from the full data
#' fit1 <- gev_mle(sdata$data_full, block_length = sdata$block_length)
#' summary(fit1)
#'
#' # An identical fit supplying the block indicator instead of block_length
#' fit1b <- gev_mle(sdata$data_full, block = sdata$block)
#' summary(fit1b)
#'
#' # Make adjustment for the numbers of non-missing values per block
#' fit2 <- gev_mle(sdata$data_miss, block_length = sdata$block_length)
#' summary(fit2)
#'
#' # Do not make the adjustment
#' fit3 <- gev_mle(sdata$data_miss, block_length = sdata$block_length,
#'                 adjust = FALSE)
#' summary(fit3)
#'
#' # Remove all block maxima with greater than 25% missing values and
#' # do not make the adjustment
#' fit4 <- gev_mle(sdata$data_miss, block_length = sdata$block_length,
#'                 adjust = FALSE, discard = 25)
#' summary(fit4)
#'
#' ## Plymouth ozone data
#'
#' # Calculate the values in Table 4 of Simpson and Northrop (2026)
#' # discard = 50 is chosen to discard data from 2001 and 2006
#' fit1 <- gev_mle(PlymouthOzoneMaxima)
#' fit2 <- gev_mle(PlymouthOzoneMaxima, adjust = FALSE)
#' fit3 <- gev_mle(PlymouthOzoneMaxima, discard = 50)
#' fit4 <- gev_mle(PlymouthOzoneMaxima, adjust = FALSE, discard = 50)
#' se <- function(x) return(sqrt(diag(vcov(x))))
#' MLEs <- cbind(coef(fit1), coef(fit2), coef(fit3), coef(fit4))
#' SEs <- cbind(se(fit1), se(fit2), se(fit3), se(fit4))
#' round(MLEs, 2)
#' round(SEs, 2)
#' @export
gev_mle <- function(data, block_length, block, adjust = TRUE, discard = 0,
                    init = "quartiles", ...) {
  # Check discard
  if (!is.numeric(discard) || any(discard < 0)) {
    stop("''discard'' must be positive number")
  }
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
  # If discard >  0 then discard any block maxima based on underlying data with
  # greater than discard% missing.
  if (discard > 0) {
    maxima_notNA <- discard_maxima(maxima_notNA, discard = discard)
  }
  # maxima_notNA is a list with 3 components
  #   maxima: the block maxima, the response
  #    notNA: the numbers of non-missing values in each block, the covariate
  #        n: the largest possible number of non-missing values in each block

  # Fit a model in which block maximum i of n has a GEV distribution with
  #     location: mu + sigma * [(n_i / n) ^ xi - 1] / xi,
  #   log(scale): log(sigma) + xi * log (n_i / n),
  #        shape: xi,
  #   where n is the block length, n_i is the number of non-missing values in
  #   block i and the block maxima are conditionally independent given the n_is
  #   (n_1, ..., n_b).

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
  fit <- try(stats::optim(par = init, fn = negated_gev_loglik, hessian = TRUE,
                          ..., maxima_notNA = maxima_notNA, adjust = adjust,
                          big_val = big_val), silent = TRUE)
  # If there is an error then it is probably because the sample size is small
  # (perhaps several block maxima have been discarded) and the log-likelihood
  # does not have a maximum that is away from the boundary of the parameter
  # space. If this happens then return NA values for the estimates etc.
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
  fit$adjust <- adjust
  fit$discard <- discard
  attr(fit, "call") <- match.call(expand.dots = TRUE)
  class(fit) <- c("evmissing", "mle", class(fit))
  return(fit)
}
