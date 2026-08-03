#' GEV ML Inference with Adjustment for Missing Data (Stationary Sequences)
#'
#' Fits a GEV distribution to block maxima using maximum likelihood estimation,
#' making an adjustment for the locations of missing raw values in each block.
#' The GEV location and scale parameters are adjusted to reflect the positions
#' of missing raw values and time series dependence in the data.
#'
#' @param data Either
#'
#'   * a numeric vector containing a time series of raw data,
#'   * an object returned from a call to [`block_maxima`] with `pseudo = TRUE`,
#'     such as [`LaPlagnePrecipMaximaSeasonal`].
#'
#' @param block_length A numeric scalar. Used calculate the maxima of
#'   disjoint blocks of `block_length` contiguous values in the vector `data`.
#'   If `length(data)` is not an integer multiple of `block_length` then
#'   the values at the end of `data` that do not constitute a complete block
#'   of length `block_length` are discarded, without warning.
#' @param block A numeric vector with the same length as `data`. The value of
#'   `block[i]` indicates the block into which `data[i]` falls. The block
#'   lengths implied by `block` may differ by at most 1. For example,
#'   `block[i]` could give the block in which observation `i` was observed,
#'   with block lengths of 366 for leap years and 365 for other years.
#' @param pseudo A logical scalar. If `pseudo = TRUE` then the pseudo-maxima
#'   returned from [`block_maxima`] are used to estimate the value of
#'   \eqn{r_i} for an incomplete, partially-observed block. See **Details**
#'   and the **Details** section of [`block_maxima`]. If `pseudo = FALSE` then
#'   the ratio \eqn{n_i/n} is used, as in [`gev_mle()`].
#' @param full,sliding,seasonal Arguments passed to [`block_maxima`] to
#'   determine how the pseudo-maxima are calculated.
#' @param init Either a character scalar, one of `"quartiles"` or `"moments"`,
#'   or a numeric vector of length 3 giving initial estimates of the GEV
#'   location, scale and shape parameters: \eqn{\mu}, \eqn{\sigma} and
#'   \eqn{\xi}. If `init = "quartiles"` then initial estimates of \eqn{\mu} and
#'   \eqn{\sigma} are based on sample quartiles of block maxima, ignoring the
#'   underlying numbers of non-missing raw data, and a value of 0 for
#'   \eqn{\xi}. If `init = "moments"` then instead we use the sample mean and
#'   variance of these maxima and an initial value of 0.1 for \eqn{\xi}.
#' @param ... Further arguments to be passed to [`stats::optim`].
#' @details If `data` is a numeric vector then exactly one of the arguments
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
#' that, \eqn{M_{n_i}} has a \eqn{\text{GEV}(\mu(r_i), \sigma(r_i), \xi)}
#' distribution, where \deqn{\mu(r_i) = \mu + \sigma [r_i^\xi -1] / \xi,}
#' \deqn{\sigma(r_i) = \sigma r_i^\xi,} for some \eqn{n_i/n \leq r_i \leq 1}.
#' These expressions are based on the \eqn{M_{n_i}} having approximately a
#' GEV distribution with distribution function \eqn{G(x)^{r_i}}.
#'
#' For a full block, \eqn{r_i = 1}. If `pseudo = TRUE`, then, for a
#' partially-observed block, the value of \eqn{r_i} is estimated using the
#' pseudo-maxima returned from [`block_maxima`] and the GEV distribution
#' function based on the current value of \eqn{(\mu, \sigma, \xi)} in the
#' optimisation routine. Suppose that we have a vector \eqn{M_i} of
#' pseudo-maxima resulting from a particular partially-observed block \eqn{i}.
#' It can be shown that the components of \eqn{V_i = -\log G(M_i)} each have an
#' exponential distribution with mean \eqn{1/r_i}. We estimate \eqn{r_i} using
#' the reciprocal of the mean of the values in \eqn{V_i}.
#'
#' The negated log-likelihood is minimised using a call to
#' [`stats::optim`] with `hessian = TRUE`. If [`stats::optim`] throws an error
#' then a warning is produced and the returned object has `NA` values for
#' the components `par`, `loglik`, `vcov` and `se` and an extra component
#' `optim_error` containing the error message. If the estimated observed
#' information matrix is singular then a warning is produced and the returned
#' object has `NA` values for the components `vcov` and `se`.
#'
#' @return A list, returned from [`stats::optim`] (the MLEs are in the
#' component `par`), with the additional components:
#'
#' * `loglik`: value of the maximised log-likelihood.
#' * `vcov`: estimated variance-covariance matrix of the parameters.
#' * `se`: estimated standard errors of the parameters.
#' * `maxima`: the vector of block maxima used to fit the model.
#' * `full_maxima,partial_maxima`: vectors of full and partial block maxima.
#' * `notNA`: the number of non-missing raw values on which the maxima in
#'   `maxima` are based.
#' * `n`: the maximal block length, that is, the largest number of values that
#'     could have been observed in each of these blocks.
#' * `rvec`: a vector of the values used for \eqn{r_1, ..., r_b}, where \eqn{b}
#'   is the number of blocks. The content depends on the argument `pseudo`.
#' * `rhats`: a vector of the subset of `rvec` for partially-observed blocks
#'    only. The attributes `"propn_notNA"` and `"unconstrained"` give,
#'    respectively, the values of \eqn{n_i/n} for these blocks and the
#'    estimates of \eqn{r_i} before they are constrained to lie in the interval
#'    \eqn{[n_i/n, 1]}.
#' * `sliding`: the input argument `sliding`.
#'
#' The call to `gev_ts` is provided in the attribute `"call"`.
#'
#' The class of the returned object is `c("evmissing", "gev_ts", "list")`.
#'
#' Objects inheriting from class `"evmissing"` have `coef`, `logLik`, `nobs`,
#' `summary`, `vcov`, `confint` and `plot` methods.  See [`evmissing_methods`].
#' @seealso [`gev_mle`] provides an adjustment for missing data in the
#' case where the raw data can be assumed to be independent and identically
#' distributed.
#'
#' @examples
#' set.seed(1632026)
#' blocks <- 50
#' block_length <- 90
#' missing_args <- list(p0miss = 0.5, min = 0, max = 0.4)
#' sdata <- sim_data(blocks = blocks, block_length = block_length,
#'                   missing_args = missing_args)
#'
#' pf <- gev_ts(sdata$data_miss, block_length = 90, pseudo = FALSE)
#' pt <- gev_ts(sdata$data_miss, block_length = 90, pseudo = TRUE)
#' pt2 <- gev_ts(sdata$data_miss, block_length = 90, pseudo = TRUE,
#'               seasonal = FALSE)

#' ptb <- gev_ts(sdata$data_miss, block = sdata$block, pseudo = TRUE)
#' ptb2 <- gev_ts(sdata$data_miss, block = sdata$block, pseudo = TRUE,
#'                seasonal = FALSE)
#'
#' ## La Plagne data
#'
#' # Precipitation
#' p <- gev_ts(LaPlagnePrecipMaximaSeasonal)
#' coef(p)
#'
#' # Snow depth
#' s <- gev_ts(LaPlagnePrecipMaximaSeasonal)
#' coef(s)
#' @export
gev_ts <- function(data, block_length, block, pseudo = TRUE, full = FALSE,
                   sliding = TRUE, seasonal = sliding, init = "quartiles",
                   ...) {
  # If data was created by block_maxima() or is a data frame that contains the
  # correct information then use it.
  # Otherwise, use block_maxima() to calculate the block maxima, the numbers
  # of non-missing values in the blocks and the largest possible number of
  # non-missing values in each block.
  from_maxima <- inherits(data, "block_maxima")
  required <- c("maxima", "notNA", "n", "whereNA", "pseudo_maxima",
                "full_maxima", "partial_maxima")
  has_components <- all(is.element(required, names(data)))
  if (from_maxima && has_components && inherits(data, "evmissing")) {
    maxima_notNA <- data
    no_NAs <- length(data$partial_maxima) == 0
  } else if (is.list(data)) {
    if (has_components) {
      maxima_notNA <- as.list(data)
    } else {
      stop("List ''data'' does not contain the required named components.")
    }
    no_NAs <- length(data$partial_maxima) == 0
  } else {
    # If data has no NAs then set pseudo = FALSE: pseudo maxima are irrelevant
    no_NAs <- !anyNA(data)
    if (no_NAs) {
      pseudo <- FALSE
    }
    message("Calculating block maxima")
    maxima_notNA <- block_maxima(data = data, block_length = block_length,
                                 block = block, pseudo = pseudo, full = full,
                                 sliding = sliding, seasonal = seasonal)
    message("Calculated block maxima")
  }
  # If there are maxima = NA, notNA = 0 entries in the data then remove them
  no_data <- which(maxima_notNA$notNA == 0)
  if (length(no_data) > 0) {
    maxima_notNA$maxima <- maxima_notNA$maxima[-no_data]
    maxima_notNA$notNA <- maxima_notNA$notNA[-no_data]
    maxima_notNA$n <- maxima_notNA$n[-no_data]
    maxima_notNA$partial_maxima <- maxima_notNA$partial_maxima[-no_data]
    maxima_notNA$pseudo_maxima <- maxima_notNA$pseudo_maxima[, -no_data]
  }
  # maxima_notNA is a list with 7 main components
  #         maxima: the block maxima, the response
  #          notNA: numbers of non-missing values in each block, the covariate
  #              n: largest possible number of non-missing values in each block
  #        whereNA: a list containing, for each block, the positions of missing
  #                 values in the block.
  #  pseudo_maxima: a numeric matrix of maxima created by applying missing value
  #                 patterns from incomplete blocks to full blocks. Each column
  #                 contains the maxima resulting from an incomplete block.
  #    full_maxima: maxima from full blocks.
  # partial_maxima: maxima from partial blocks.

  # If init is a has not been supplied then calculate initial estimates of mu
  # and sigma assuming that xi = 0
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
  message("Calling optim()")
  # If pseudo = FALSE or the data have no NAs then call the same log-likelihood
  # as gev_mle()
  if (!pseudo || no_NAs) {
    fit <- try(stats::optim(par = init, fn = negated_gev_loglik, hessian = TRUE,
                            ..., maxima_notNA = maxima_notNA, adjust = TRUE,
                            big_val = big_val), silent = TRUE)
  } else {
    fit <- try(stats::optim(par = init, fn = negated_gev_loglik_ts,
                            hessian = TRUE, ..., maxima_notNA = maxima_notNA,
                            pseudo = pseudo, big_val = big_val),
               silent = TRUE)
  }
  # If there is an error then it is probably because the sample size is small
  # and the log-likelihood does not have a maximum that is away from the
  # boundary of the parameter space. If this happens then return NA values for
  # the estimates etc.
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
    warning("MLEs not found: NA values returned")
    warning("optim() error: ", optim_error)
  } else {
    # Include the final values of r used in the fitting (if pseudo = TRUE)
    # Also a vector of the values of r for all blocks (=1 for a full block)
    if (pseudo && !no_NAs) {
      fit$rhats <- rhat(parameters = fit$par, maxima_notNA = maxima_notNA)
      rvec <- rep(1, length(maxima_notNA$maxima))
      names(rvec) <- names(maxima_notNA$maxima)
      rvec[names(fit$rhats)] <- fit$rhats[names(fit$rhats)]
      fit$rvec <- rvec
      # Call stats::optimHess(), supplying fixed_r, so that the values of r do
      # not vary as the GEV parameter values are changed
      hessian <- stats::optimHess(par = fit$par,
                                  fn = negated_gev_loglik_ts,
                                  maxima_notNA = maxima_notNA,
                                  pseudo = pseudo,
                                  fixed_r = as.numeric(fit$rhats),
                                  big_val = big_val)
    } else {
      fit$rhats <- NULL
      fit$rvec <- maxima_notNA$notNA / maxima_notNA$n
      rm_full <- which(is.element(names(fit$rvec),
                                  names(maxima_notNA$full_maxima)))
      fit$rhats <- fit$rvec[-rm_full]
      hessian <- fit$hessian
    }
    var_cov <- qr(hessian)
    if (var_cov$rank != ncol(var_cov$qr)) {
      fit$vcov <- matrix(NA, 3, 3)
      dimnames(fit$vcov) <- list(names(fit$par), names(fit$par))
      warning("The observed information matrix is singular.")
    } else {
      fit$vcov <- solve(hessian)
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
  fit$full_maxima <- maxima_notNA$full_maxima
  fit$partial_maxima <- maxima_notNA$partial_maxima
  fit$notNA <- maxima_notNA$notNA
  fit$n <- maxima_notNA$n
  fit$sliding <- sliding
  attr(fit, "call") <- match.call(expand.dots = TRUE)
  class(fit) <- c("evmissing", "gev_ts", class(fit))
  return(fit)
}
