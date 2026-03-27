#' GEV ML Inference with Adjustment for Missing Data (Stationary Sequences)
#'
#' Fits a GEV distribution to block maxima using maximum likelihood estimation,
#' making an adjustment for the locations of missing raw values in each block.
#' The GEV location and scale parameters are adjusted to reflect the proportion
#' of raw values that are missing and the time series dependence in the data.
#'
#' @param data Either
#'
#'   * a numeric vector containing a time series of raw data,
#'   * an object returned from [`block_maxima_ts`] or
#'     [`sliding_block_maxima_ts`], a list with components `maxima`, `notNA`,
#'     `n`, `whereNA`, `pseudo_maxima`, `full_maxima` and `partial_maxima`,
#'   * a named list containing the same information, that is, the variables
#'     `maxima`, `notNA`, `n`, `whereNA`, `pseudo_maxima`, `full_maxima` and
#'     `partial_maxima` as an object returned from [`block_maxima_ts`] or
#'     [`sliding_block_maxima_ts`].
#'
#'   There must be at least one full block of data, that is, at least one block
#'   for which no data are missing.
#' @param block_length A numeric scalar. Used calculate the maxima of disjoint
#'   blocks of `block_length` contiguous values in the vector `data`.
#'   If `sliding = FALSE` and if `length(data)` is not an integer multiple of
#'   `block_length`, then the values at the end of `data` that do not constitute
#'   a complete block of length `block_length` are discarded, without warning.
#' @param block A numeric vector with the same length as `data`. The value of
#'   `block[i]` indicates the block into which `data[i]` falls. For example,
#'   `block` could provide the year in which observation `i` was observed.
#'   Not applicable if `sliding = TRUE`. If `sliding = TRUE`, then
#'   `block_length` must be supplied.
#' @param pseudo A logical scalar. If `pseudo = TRUE` then the pseudo-maxima
#'   returned from [`block_maxima_ts`] are used to estimate the value of
#'   \eqn{r_i} for an incomplete, partially-observed block. See **Details**.
#'   If `pseudo = FALSE` then the ratio \eqn{n_i/n} is used, as in
#'   [`gev_mle()`].
#' @param sliding A logical scalar. If `sliding = TRUE` then inferences are
#'   based on sliding block maxima returned by [`sliding_block_maxima_ts`]
#'   and `block_length` must be supplied. If `sliding = FALSE` then they are
#'   based on disjoint block maxima returned from [`block_maxima_ts`].
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
#'   `block_length` or `block` must be supplied if `sliding = FALSE` and only
#'   `block_length` can be supplied if `sliding = TRUE`. The parameters are
#'   fitted using maximum likelihood estimation.
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
#' For a full block, \eqn{r_i = 1}. If `pseudo = TRUE`, then, for an
#' incomplete, partially-observed block, the value of \eqn{r_i} is estimated
#' using the pseudo-maxima returned from [`block_maxima_ts`] and the GEV
#' distribution function based on the current value of \eqn{(\mu, \sigma, \xi)}
#' in the optimisation routine. Suppose that we have a vector \eqn{M_i} of
#' pseudo-maxima resulting from a particular incomplete block \eqn{i}.  It can
#' be shown that the components of \eqn{V_i = -\log G(M_i)} each have an
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
#' * `notNA`: the number of non-missing raw values on which the maxima in
#'   `maxima` are based.
#' * `n`: the maximal block length, that is, the largest number of values that
#'     could have been observed in each of these blocks.
#' * `rvec`: a vector of the values used for \eqn{r_1, ..., r_b}, where \eqn{b}
#'   is the number of blocks. The content depends on the argument `pseudo`.
#' * `rhats`: if `pseudo = TRUE`, a vector of the subset of `rvec` for
#'    partially-observed blocks only. The attributes `"propn_notNA"` and
#'    `"unconstrained"` give, respectively, the values of \eqn{n_i/n} for these
#'    blocks and the estimates of \eqn{r_i} before they are constrained to lie
#'    in the interval \eqn{[n_i/n, 1]}.
#' * `sliding`: the input argument `sliding`.
#'
#' The call to `gev_ts` is provided in the attribute `"call"`.
#'
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
#' # Using disjoint blocks
#' pt <- gev_ts(sdata$data_miss, block_length = 90, pseudo = TRUE)
#' pf <- gev_ts(sdata$data_miss, block_length = 90, pseudo = FALSE)
#' pf2 <- gev_mle(sdata$data_miss, block_length = 90)
#'
#' # Using sliding blocks
#' pts <- gev_ts(sdata$data_miss, block_length = 90, pseudo = TRUE,
#'               sliding = TRUE)
#' pfs <- gev_ts(sdata$data_miss, block_length = 90, pseudo = FALSE,
#'               sliding = TRUE)
#' @export
gev_ts <- function(data, block_length, block, pseudo = TRUE, sliding = FALSE,
                   init = "quartiles", ...) {
  # If sliding = TRUE then check that only block_length is supplied
  block_length_supplied <- !missing(block_length)
  block_supplied <- !missing(block)
  if (sliding && (block_supplied || !block_length_supplied)) {
    stop("If ''sliding = TRUE'' then (only) ''block_length'' must be supplied")
  }
  # If data was created by block_maxima_ts() or sliding_block_maxima_ts() or is
  # a data frame that contains the correct information then use it.
  # Otherwise, use block_maxima_ts() to calculate the block maxima, the numbers
  # of non-missing values in the blocks and the largest possible number of
  # non-missing values in each block.
  from_maxima_ts <- inherits(data, "block_maxima_ts") ||
    inherits(data, "sliding_block_maxima_ts")
  required <- c("maxima", "notNA", "n", "whereNA", "pseudo_maxima",
                "full_maxima", "partial_maxima")
  if (from_maxima_ts && inherits(data, "evmissing")) {
    maxima_notNA <- data
  } else if (is.list(data)) {
    if (all(is.element(required, names(data)))) {
      maxima_notNA <- as.list(data)
    } else {
      stop("List ''data'' does not contain the required named components.")
    }
  } else {
    if (sliding) {
      print("Calculating sliding block maxima")
      maxima_notNA <- sliding_block_maxima_ts(data, block_length)
      print("Calculated sliding block maxima")
    } else {
      maxima_notNA <- block_maxima_ts(data, block_length, block)
    }
  }
  # If there are maxima = NA, notNA = 0 entries in the data then remove them
  no_data <- which(maxima_notNA$notNA == 0)
  if (length(no_data) > 0) {
    maxima_notNA <- lapply(maxima_notNA, function(x) x[-no_data])
  }
  # maxima_notNA is a list with 7 components
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
  print("Calling optim()")
  fit <- try(stats::optim(par = init, fn = negated_gev_loglik_ts,
                          hessian = TRUE, ..., maxima_notNA = maxima_notNA,
                          pseudo = pseudo, big_val = big_val),
             silent = TRUE)
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
    warning("MLEs not found: NA values returned.")
  } else {
    # Include the final values of r used in the fitting (if pseudo = TRUE)
    # Also a vector of the values of r for all blocks (=1 for a full block)
    if (pseudo) {
      fit$rhats <- rhat(parameters = fit$par, maxima_notNA = maxima_notNA)
      rvec <- rep(1, length(maxima_notNA$maxima))
      rvec[as.numeric(names(fit$rhats))] <- fit$rhats
      fit$rvec <- rvec
      # Call stats::optimHess(), suppling fixed_r, so that the values of r do
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
  fit$notNA <- maxima_notNA$notNA
  fit$n <- maxima_notNA$n
  fit$sliding <- sliding
  attr(fit, "call") <- match.call(expand.dots = TRUE)
  class(fit) <- c("evmissing", "mle", class(fit))
  return(fit)
}
