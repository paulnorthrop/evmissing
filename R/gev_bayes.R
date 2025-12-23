#' GEV Bayesian Inference with Adjustment for Missing Data
#'
#' Performs Bayesian inference using a GEV distribution using block maxima,
#' with the option to make an adjustment for the numbers of non-missing raw
#' values in each block.
#'
#' @inheritParams gev_mle
#' @param prior Specifies a prior distribution for the GEV parameters. This is
#'   most easily set using [`revdbayes::set_prior`]. The default is a prior
#'   \eqn{\pi(\mu, \sigma, \xi) \propto \sigma^{-1}} for \eqn{\sigma > 0}. See
#'   [`revdbayes::set_prior`] for details.
#' @param n A non-negative integer. The number of values to simulate from the
#'   posterior distribution for \eqn{(\mu, \sigma, \xi)}.
#' @param ... Further arguments to be passed to [`rust::ru`].
#' @details The likelihood described in [`gev_mle`] is combined with the prior
#'   density provided by `prior` to produce, up to proportionality, a
#'   posterior density for \eqn{(\mu, \sigma, \xi)}.
#'
#'   A function to evaluate the log-posterior is passed to [`rust::ru`] to
#'   simulate a random sample from this posterior distribution using the
#'   generalised ratio-of-uniforms method, using relocation of the mode of the
#'   density to the origin to increase efficiency. The value of `init` is used
#'   as an initial estimate in a search for the posterior mode. Arguments to
#'   [`rust::ru`] can be passed via `...`. The default setting is
#'   `trans = "none"`, that is, no transformation of the margins, and
#'   `rotate = TRUE`, rotation of the parameter axes to improve isotropy
#'   with a view to increasing efficiency.
#' @return An object returned from [`rust::ru`]. The following components are
#'   added to this list
#'
#' * `model`: = `"gev"`.
#' * `data,prior`: the inputs `data` and `prior`.
#' * `call`: the call to `gev_bayes`.
#' * `maxima`: the vector of block maxima used to fit the model.
#' * `notNA`: the number of non-missing raw values on which the maxima in
#'   `maxima` are based.
#' * `n`: the maximal block length, that is, the largest number of values that
#'    could have been observed in each of these blocks.
#' * `adjust`: a logical scalar indicating whether or not the adjustment in
#'   the **Details** section of [`gev_mle`] was performed. This is `TRUE`
#'   only if the input argument `adjust` was `TRUE`.
#' * `adjust,discard` : the values of these input arguments.
#'
#' The class of the returned object is
#' `c("evpost", "ru", "bayes", "evmissing")`.
#' Objects of class `"evpost"` have [`print`][revdbayes::print.evpost],
#' [`summary`][revdbayes::summary.evpost] and [`plot`][revdbayes::plot.evpost]
#' S3 methods.
#' @examples
#' ## Simulate data with missing values
#'
#' set.seed(24032025)
#' blocks <- 50
#' block_length <- 365
#'
#' # Simulate raw data from an exponential distribution
#' sdata <- sim_data(blocks = blocks, block_length = block_length)
#'
#' block_length <- sdata$block_length
#' # Sample from the posterior based on block maxima from full data
#' post1 <- gev_bayes(sdata$data_full, block_length = block_length)
#' summary(post1)
#'
#' # Sample with adjustment for the number of non-missing values per block
#' post2 <- gev_bayes(sdata$data_miss, block_length = block_length)
#' summary(post2)
#'
#' # Do not make the adjustment
#' post3 <- gev_bayes(sdata$data_miss, block_length = block_length,
#'                    adjust = FALSE)
#' summary(post3)
#'
#' # Remove all block maxima with greater than 25% missing values and
#' # do not make the adjustment
#' post4 <- gev_bayes(sdata$data_miss, block_length = block_length,
#'                    adjust = FALSE, discard = 25)
#' summary(post4)
#'
#' ## Brest sea surge data
#'
#' post <- gev_bayes(BrestSurgeMaxima)
#' summary(post)
#' plot(post)
#' @export
gev_bayes <- function(data, block_length, block, adjust = TRUE, discard = 0,
                      init = "quartiles",
                      prior = revdbayes::set_prior(prior = "flat",
                                                   model = "gev"),
                      n = 1000, ...) {
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

  # Sample from a posterior distribution for a model in which block maximum i
  # of n has a GEV distribution with
  #     location: mu + sigma * [(n_i / n) ^ xi - 1] / xi,
  #   log(scale): log(sigma) + xi * log (n_i / n),
  #        shape: xi,
  #   where n is the block length, n_i is the number of non-missing values in
  #   block i and the block maxima are conditionally independent given
  #   (n_1, ..., n_b).

  # If init is a has not been supplied then calculate initial estimates of mu and
  # sigma for assuming that xi = 0
  if (is.character(init)) {
    init_method <- match.arg(init, c("quartiles", "moments"))
    init <- gev_init(maxima_notNA, init_method = init_method)
  } else {
    names(init) <- c("mu", "sigma", "xi")
  }
  # Create a function to return the log-posterior
  logpost <- function(parameters, log_post_args) {
    loglik <- -do.call(negated_gev_loglik ,
                       c(list(parameters = parameters), logpost_args))
    if (is.infinite(loglik)) return(loglik)
    logprior <- do.call(prior$prior, c(list(parameters), prior[-1]))
    return(loglik + logprior)
  }
  # Arguments to logpost, in addition to the parameters
  logpost_args <- list(maxima_notNA = maxima_notNA, adjust = adjust)
  # Extract min_xi and max_xi from prior (if supplied)
  min_xi <- ifelse(is.null(prior$min_xi), -Inf, prior$min_xi)
  max_xi <- ifelse(is.null(prior$max_xi), +Inf, prior$max_xi)
  # Extract arguments to be passed to ru()
  ru_args <- list(...)
  # Set default values for trans and rotate if they have not been supplied.
  # This isn't necessary, because these are the current defaults in rust::ru(),
  # but it can be useful to be explicit.
  if (is.null(ru_args$trans)) ru_args$trans <- "none"
  if (is.null(ru_args$rotate)) ru_args$rotate <- TRUE
  # Create list of objects to send to function ru()
  # d is the number of para,=meters
  # lower and upper gives the respective supports for mu, sigma and xi
  # var_names sets variable (parameter) names
  fr <- list(d = 3, lower = c(-Inf, 0, min_xi), upper = c(Inf, Inf, max_xi),
             var_names = c("mu","sigma", "xi"))
  # Create a list of argument for rust::ru()
  for_ru <- c(list(logf = logpost), fr, list(init = init, n = n),
              ru_args)
  post <- do.call(rust::ru, for_ru)
  # Add the components that are in an object returned from revdbayes::rpost()
  post$model <- "gev"
  post$data <- data
  post$prior <- prior
  # Replace the (verbose) call to rust::ru() with the call to gev_bayes()
  post$call <- match.call(expand.dots = TRUE)
  post$maxima <- maxima_notNA$maxima
  post$notNA <- maxima_notNA$notNA
  post$n <- maxima_notNA$n
  post$adjust <- adjust
  post$discard <- discard
  class(post) <- c("evpost", class(post), "bayes", "evmissing")
  return(post)
}
