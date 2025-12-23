#' Return Level Inferences
#'
#' Calculates point estimates of \eqn{m}-year return levels for fitted model
#' objects returned from [`gev_mle`].
#'
#' @param x An object inheriting from class `evmissing` returned from
#'   [`gev_mle`].
#' @param m A numeric vector.  Values of \eqn{m}, the return periods of
#'   interest, in years.
#' @param npy A numeric scalar.  The number \eqn{n_{py}} of block maxima per
#'   year. If the blocks are of length 1 year then `npy = 1`.
#' @details For \eqn{\xi \neq 0}, the \eqn{m}-year return level is given by
#'   \eqn{z_m = \mu + \sigma (y_p ^ {-\xi} - 1) / \xi}, where
#'   \eqn{y_p = -\log(1 - p)} and \eqn{p = 1 - (1 - 1 / m) ^ {1 / n_{py}}}.
#'   For \eqn{\xi = 0}, \eqn{z_m = \mu - \sigma \log y_p}.
#'   Equivalently, we could note that \eqn{z_m = \mu - \sigma BC(y_p, -\xi)},
#'   where \eqn{BC(x, \lambda)} is a Box-Cox transformation.
#' @return An object with class `c("return_level", "numeric", "evmissing")`.
#'   A numeric vector containing the MLEs of the required return levels, with
#'   names indicating the return period. The fitted model object returned
#'   from [`gev_mle`] is included as an attribute called `"gev_mle"`.
#'   The input arguments `m` and `npy` are also included as attributes as is
#'   the call to `gev_return`.
#' @references Coles, S. G. (2001) *An Introduction to Statistical
#'   Modeling of Extreme Values*, Springer-Verlag, London.
#'   \doi{10.1007/978-1-4471-3675-0_3}
#' @seealso [`return_level_methods`] for `print, summary, coef, vcov` and
#'   `confint` methods.
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
#' # Make adjustment for the numbers of non-missing values per block
#' fit2 <- gev_mle(sdata$data_miss, block_length = sdata$block_length)
#' summary(fit2)
#'
#' gev_return(fit1, m = c(100, 1000))
#' gev_return(fit2, m = c(100, 1000))
#'
#' ## Plymouth ozone data
#'
#' fit <- gev_mle(PlymouthOzoneMaxima)
#' rl <- gev_return(fit, m = c(100, 200))
#'
#' # Symmetric confidence intervals
#' sym <- confint(rl)
#'
#' # Profile-based confidence intervals
#'
#' prof <- confint(rl, profile = TRUE)
#' prof
#' plot(prof, digits = 4)
#' plot(prof, parm = 2, digits = 3)
#'
#' # Doing this more quickly when we only care about the confidence limits
#' prof <- confint(rl, profile = TRUE, mult = 32, faster = TRUE)
#' plot(prof, digits = 3, type = "b")
#' plot(prof, parm = 2, digits = 3, type = "b")
#' @export
gev_return <- function(x, m = 100, npy = 1) {
  # Extract the MLEs of the GEV parameters
  mles <- coef(x)
  mu <- mles["mu"]
  sigma <- mles["sigma"]
  xi <- mles["xi"]
  # Set the annual probability of exceedance based on m and npy
  p <- 1 - (1 - 1 / m) ^ (1 / npy)
  return_levels <- nieve::qGEV(p, loc = mu, scale = sigma, shape = xi,
                               lower.tail = FALSE)
  names(return_levels) <- paste0(round(m, 2), "-year level")
  attr(return_levels, "gev_mle") <- x
  attr(return_levels, "m") <- m
  attr(return_levels, "npy") <- npy
  attr(return_levels, "call") <- match.call(expand.dots = TRUE)
  class(return_levels) <- c("return_level", class(return_levels), "evmissing")
  return(return_levels)
}
