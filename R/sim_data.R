#' Simulate raw data
#'
#' Simulates data from a user-supplied distribution and creates missing values
#' artificially. Functions `mcar` and `mcar2` provides an example mechanisms
#' for doing this based on a Missing Completely At Random (MCAR) assumption.
#'
#' @param blocks A numeric scalar. The number of blocks of data required.
#'   Usually, this will be a positive integer, but `blocks = 0` returns a list
#'   containing in the input arguments, in particular, `distn`, `distn_args`
#'   and `block_length`. This feature is provided so that a simulation setup
#'   could be replicated in without actually simulating data.
#' @param block_length A numeric scalar. The number of raw observations per
#'   block.
#' @param distn A character scalar. Specifies the distribution from which raw
#'   data are simulated. The name in the `xxx` part of the `dxxx, pxxx, qxxx`
#'   and `rxxx` distributional functions in the `stats` package. See
#'   [`stats::Distributions`].
#' @param missing_fn A function to simulate the positions of the missing values
#'   within each block year. See **Details**.
#' @param missing_args Arguments to be passed to `missing_fn`. If `missing_fn`
#'   is `mcar` then a subset of `p0miss`, `min` and `max` may be supplied
#'   in the list `missing_args`. The values of the remaining components will be
#'   set at their default values.
#' @param ... Further arguments to the function `stats::rxxx`. The argument `n`
#'   is set within `sim_data` to be equal to `block_length * blocks`.
#' @param sim_data A numeric vector of raw observations into, some of which
#'   will be made missing.
#'
#' @details The function `missing_fn` must return a, possibly empty,
#'   subset of `c(1, 2, ..., block_length)`. This function is applied within
#'   each simulated block, independently of other blocks.
#'
#'   The default function `mcar` simulates the numbers of missing values in the
#'   blocks as follows.
#'
#'   * A proportion `p0miss` of the blocks have **no** missing values.
#'   * In the other blocks, the number of missing values is
#'     `ceiling(prop_miss * block_length)`, where `prob_miss` is a value
#'     simulated from a Uniform(`min`, `max`) distribution. The positions of
#'     these missing values within the block is random.
#'
#'   The function `mcar2` identifies at random a proportion `pmiss` of the
#'   simulated raw observations to become missing.
#'
#'   Care may need to be taken if these simulated data are used as input to
#'   [`gev_mle`] using an approach that discards block maxima based on more
#'   than a certain percentage of  missing values, that is, with `discard > 0`.
#'   For example, using the default argument `blocks = 50` and
#'   `missing_fn = mcar`, with its default `missing_args`, may result
#'   in a sample size of retained block maxima that contains insufficient
#'   information to make reliable inferences, leading to difficulties finding
#'   an appropriate MLE for the shape parameter \eqn{\xi} and/or a singular
#'   observed information matrix.
#'
#' @return If `blocks > 0`, a list with the following components:
#'
#' * `data_full`: simulated raw data with no missing values.
#' * `data_miss`: simulated data after missing values have been created.
#' * `blocks, block_length`: the respective input values of `blocks` and
#'   `block_length`.
#' * `block`: a block indicator vector, suitable as an argument to [`gev_mle`].
#' * `distn`: the input argument `distn`.
#' * `distn_args`: further arguments to `stats::rxxx` supplied via `...`.
#'
#' If `blocks = 0`, a list containing all the inputs arguments.
#' @seealso [`gev_mle`] and [`gev_bayes`].
#' @examples
#' # Using missing_fn = mcar
#' sdata <- sim_data()
#'
#' # Using missing_fn = mcar2
#' sdata2 <- sim_data(missing_fn = mcar2)
#' @name sim_data
NULL
## NULL

#' @rdname sim_data
#' @export
sim_data <- function(blocks = 50, block_length = 365, distn = "exp",
                     missing_fn = mcar,
                     missing_args = formals(missing_fn)$missing_args, ...) {
  # If missing_args has been supplied then fill in missing default arguments
  if (!missing(missing_args)) {
    default_args <- as.list(formals(missing_fn)$missing_args)[-1]
    missing_args <- merge_two_lists(missing_args, default_args)
  }
  if (blocks > 0) {
    # Simulate raw data with no missings
    # Add to sim_args the block length and any arguments in ...
    sim_args <- c(list(n = block_length * blocks), list(...))
    sim_fn <- paste0("r", distn)
    data_full <- do.call(sim_fn, sim_args)
    # Create missing values in data_full
    for_missing_fn <- list(sim_data = data_full, blocks = blocks,
                           block_length = block_length,
                           missing_args = missing_args)
    # Retain only components that are formal arguments of missing_fn()
    for_missing_fn <- for_missing_fn[names(formals(missing_fn))]
    # Call missing_fn()
    data_miss <- do.call(missing_fn, for_missing_fn)
    val <- list(data_full = data_full, data_miss = data_miss,
                blocks = blocks, block_length = block_length,
                block = rep(c(1:blocks), each = block_length),
                distn = distn, distn_args = list(...))
  } else {
    val <- list(blocks = blocks, block_length = block_length, distn = distn,
                distn_args = list(...), missing_fn = missing_fn,
                missing_args = missing_args)
  }
  return(val)
}

#' @rdname sim_data
#' @export
mcar <- function(sim_data, blocks, block_length,
                 missing_args = list(p0miss = 0, min = 0, max = 0.5)) {
  # Add to missing_args the n = 1
  missing_args <- c(list(n = 1), missing_args)
  # A function to return the locations of the missing values in a block
  block_fun <- function(block) {
    prop_miss <- do.call(stats::runif, missing_args[-2])
    num_miss <- ceiling(prop_miss * block_length)
    block_miss <- sample(block_length, num_miss, replace = FALSE)
    missing_values <- block_length * (block - 1) + block_miss
    return(missing_values)
  }
  x <- sapply(1:blocks, block_fun)
  # Select a proportion p0miss of the blocks to revert back no missings
  n0 <- round(missing_args$p0miss * blocks)
  revert_to_no_missings <- sample(x = blocks, size = n0)
  x[revert_to_no_missings] <- list(numeric(0))
  # Create the missing values
  sim_data[unlist(x)] <- NA
  return(sim_data)
}

#' @rdname sim_data
#' @export
mcar2 <- function(sim_data, missing_args = list(pmiss = 0.5)) {
  nsim <- length(sim_data)
  nmiss <- round(missing_args$pmiss * nsim)
  missing <- sample(nsim, nmiss)
  sim_data[missing] <- NA
  return(sim_data)
}
