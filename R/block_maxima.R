#' Block Maxima and Missing Information
#'
#' Extracts disjoint block maxima from a time series of raw data. Can also
#' provide information about the effect of missing value patterns on block
#' maxima via pseudo-maxima created by applying blockwise missing value
#' patterns in partially-observed (partial) blocks to fully-observed (full)
#' blocks.
#'
#' @param data A numeric vector containing a time series of raw data. If
#'   `sliding = FALSE` then `data` must contain at least one full disjoint
#'   block. If `sliding = TRUE` then `data` must contain at least one full
#'   sliding block.
#' @param block_length A numeric scalar. Used calculate the maxima of
#'   disjoint blocks of `block_length` contiguous values in the vector `data`.
#'   If `length(data)` is not an integer multiple of `block_length` then
#'   the values at the end of `data` that do not constitute a complete block
#'   of length `block_length` are discarded, without warning.
#' @param block A numeric vector with the same length as `data`. The value of
#'   `block[i]` indicates the block into which `data[i]` falls. For example,
#'   `block` could provide the year in which observation `i` was observed.
#'   The block lengths implied by `block` should have similar values, for
#'   example, 366 for leap years and 365 for other years.
#' @param pseudo A logical scalar. If `pseudo = TRUE` then pseudo-maxima are
#'   calculated as the block maxima obtained by applying the missing value
#'   patterns from partial blocks to all full blocks.
#' @param sliding A logical scalar. Only relevant if `pseudo = TRUE`. If
#'   `sliding = TRUE` then pseudo-maxima are calculated for **all** full blocks
#'   of the relevant length, rather than only a set of disjoint blocks.
#'   See **Details**.
#' @param season A logical scalar. Only relevant if `pseudo = TRUE` and
#'   `sliding = TRUE`. If `season = TRUE` then the way in which the
#'   pseudo-maxima are calculated respects the seasonality that may be
#'   exhibited over the duration of a block. If, for example, a block covers a
#'   single year, then the missing values applied to a full block occur at the
#'   same time of year as in the originating partial block. If `season = FALSE`
#'   then the missing values applied have the same positions in the partial and
#'   full blocks with respect to the start of each block.
#' @details Exactly one of the arguments `block_length` or `block` must be
#'   supplied.
#'
#'   If `block_length` is supplied and `sliding = TRUE` then the maxima are
#'   calculated for **all** blocks of length `block_length` present in `data`,
#'   starting with the first block `data[1:block_length]` and sliding the block
#'   repeatedly by one observation until reaching the final block
#'   `data[(length(data) - block_length + 1):length(data)]`.
#'
#'   Also explain for `block`
#'
#' @return A list, with class
#'   `c("list", "block_maxima", "disjoint", "evmissing")`, containing the
#'   following numeric vectors:
#'
#'  * `maxima`: the block maxima.
#'  * `notNA`: the numbers of non-missing observations in each block.
#'  * `n`: the maximal block length, that is, the largest number of values that
#'     could have been observed in each block.
#'
#' If `pseudo = TRUE` then the returned list also contains the following:
#'
#'  * `whereNA`: a named list containing, for each block, the positions of any
#'    missing values in the block. For example, if (only) the first and fifth
#'    observations in block 3 are missing then the third component (named
#'    `block3`) of `whereNA` is `c(1, 5)`. If a block has no missing values
#'    then its component in `whereNA` is `integer(0)`.
#'  * `pseudo_maxima`: a numeric matrix containing block maxima created by
#'    applying the missing value patterns from partial blocks to all full
#'    blocks. Each column contains the pseudo-maxima resulting from a
#'    particular partial block. The columns are labelled by the number of the
#'    partial block and the columns by the number of the full block. If a
#'    partial block contains all missing values then its entry in
#'    `pseudo_maxima` is `NA`. If there are no full blocks or no partial blocks
#'    then `pseudo_maxima` is `NA`.
#'  * `full_maxima`: a numeric vector of maxima from full blocks.
#'  * `partial_maxima`: a numeric vector of maxima from partial blocks.
#'
#' The input arguments `pseudo` and `sliding` are also included.
#'
#' If a block contains only missing values then its value of `maxima` is `NA`
#' and its value of `notNA` is `0`.
#'
#' If `block` is supplied then these vectors are named using the values in
#' `block`. Otherwise, these vectors do not have names.
#' @examples
#' ## Simulate example data
#' set.seed(7032025)
#' data <- rexp(15)
#'
#' # Create some missing values
#' data[c(5, 7:8)] <- NA
#' # 5 blocks (columns), each with 3 observations
#' matrix(data, ncol = 5)
#'
#' # Supplying block_length, disjoint maxima
#' block_length <- 3
#' block_maxima(data, block_length = block_length)
#' # Supplying block_length, sliding maxima
#' block_maxima(data, block_length = block_length, sliding = TRUE)
#'
#' # Supplying block
#' block <- rep(1:5, each = 3)
#' block_maxima(data, block = block)
#'
#' ## Data with a partially-observed block
#' data <- c(data, 1:2)
#'
#' # Supplying block_length (the extra 2 observations are ignored)
#' block_length <- 3
#' block_maxima(data, block_length = block_length)
#' # Supplying block (with an extra group indicator)
#' block <- c(block, 7, 7)
#' block_maxima(data, block = block)
#'
#' @seealso Plot method [`plot.block_maxima`].
#' @export
block_maxima <- function(data, block_length, block, pseudo = FALSE,
                         sliding = FALSE, season = FALSE) {
  # Check that data is a numeric vector
  if (!is.numeric(data)) {
    stop("''data'' must be a numeric vector.")
  }
  block_length_supplied <- !missing(block_length)
  block_supplied <- !missing(block)
  # Check that exactly one of block_length or block has been supplied
  if (block_length_supplied & block_supplied) {
    stop("Only one of ''block_length'' or ''block'' may be supplied.")
  }
  if (!block_length_supplied && !block_supplied) {
    stop("''block_length'' or ''block'' must be supplied.")
  }
  # Find block maxima depending on whether block_length or block is supplied
  if (block_length_supplied) {
    r <- maxima_block_length(data = data, block_length = block_length,
                             pseudo = pseudo)
  } else {
    r <- maxima_block(data = data, block = block,
                      pseudo = pseudo)
  }
  # If pseudo = TRUE then, for each partial block, apply its missing value
  # pattern to each full block and calculate the resulting pseudo maxima
  if (pseudo) {
    if (block_length_supplied) {
      pseudo_maxima <- pseudo_maxima_block_length(maxima_notNA = r, data = data,
                                                  block_length = block_length)
    } else {
      pseudo_maxima <- pseudo_maxima_block(maxima_notNA = r, data = data,
                                           block = block)
    }
    r <- c(r, list(pseudo_maxima = pseudo_maxima))
    # Create vectors that contain the full maxima and partial maxima
    r$full_maxima <- r$maxima[as.numeric(rownames(r$pseudo_maxima))]
    r$partial_maxima <- r$maxima[as.numeric(colnames(r$pseudo_maxima))]
  }
  # Add the input arguments sliding and pseudo
  r$pseudo <- pseudo
  r$sliding <- sliding
  # Give the returned object a class, so that we can detect block maxima data
  # created by block_maxima()
  class(r) <- c("list", "block_maxima", "disjoint", "evmissing")
  return(r)
}
