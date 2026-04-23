#' Block Maxima and missing information
#'
#' Extracts disjoint or sliding block maxima and information relating to the
#' pattern of missing values. The (pseudo-)maxima created by applying blockwise
#' missing value patterns in incomplete block to full blocks can be calculated.
#' To be useful, the input data, `data` must contain at least one full block.
#'
#' @param data A numeric vector containing a time series of raw data.
#' @param block_length A numeric scalar. Used calculate the maxima of
#'   blocks of `block_length` contiguous values in the vector `data`.
#'   If `length(data)` is not an integer multiple of `block_length` then
#'   the values at the end of `data` that do not constitute a complete block
#'   of length `block_length` are discarded, without warning.
#' @param block A numeric vector with the same length as `data`. The value of
#'   `block[i]` indicates the block into which `data[i]` falls. For example,
#'   `block` could provide the year in which observation `i` was observed.
#' @param sliding A logical scalar. This is only allowed if `block_length` is
#'   supplied. If `sliding = TRUE` then maxima are calculated for **all**
#'   blocks of length `block_length`. See **Details**.
#' @param pseudo A logical scalar. If `pseudo = TRUE` then
#' @details Exactly one of the arguments `block_length` or `block` must be
#'   supplied.
#'
#'   If `block_length` is supplied and `sliding = TRUE` then the maxima are
#'   calculated for **all** blocks of length `block_length` present in `data`,
#'   starting with the first block `data[1:block_length]` and sliding the block
#'   repeatedly by one observation until reaching the final block
#'   `data[(length(data) - block_length + 1):length(data)]`.
#'
#' @return A list, with class `c("list", "block_maxima", "evmissing")`,
#'   containing the following numeric vectors:
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
#'  * `pseudo_maxima`: a numeric matrix containing (pseudo) block maxima
#'    created by applying the missing value patterns from incomplete blocks to
#'    all full blocks, that is, blocks without any missing values. Each column
#'    contains the pseudo-maxima resulting from a particular incomplete block.
#'    The columns are labelled by the number of the incomplete block and the
#'    columns by the number of the full block. If an incomplete block contains
#'    all missing values then its entry in `pseudo_maxima` is `NA`. If there
#'    are no full blocks or no incomplete blocks then `pseudo_maxima` is `NA`.
#'  * `full_maxima`: a numeric vector of maxima from full blocks.
#'  * `partial_maxima`: a numeric vector of maxima from partial blocks.
#'
#' The input arguments `sliding` and `pseudo` are also included.
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
#' ## Data with an incomplete block
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
block_maxima <- function(data, block_length, block, sliding = FALSE,
                         pseudo = FALSE) {
  # Check that data is a numeric vector
  if (!is.numeric(data)) {
    stop("''data'' must be a numeric vector.")
  }
  block_length_supplied <- !missing(block_length)
  block_supplied <- !missing(block)
  # If sliding = TRUE then check that block_length is supplied
  if (block_supplied && sliding) {
    stop("''sliding = TRUE'' is only allowed if ''block_length'' is supplied.")
  }
  # Check that exactly one of block_length or block has been supplied
  if (block_length_supplied & block_supplied) {
    stop("Only one of ''block_length'' or ''block'' may be supplied.")
  }
  if (!block_length_supplied && !block_supplied) {
    stop("''block_length'' or ''block'' must be supplied.")
  }
  # Find block maxima depending on whether block_length or block is supplied
  if (block_length_supplied) {
    # Set the number of blocks and function to select data from a block
    if (sliding) {
      number_of_blocks <- length(data) - block_length + 1
      find_block_i <- function(i) {
        return(data[i:(i + block_length - 1)])
      }
    } else {
      number_of_blocks <- floor(length(data) / block_length)
      find_block_i <- function(i) {
        return(data[(1 + block_length * (i - 1)):(block_length * i)])
      }
    }
    # Function to calculate the block maxima
    find_maxima_notNA <- function(i) {
      temp <- find_block_i(i)
      not_na <- !is.na(temp)
      if (any(not_na)) {
         val <- c(max(temp, na.rm = TRUE), sum(not_na))
      } else {
         val <- c(NA, 0)
      }
      return(val)
    }
    r <- vapply(seq_len(number_of_blocks), find_maxima_notNA, c(0.0, 0.0))
    r <- list(maxima = r[1, ], notNA = r[2, ],
              n = rep_len(block_length, number_of_blocks))
  } else {
    if (length(block) != length(data)) {
      stop("''block'' must have the same length as ''data''.")
    }
    FUN <- function(data) {
      not_na <- !is.na(data)
      if (any(not_na)) {
        val <- c(max(data, na.rm = TRUE), sum(not_na), length(data))
      } else {
        val <- c(NA, 0, length(data))
      }
      return(val)
    }
    r <- tapply(data, block, FUN = FUN)
    r <- .mapply(c, r, NULL)
    names(r) <- c("maxima", "notNA", "n")
  }
  # Add the input arguments sliding and pseudo
  r$sliding <- sliding
  r$pseudo <- pseudo
  # Give the returned object a class, so that we can detect block maxima data
  # created by block_maxima()
  class(r) <- c("list", "block_maxima", "evmissing")
  return(r)
}
