#' Block maxima
#'
#' Extracts block maxima and the number of non-missing observations per block.
#'
#' @param data A numeric vector containing a time series of raw data.
#' @param block_length A numeric scalar. Used calculate the maxima of disjoint
#'   blocks of `block_length` contiguous values in the vector `data`.
#'   If `length(data)` is not an integer multiple of `block_length` then
#'   the values at the end of `data` that do not constitute a complete block
#'   of length `block_length` are discarded, without warning.
#' @param block A numeric vector with the same length as `data`. The value of
#'   `block[i]` indicates the block into which `data[i]` falls. For example,
#'   `block` could provide the year in which observation `i` was observed.
#' @details Exactly one of the arguments `block_length` or `block` must be
#'   supplied.
#'
#' @return A list, with class `c("list", "block_maxima", "evmissing")`,
#'   containing the following numeric vectors:
#'
#'  * `maxima`: the block maxima.
#'  * `notNA`: the numbers of non-missing observations in each block.
#'  * `n`: the maximal block length, that is, the largest number of values that
#'     could have been observed in each block.
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
#' # Supplying block_length
#' block_length <- 3
#' block_maxima(data, block_length = block_length)
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
#' @export
block_maxima <- function(data, block_length, block) {
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
    number_of_blocks <- floor(length(data) / block_length)
    find_maxima_notNA <- function(i) {
      temp <- data[(1 + block_length * (i - 1)):(block_length * i)]
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
  # Give the returned object a class, so that we can detect block maxima data
  # created by block_maxima()
  class(r) <- c("list", "block_maxima", "evmissing")
  return(r)
}
