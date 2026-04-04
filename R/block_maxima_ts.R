#' Block Maxima for a Time Series
#'
#' Extracts block maxima and missing value information for each block.
#' Works like [`block_maxima`] but returns extra components, including:
#' `whereNA`, the positions of the missing values within each block, and
#' `pseudo_maxima`, the maxima created by applying blockwise missing value
#' patterns in incomplete blocks to full blocks, that is, blocks without any
#' missing values. To be useful, the input data, `data`, must contain at least
#' one full block.
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
#'   supplied. If the block sizes implied by `block` are unequal then an
#'   incomplete block and a full block may have different lengths. If this
#'   occurs when `pseudo_maxima` are calculated, then the longer block is
#'   trimmed, by discarding trailing values, so that the lengths match.
#'
#' @return A list, with class `c("list", "block_maxima_ts", "evmissing")`,
#'   containing the following components:
#'
#'  * `maxima`: the block maxima.
#'  * `notNA`: the numbers of non-missing observations in each block.
#'  * `n`: the maximal block length, that is, the largest number of values that
#'     could have been observed in each block.
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
#' If a block contains only missing values then its value of `maxima` is `NA`,
#' its value of `notNA` is `0` and `whereNA` contains the positions of all the
#' observations in the block.
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
#' block_maxima_ts(data, block_length = block_length)
#' # Supplying block
#' block <- rep(1:5, each = 3)
#' block_maxima_ts(data, block = block)
#'
#' ## Data with an incomplete block
#' data <- c(data, 1:2)
#'
#' # Supplying block_length (the extra 2 observations are ignored)
#' block_length <- 3
#' block_maxima_ts(data, block_length = block_length)
#' # Supplying block (with an extra group indicator)
#' block <- c(block, 7, 7)
#' block_maxima_ts(data, block = block)
#' @export
block_maxima_ts <- function(data, block_length, block) {
  # Check that data is a numeric vector
  if (!is.numeric(data)) {
    stop("''data'' must be a numeric vector.")
  }
  block_length_supplied <- !missing(block_length)
  block_supplied <- !missing(block)
  # Check that exactly one of block_length or block has been supplied
  if (block_length_supplied && block_supplied) {
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
         val <- list(maxima = max(temp, na.rm = TRUE),
                     notNA = sum(not_na),
                     whereNA = which(!not_na))
      } else {
         val <- list(maxima = NA,
                     notNA = 0,
                     whereNA = seq_len(length(temp)))
      }
      return(val)
    }
    r <- sapply(seq_len(number_of_blocks), find_maxima_notNA)
    r <- list(maxima = unlist(r[1, ]), notNA = unlist(r[2, ]),
              n = rep_len(block_length, number_of_blocks),
              whereNA = r[3, ])
    names(r[[4]]) <- paste0("block", 1:length(r[[4]]))
  } else {
    if (length(block) != length(data)) {
      stop("''block'' must have the same length as ''data''.")
    }
    FUN <- function(data) {
      not_na <- !is.na(data)
      if (any(not_na)) {
        val <- list(maxima = max(data, na.rm = TRUE), notNA = sum(not_na),
                    n = length(data),
                    whereNA = which(!not_na))
      } else {
        val <- list(maxima = NA, notNA = 0, n = length(data),
                    whereNA = seq_len(length(data)))
      }
      return(val)
    }
    r <- tapply(data, block, FUN = FUN)
    # Extract all the whereNA components (4th in the list)
    whereNA <- sapply(r, `[`, 4)
    names(whereNA) <- paste0("block", 1:length(whereNA))
    # Extract all the other components (1st, 2nd and 3rd in the list)
    others <- sapply(r, `[`, 1:3)
    r <- list(maxima = unlist(others[1, ]), notNA = unlist(others[2, ]),
              n = unlist(others[3, ]), whereNA = whereNA)
    names(r) <- c("maxima", "notNA", "n", "whereNA")
  }
  # For each incomplete block, that is, a block with at least one missing
  # value, apply its missing value pattern to each full block and calculate
  # the resulting pseudo maxima
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
  # Give the returned object a class, so that we can detect block maxima data
  # created by block_maxima_ts()
  class(r) <- c("list", "block_maxima_ts", "evmissing")
  return(r)
}
