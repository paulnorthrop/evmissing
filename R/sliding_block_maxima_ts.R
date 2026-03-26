#' Sliding block maxima for a Time Series
#'
#' Extracts sliding block maxima and missing value information for each block.
#' Works like [`block_maxima_ts`] but returns information for sliding
#' (overlapping) blocks rather than disjoint blocks and defines blocks using
#' only block length, that is, there is no `block` argument.
#'
#' @param data A numeric vector containing a time series of raw data.
#' @param block_length A numeric scalar. Used calculate the maxima of sliding
#'   blocks of `block_length` contiguous values in the vector `data`.
#' @details Add details.
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
#'
#' ## Data with an incomplete block
#' data <- c(data, 1:2)
#'
#' # Supplying block_length (the extra 2 observations are ignored)
#' block_length <- 3
#' sliding_block_maxima_ts(data, block_length = block_length)
#' @export
sliding_block_maxima_ts <- function(data, block_length) {
  # Check that data is a numeric vector
  if (!is.numeric(data)) {
    stop("''data'' must be a numeric vector.")
  }
  # Find block maxima
  number_of_blocks <- length(data) - block_length + 1
  find_maxima_notNA <- function(i) {
    temp <- data[i:(i + block_length - 1)]
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
  # For each incomplete block, that is, a block with at least one missing
  # value, apply its missing value pattern to each full block and calculate
  # the resulting pseudo maxima.
  # Pass sliding = TRUE to pseudo_maxima_block_length() to ensure that sliding
  # blocks are identified correctly
  pseudo_maxima <- pseudo_maxima_block_length(maxima_notNA = r, data = data,
                                              block_length = block_length,
                                              sliding = TRUE)
  r <- c(r, list(pseudo_maxima = pseudo_maxima))
  # Create vectors that contain the full maxima and partial maxima
  r$full_maxima <- r$maxima[as.numeric(rownames(r$pseudo_maxima))]
  r$partial_maxima <- r$maxima[as.numeric(colnames(r$pseudo_maxima))]
  # Give the returned object a class, so that we can detect block maxima data
  # created by block_maxima()
  class(r) <- c("list", "block_maxima_ts", "evmissing")
  return(r)
}
