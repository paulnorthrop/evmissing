#' Block Maxima and Missing Information
#'
#' Extracts disjoint block maxima from a time series of raw data. Can also
#' provide information about the effect of missing value patterns on block
#' maxima via pseudo-maxima created by applying blockwise missing value
#' patterns in partially-observed (partial) disjoint blocks to fully-observed
#' (full) blocks or to full blocks and other suitable partially-observed
#' blocks.
#'
#' @param data A numeric vector containing a time series of raw data.
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
#' @param pseudo A logical scalar. If `pseudo = TRUE` then pseudo-maxima are
#'   calculated as the block maxima obtained by applying the missing value
#'   patterns from partial disjoint blocks to all suitable other blocks, which
#'   we call **donor** blocks. See `full`.
#' @param full If `full = FALSE` then a donor block for a given partial
#'   disjoint block is any block that has non-missing values in all the
#'   comparable positions where the partial block has non-missing values.
#'   The comparable positions depend on whether we respect the seasonality that
#'   may be exhibited within a block. See `seasonal`.
#'   If `full = TRUE` then only full blocks are used as donor blocks:
#'   if `sliding = FALSE` then `data` must contain at least one full disjoint
#'   block and if `sliding = TRUE` then `data` must contain at least one full
#'   sliding block.
#' @param sliding A logical scalar. Only relevant if `pseudo = TRUE`. If
#'   `sliding = TRUE` then pseudo-maxima are calculated for **all** donor blocks
#'   of the relevant length, rather than only a set of disjoint donor blocks.
#' @param seasonal A logical scalar. Only relevant if `pseudo = TRUE` and
#'   `sliding = TRUE`. If `seasonal = TRUE` then the way in which the
#'   pseudo-maxima are calculated respects the seasonality that may be
#'   exhibited over the duration of a block. If, for example, a block covers a
#'   single year, then the missing values applied to a donor block occur at the
#'   same time of year as in the originating partial block. If
#'   `seasonal = FALSE` then the missing values applied have the same positions
#'   in the partial and donor blocks with respect to the start of each block.
#' @details Exactly one of the arguments `block_length` or `block` must be
#'   supplied.
#'
#'   If `block_length` is supplied and `sliding = TRUE` then the pseudo-maxima
#'   are calculated for **all** blocks of length `block_length` present in
#'   `data`, starting with the first block `data[1:block_length]` and sliding
#'   the block repeatedly by one observation until reaching the final block
#'   `data[(length(data) - block_length + 1):length(data)]`.
#'
#'   **Also explain for `block`**
#'
#'   ## Calculation of pseudo-maxima
#'
#'   Consider example data `c(1, 2, NA, 3, 6, 5, 8, 7, 4, 3, NA, NA)` and a
#'   block length of 4. There are 3 disjoint blocks, with block maxima
#'   `c(3, 8, 4)`. Blocks 1 and 3 are partially-observed. We consider the 6
#'   valid combinations of `full`, `sliding` and `seasonal`.
#'
#'   ### Full donor blocks only
#'
#'   * `sliding = FALSE`. Disjoint block 2, `c(6, 5, 8, 7)`, is the only donor
#'     block. Applying the missing value patterns from disjoint blocks 1 and 3
#'     to block 2 leads to respective pseudo-maxima `c(7, 6)`.
#'   * `sliding = TRUE` and `seasonal = FALSE`. There are 4 full sliding donor
#'     blocks. For example, sliding block 4, `c(3, 6, 5, 8)`, leads to
#'     pseudo-maxima `c(8, 6)` for disjoint blocks 1 and 3. The full sets of
#'     pseudo-maxima are `c(8, 7, 8, 8)` for block 1 and `c(6, 6, 8, 8)` for
#'     block 3.
#'   * `sliding = TRUE` and `seasonal = TRUE`. The vector of seasons is
#'     `c(1:4, 1;4, 1:4)`. Sliding block 4, `c(3, 6, 5, 8)` has season vector
#'     `c(4, 1, 2, 3)` leads to pseudo-maxima `c(6, 6)` for disjoint blocks 1
#'     and 3 because the `8` becomes `NA` for block 1 and both `3` and `8`
#'     become `NA` for block 3. The full sets of pseudo-maxima are
#'     `c(6, 6, 6, 7)` for block 1 and `c(6, 6, 5, 4)` for block 3.
#'
#'   ### All suitable donor blocks
#'
#'   * `sliding = FALSE`. Disjoint block 1, `c(1, 2, NA, 3)`, can donate to
#'     disjoint block 3, `c(4, 3, NA, NA)`, leading to a pseudo-maximum of
#'     `2`. The full sets of pseudo-maxima are `c(7)` for block 1 and `c(2, 6)`
#'     for block 3.
#'   * `sliding = TRUE` and `seasonal = FALSE`. In addition to the full sliding
#'     blocks, sliding blocks 1 `c(1, 2, NA, 3)` and 8, `c(7, 4, 3, NA)`, can
#'     donate to disjoint block 3, `c(4, 3, NA, NA)`. The full sets of
#'     pseudo-maxima are `c(8, 7, 8, 7, 8)` for block 1 and
#'     `c(2, 6, 6, 8, 8, 7)` for block 3.
#'   * `sliding = TRUE` and `seasonal = TRUE`. In addition to the full sliding
#'     blocks, sliding block 1 can donate to disjoint block 3 and sliding
#'     blocks 2, 3 and 8 donate to both disjoint blocks 1 and 3. The full sets
#'     of pseudo-maxima are `c(6, 6, 6, 7, 7, 7, 7)` for block 1 and
#'     `c(2, 6, 6, 6, 6, 5, 4, 4)` for block 3.
#'
#'   See [`gev_ts`] for an explanation of how the pseudo-maxima are used.
#'
#' @return A list, with class
#'   `c("list", "block_maxima", "disjoint", "evmissing")`,
#'   if `sliding = FALSE`, or
#'   `c("list", "block_maxima", "sliding", "evmissing")`,
#'   if `sliding = TRUE`,
#'   containing the following numeric vectors:
#'
#'  * `maxima`: the block maxima.
#'  * `notNA`: the numbers of non-missing observations in each block.
#'  * `n`: the maximal block length, that is, the largest number of values that
#'     could have been observed in each block.
#'
#' If `block` is supplied then these vectors are named using the values in
#' `block`. Otherwise, they do not have names.
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
#'    partial block and the rows by the number of the full (disjoint or
#'    sliding) block. If a partial block contains all missing values then its
#'    entry in `pseudo_maxima` is `NA`. If there are no full blocks or no
#'    partial blocks then `pseudo_maxima` is `NA`.
#'  * `full_maxima`: a numeric vector of maxima from full blocks.
#'  * `partial_maxima`: a numeric vector of maxima from partial blocks.
#'
#' The input arguments `pseudo`, `full`, `sliding` and `seasonal` are also
#' included.
#'
#' If a block contains only missing values then its value of `maxima` is `NA`
#' and its value of `notNA` is `0`.
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
#'
#' # Supplying block_length, sliding maxima
#' res1 <- block_maxima(data, block_length = block_length, pseudo = TRUE)
#' res2 <- block_maxima(data, block_length = block_length, pseudo = TRUE,
#'                      sliding = TRUE)
#' # res1 and res2 have the same missing positions, full and partial maxima
#' res1$whereNA
#' res1$full_maxima
#' res1$partial_maxima
#' # res2 has more pseudo maxima than res1
#' res1$pseudo_maxima
#' res2$pseudo_maxima
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
#' @seealso Plot method [`plot.block_maxima`].
#' @export
block_maxima <- function(data, block_length, block, pseudo = FALSE,
                         full = FALSE, sliding = FALSE, seasonal = sliding) {
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
  # Check that the block lengths in block differ by at most 1
  if (block_supplied) {
    block_length_range <- diff(range(table(block)))
    if (block_length_range > 1) {
      stop("The block lengths in ''block'' my differ by at most 1")
    }
  }
  # The seasonal option is only relevant if sliding = TRUE
  if (pseudo && !sliding && seasonal) {
    stop("If sliding = FALSE is then seasonal = TRUE has no relevance")
  }
  # Find block maxima depending on whether block_length or block is supplied
  # If pseudo = TRUE then we also find the positions of missing values
  if (block_length_supplied) {
    r <- maxima_block_length(data = data, block_length = block_length,
                             pseudo = pseudo)
  } else {
    r <- maxima_block(data = data, block = block,
                      pseudo = pseudo)
  }
  # Identify full and incomplete blocks
  nmissing <- r$n - r$notNA
  full_blocks <- which(nmissing == 0)
  incomplete_blocks <- which(nmissing > 0)
  n_full <- length(full_blocks)
  n_incomplete <- length(incomplete_blocks)
  # If full = TRUE and there are no full blocks, or if there are no incomplete
  # blocks then return NA
  if ((full && n_full == 0) || n_incomplete == 0) {
    pseudo_maxima <- NA
  }
  # If pseudo = TRUE then, for each partial block, apply its missing value
  # pattern to each full block and calculate the resulting pseudo maxima
  if (pseudo) {
    if (block_length_supplied) {
#      pseudo_maxima <- pseudo_maxima_block_length(maxima_notNA = r, data = data,
#                                                  block_length = block_length,
#                                                  sliding = sliding,
#                                                  season = seasonal)
      pseudo_maxima <- find_pseudo_maxima_block_length(data = data,
                                                       block_length =
                                                         block_length,
                                      full = full, sliding = sliding,
                                      seasonal = seasonal)
    } else {
#      pseudo_maxima <- pseudo_maxima_block(maxima_notNA = r, data = data,
#                                           block = block)
      pseudo_maxima <- find_pseudo_maxima_block(data = data, block = block,
                                                full = full, sliding = sliding,
                                                seasonal = seasonal)
    }
#    # Name columns and rows according to the positions of block in the raw data
#    #   columns: index of the disjoint receiver block
#    #   rows: index of the (disjoint or sliding) donor block
#    colnames(pseudo_maxima) <- incomplete_blocks
#    rownames(pseudo_maxima) <- seq_len(nrow(pseudo_maxima))
    # Remove rows of all NA from pseudo_maxima
    pseudo_maxima <- rm_NA_rows(pseudo_maxima)
    # How we create the row names depends on whether sliding = TRUE or FALSE

    r <- c(r, list(pseudo_maxima = pseudo_maxima))
    # Create vectors that contain the full maxima and partial maxima
    r$full_maxima <- r$maxima[r$notNA == r$n]
    r$partial_maxima <- r$maxima[r$notNA < r$n]
  }
  # Add the input arguments pseudo, full, sliding and seasonal
  r$pseudo <- pseudo
  r$full <- full
  r$sliding <- sliding
  r$seasonal <- seasonal
  # Give the returned object a class, so that we can detect block maxima data
  # created by block_maxima()
  if (sliding) {
    class(r) <- c("list", "block_maxima", "sliding", "evmissing")
  } else {
    class(r) <- c("list", "block_maxima", "disjoint", "evmissing")
  }
  return(r)
}
