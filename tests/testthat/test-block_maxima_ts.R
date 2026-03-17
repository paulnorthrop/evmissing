# Check that block_maxima() works as intended

# block_maxima() is data is a list, such as the sdata created in setup.R

test_that("block_maxima_ts() when data is a list", {
  testthat::expect_error(block_maxima_ts(sdata))
})

# A very simple example
data <- c(1:10, 10:1)
# Add some missing values
data[c(3, 8, 9, 19, 20)] <- NA
# Set block_length and block to give the same output (5 blocks of 4 values)
# Call these a_block_length and a_block to avoid over-writing values in setup.R
a_block_length <- 4
a_block <- rep(1:5, each = 4)

# block_maxima_ts() errors unless exactly one of block_length or block is supplied
test_that("block_maxima_ts() errors when neither block_length or block are supplied", {
  testthat::expect_error(block_maxima_ts(data))
})
test_that("block_maxima() errors when both block_length and block are supplied", {
  testthat::expect_error(block_maxima_ts(data, a_block_length, a_block))
})

# block_maxima() errors when block does no have the same length as data
test_that("block_maxima() errors when length(block) != length(data)", {
  testthat::expect_error(block_maxima_ts(data, block = a_block[1:5]))
})

maxima <- c(4, 7, 10, 8, 4)
notNA <- c(3, 3, 3, 4, 2)
n <- c(4, 4, 4, 4, 4)
whereNA <- list(3, 4, 1, integer(0), 3:4)
names(whereNA) <- paste0("block", 1:5)
pseudo_maxima <- matrix(c(8, 8, 7, 8), ncol = 4, nrow = 1)
colnames(pseudo_maxima) <- c(1, 2, 3, 5)
rownames(pseudo_maxima) <- 4
full_maxima <- c("4" = 8)
partial_maxima <- c(4, 7, 10, 4)
names(partial_maxima) <- c("1", "2", "3", "5")
results <- list(maxima = maxima, notNA = notNA, n = n,
                whereNA = whereNA, pseudo_maxima = pseudo_maxima,
                full_maxima = full_maxima, partial_maxima = partial_maxima)

test_that("block_maxima_ts(): example data 1, block gives correct result", {
  testthat::expect_equal(block_maxima_ts(data, block = a_block),
                         results, ignore_attr = TRUE)
})
test_that("block_maxima_ts(): example data 1, block_length gives correct result", {
  testthat::expect_equal(block_maxima_ts(data, block_length = a_block_length),
                         results, ignore_attr = TRUE)
})

# Start again with the same example data but modify them

# A simple example
data <- c(1:10, 10:1)
# Add some missing values
data[c(3, 8, 9, 19, 20)] <- NA
# Create data with 2 extra blocks:
#   One with all (4) NA
#   Another with an incomplete block of length 3
#data <- c(data, rep(NA, 4), c(1, 2, NA))
data <- c(data, rep(NA, 4), c(NA, 1, 2))
# Remove data[3] and set block so that the first block full and length 3
data <- data[-3]
# If block is supplied (correctly) then the incomplete block is not ignored
a_new_block <- c(rep(1, 3),
                 rep(2, 4), rep(3, 4), rep(4, 4), rep(5, 4), rep(6, 4),
                 rep(7, 3))

# (a) the final block is incomplete and shorter (length 3) than the
# lengths (4) of many other blocks, including the full block (block 4)
# (b) the first block is full and shorter (length 3) than the
# lengths (4) of many of the other blocks, including some incomplete blocks

# If block_length is supplied then the final block should be ignored
maxima <- c(5, 7, 10, 7, 3, NA)
notNA <- c(4, 2, 4, 4, 1, 0)
n <- rep(4, 6)
whereNA1 <- list(integer(0), 3:4, integer(0), integer(0), 2:4, 1:4)
names(whereNA1) <- paste0("block", 1:6)
pseudo_maxima <- matrix(c(2, 10, 7, 1, 10, 7, NA, NA, NA), ncol = 3, nrow = 3)
colnames(pseudo_maxima) <- c(2, 5, 6)
rownames(pseudo_maxima) <- c(1, 3, 4)
full_maxima <- c(5, 10, 7)
partial_maxima <- c(7, 3, NA)
block_length_results <- list(maxima = maxima, notNA = notNA, n = n,
                             whereNA = whereNA1, pseudo_maxima = pseudo_maxima,
                             full_maxima = full_maxima,
                             partial_maxima = partial_maxima)

maxima <- c(4, 7, 10, 8, 4, NA, 2)
notNA <- c(3, 3, 3, 4, 2, 0, 2)
n <- c(3, 4, 4, 4, 4, 4, 3)
whereNA2 <- list(integer(0), 4, 1, integer(0), 3:4, 1:4, 1)
names(whereNA2) <- paste0("block", 1:7)
pseudo_maxima <- matrix(c(4, 8, 4, 7, 2, 8, NA, NA, 4, 7), ncol = 5, nrow = 2)
colnames(pseudo_maxima) <- c(2, 3, 5, 6, 7)
rownames(pseudo_maxima) <- c(1, 4)
full_maxima <- c(4, 8)
names(full_maxima) <- c(1, 4)
partial_maxima <- c(7, 10, 4, NA, 2)
names(partial_maxima) <- c(2, 3, 5, 6, 7)
block_results <- list(maxima = maxima, notNA = notNA, n = n,
                      whereNA = whereNA2, pseudo_maxima = pseudo_maxima,
                      full_maxima = full_maxima,
                      partial_maxima = partial_maxima)

test_that("block_maxima(): example data 2, block gives correct result", {
  testthat::expect_equal(block_maxima_ts(data, block = a_new_block),
                         block_results, ignore_attr = TRUE)
})
test_that("block_maxima(): example data 2, block_length gives correct result", {
  testthat::expect_equal(block_maxima_ts(data, block_length = a_block_length),
                         block_length_results, ignore_attr = TRUE)
})

# Simulate some example data
set.seed(7032025)
data <- stats::rexp(15)
# Set block_length and block to give the same output (5 blocks of 3 values)
# Call these a_block_length and a_block to avoid over-writing values in setup.R
a_block_length <- 3
a_block <- rep(1:5, each = 3)

# block_maxima_ts() gives the same output for equivalent block_length and block

test_that("block_maxima(): simulated data, block and block_length agree", {
  testthat::expect_equal(block_maxima_ts(data, block = a_block),
                         block_maxima_ts(data, block_length = a_block_length),
                         ignore_attr = TRUE)
})
