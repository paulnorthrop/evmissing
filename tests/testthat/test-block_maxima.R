# Check that block_maxima() works as intended

# block_maxima() is data is a list, such as the sdata created in setup.R

test_that("block_maxima() when data is a list", {
  testthat::expect_error(block_maxima(sdata))
})

# A very simple example
data <- c(1:10, 10:1)
# Add some missing values
data[c(3, 8, 9, 19, 20)] <- NA
# Set block_length and block to give the same output (5 blocks of 4 values)
# Call these a_block_length and a_block to avoid over-writing values in setup.R
a_block_length <- 4
a_block <- rep(1:5, each = 4)

# block_maxima() errors unless exactly one of block_length or block is supplied
test_that("block_maxima() errors when neither block_length or block are supplied", {
  testthat::expect_error(block_maxima(data))
})
test_that("block_maxima() errors when both block_length and block are supplied", {
  testthat::expect_error(block_maxima(data, a_block_length, a_block))
})

# block_maxima() errors when block does no have the same length as data
test_that("block_maxima() errors when length(block) != length(data)", {
  testthat::expect_error(block_maxima(data, block = a_block[1:5]))
})


maxima <- c(4, 7, 10, 8, 4)
notNA <- c(3, 3, 3, 4, 2)
n <- c(4, 4, 4, 4, 4)
results <- list(maxima = maxima, notNA = notNA, n = n)

test_that("block_maxima(): example data 1, block gives correct result", {
  testthat::expect_equal(block_maxima(data, block = a_block),
                         results, ignore_attr = TRUE)
})

test_that("block_maxima(): example data 1, block_length gives correct result", {
  testthat::expect_equal(block_maxima(data, block_length = a_block_length),
                         results, ignore_attr = TRUE)
})

# Create data with 2 extra blocks:
#   One with all (4) NA
#   Another with an incomplete block of length 3
data <- c(data, rep(NA, 4), c(1, 2, NA))
# If block_length is supplied then the incomplete block should be ignored
block_length_results <- list(maxima = c(maxima, NA), notNA = c(notNA, 0),
                             n = c(n, 4))
# If block is supplied (correctly) then the incomplete block is not ignored
a_new_block <- c(a_block, rep(6, 4), rep(7, 3))
block_results <- list(maxima = c(maxima, NA, 2), notNA = c(notNA, 0, 2),
                      n = c(n, 4, 3))

test_that("block_maxima(): example data 2, block gives correct result", {
  testthat::expect_equal(block_maxima(data, block = a_new_block),
                         block_results, ignore_attr = TRUE)
})

test_that("block_maxima(): example data 2, block_length gives correct result", {
  testthat::expect_equal(block_maxima(data, block_length = a_block_length),
                         block_length_results, ignore_attr = TRUE)
})

# Simulate some example data
set.seed(7032025)
data <- stats::rexp(15)
# Set block_length and block to give the same output (5 blocks of 3 values)
# Call these a_block_length and a_block to avoid over-writing values in setup.R
a_block_length <- 3
a_block <- rep(1:5, each = 3)

# block_maxima() gives the same output for equivalent block_length and block

test_that("block_maxima(): simulated data, block and block_length agree", {
  testthat::expect_equal(block_maxima(data, block = a_block),
                         block_maxima(data, block_length = a_block_length),
                         ignore_attr = TRUE)
})
