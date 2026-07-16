# Check that block_maxima() works as intended in cases where pseudo = TRUE and
# sliding = FALSE or sliding = TRUE, and
# seasonal = FALSE or seasonal = TRUE

# Uses argument block_length
# To do: argument block

# A very simple example
data <- c(1:10, 10:1)
# Add some missing values
data[c(3, 8, 9, 19, 20)] <- NA
# Set block_length and block to give the same output (5 blocks of 4 values)
# Call these a_block_length and a_block to avoid over-writing values in setup.R
a_block_length <- 4
a_block <- rep(1:5, each = 4)

maxima <- c(4, 7, 10, 8, 4)
notNA <- c(3, 3, 3, 4, 2)
n <- c(4, 4, 4, 4, 4)
whereNA <- list(3, 4, 1, integer(0), 3:4)
names(whereNA) <- paste0("block", 1:5)
pseudo_maxima <- matrix(c(8, 8, 7, 8), ncol = 4, nrow = 1)
colnames(pseudo_maxima) <- c(1, 2, 3, 5)
rownames(pseudo_maxima) <- 4
full_maxima <- c(8)
partial_maxima <- c(4, 7, 10, 4)
names(partial_maxima) <- c("1", "2", "3", "5")
pseudo <- TRUE
full <- TRUE

# ============================= sliding = FALSE ============================= #

sliding <- FALSE
correctSlidingFALSE <- list(maxima = maxima, notNA = notNA, n = n,
                            whereNA = whereNA, pseudo_maxima = pseudo_maxima,
                            full_maxima = full_maxima,
                            partial_maxima = partial_maxima, pseudo = pseudo,
                            full = full, sliding = sliding, seasonal = FALSE)
resultSlidingFALSE <- block_maxima(data, block_length = a_block_length,
                                   pseudo = pseudo, full = full,
                                   sliding = sliding, seasonal = FALSE)

test_that("block_maxima(): example data 1, pseudo = TRUE, block_length", {
  testthat::expect_equal(resultSlidingFALSE,
                         correctSlidingFALSE, ignore_attr = TRUE)
})

# ============================= sliding = TRUE ============================== #

sliding <- TRUE

#------------------------------ seasonal = TRUE ----------------------------- #

seasonal <- TRUE
row4 <- c(6, 7, 7, 6)
row10 <- c(10, 10, 10, 10)
row11 <- c(9, 10, 10, 8)
row12 <- c(9, 8, 9, 8)
row13 <- c(8, 8, 7, 8)
row14 <- c(7, 7, 7, 7)
row15 <- c(5, 6, 6, 4)
pmSeasonTRUE <- matrix(c(row4, row10, row11, row12, row13, row14, row15),
                       ncol = 4, nrow = 7, byrow = TRUE)

correctSeasonTRUE <- list(maxima = maxima, notNA = notNA, n = n,
                          whereNA = whereNA, pseudo_maxima = pmSeasonTRUE,
                          full_maxima = full_maxima,
                          partial_maxima = partial_maxima, pseudo = pseudo,
                          ful = full, sliding = sliding, seasonal = TRUE)
resultSeasonTRUE <- block_maxima(data, block_length = a_block_length,
                                 pseudo = pseudo, full = full,
                                 sliding = sliding, seasonal = seasonal)

test_that("block_maxima(): example data 1, pseudo = TRUE, block_length", {
  testthat::expect_equal(resultSeasonTRUE,
                         correctSeasonTRUE, ignore_attr = TRUE)
})

#------------------------------ seasonal = FALSE ---------------------------- #

seasonal <- FALSE
row4 <- c(7, 6, 7, 5)
row10 <- c(10, 10, 10, 10)
row11 <- c(10, 10, 9, 10)
row12 <- c(9, 9, 8, 9)
row13 <- c(8, 8, 7, 8)
row14 <- c(7, 7, 6, 7)
row15 <- c(6, 6, 5, 6)
pmSeasonFALSE <- matrix(c(row4, row10, row11, row12, row13, row14, row15),
                        ncol = 4, nrow = 7, byrow = TRUE)

correctSeasonFALSE <- list(maxima = maxima, notNA = notNA, n = n,
                          whereNA = whereNA, pseudo_maxima = pmSeasonFALSE,
                          full_maxima = full_maxima,
                          partial_maxima = partial_maxima, pseudo = pseudo,
                          full = full, sliding = sliding, seasonal = FALSE)
resultSeasonFALSE <- block_maxima(data, block_length = a_block_length,
                                  pseudo = pseudo, full = full, sliding = sliding,
                                  seasonal = FALSE)

test_that("block_maxima(): example data 1, pseudo = TRUE, block_length", {
  testthat::expect_equal(resultSeasonFALSE,
                         correctSeasonFALSE, ignore_attr = TRUE)
})
