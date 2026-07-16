# Check that find_pseudo_maxima_block_length() works as intended
# Also check against block_maxima()$pseudo_maxima

# A very simple example
data <- c(1:10, 10:1)
# Add some missing values
data[c(3, 8, 9, 19, 20)] <- NA
# Set block_length and block to give the same output (5 blocks of 4 values)
# Call these a_block_length and a_block to avoid over-writing values in setup.R
a_block_length <- 4

# 1. full = FALSE, sliding = TRUE, seasonal = FALSE

correct1 <- matrix(NA, nrow = 17, ncol = 4)
correct1[1, ] <- c(NA, NA, NA, 2)
correct1[2, ] <- c(NA, NA, NA, NA)
correct1[3, ] <- c(NA, NA, 6, NA)
correct1[4, ] <- c(7, 6, 7, 5)
correct1[5, ] <- c(NA, NA, NA, 6)
correct1[6, ] <- c(NA, NA, NA, 7)
correct1[7, ] <- c(NA, NA, NA, NA)
correct1[8, ] <- c(NA, NA, NA, NA)
correct1[9, ] <- c(NA, NA, NA, NA)
correct1[10, ] <- c(10, 10, 10, 10)
correct1[11, ] <- c(10, 10, 9, 10)
correct1[12, ] <- c(9, 9, 8, 9)
correct1[13, ] <- c(8, 8, 7, 8)
correct1[14, ] <- c(7, 7, 6, 7)
correct1[15, ] <- c(6, 6, 5, 6)
correct1[16, ] <- c(NA, 5, NA, 5)
correct1[17, ] <- c(NA, NA, NA, NA)

results1 <- find_pseudo_maxima_block_length(data,
                                            block_length = a_block_length,
                                            full = FALSE, sliding = TRUE,
                                            seasonal = FALSE)
test_that(
  "find_pseudo_maxima_block_length(): not full, sliding, not seasonal", {
  testthat::expect_equal(results1, correct1, ignore_attr = TRUE)
})
bm1 <- block_maxima(data, block_length = a_block_length, pseudo = TRUE,
                    full = FALSE, sliding = TRUE, seasonal = FALSE)
test_that(
  "block_maxima(): not full, sliding, not seasonal", {
    testthat::expect_equal(bm1$pseudo_maxima, rm_NA_rows(results1),
                           ignore_attr = TRUE)
  })

# 2. full = FALSE, sliding = TRUE, seasonal = TRUE

correct2 <- matrix(NA, nrow = 17, ncol = 4)
correct2[1, ] <- c(NA, NA, NA, 2)
correct2[2, ] <- c(5, NA, NA, 5)
correct2[3, ] <- c(6, NA, NA, 6)
correct2[4, ] <- c(6, 7, 7, 6)
correct2[5, ] <- c(NA, NA, NA, 6)
correct2[6, ] <- c(NA, NA, NA, NA)
correct2[7, ] <- c(NA, NA, NA, NA)
correct2[8, ] <- c(NA, NA, NA, NA)
correct2[9, ] <- c(NA, NA, NA, NA)
correct2[10, ] <- c(10, 10, 10, 10)
correct2[11, ] <- c(9, 10, 10, 8)
correct2[12, ] <- c(9, 8, 9, 8)
correct2[13, ] <- c(8, 8, 7, 8)
correct2[14, ] <- c(7, 7, 7, 7)
correct2[15, ] <- c(5, 6, 6, 4)
correct2[16, ] <- c(5, NA, NA, 4)
correct2[17, ] <- c(NA, NA, NA, NA)

results2 <- find_pseudo_maxima_block_length(data,
                                            block_length = a_block_length,
                                            full = FALSE, sliding = TRUE,
                                            seasonal = TRUE)
test_that(
  "find_pseudo_maxima_block_length(): not full, sliding, seasonal", {
  testthat::expect_equal(results2, correct2, ignore_attr = TRUE)
})
bm2 <- block_maxima(data, block_length = a_block_length, pseudo = TRUE,
                    full = FALSE, sliding = TRUE, seasonal = TRUE)
test_that(
  "block_maxima(): not full, sliding, seasonal", {
    testthat::expect_equal(bm2$pseudo_maxima, rm_NA_rows(results2),
                           ignore_attr = TRUE)
  })

# 3. full = FALSE, sliding = FALSE (seasonal = FALSE)

correct3a <- correct1[c(1, 5, 9, 13, 17), ]
correct3b <- correct2[c(1, 5, 9, 13, 17), ]
results3 <- find_pseudo_maxima_block_length(data,
                                            block_length = a_block_length,
                                            full = FALSE, sliding = FALSE,
                                            seasonal = FALSE)
test_that(
  "find_pseudo_maxima_block_length(): not full, not sliding, check 1", {
  testthat::expect_equal(results3, correct3a, ignore_attr = TRUE)
})
test_that(
  "find_pseudo_maxima_block_length(): not full, not sliding, check 2", {
  testthat::expect_equal(results3, correct3b, ignore_attr = TRUE)
})

# 4. full = TRUE, sliding = TRUE, seasonal = FALSE

results4 <- find_pseudo_maxima_block_length(data,
                                            block_length = a_block_length,
                                            full = TRUE, sliding = TRUE,
                                            seasonal = FALSE)
correct4 <- correct1
correct4[c(1:3, 5:9, 16:17), ] <- NA
test_that(
  "find_pseudo_maxima_block_length(): full, sliding, not seasonal", {
  testthat::expect_equal(results4, correct4, ignore_attr = TRUE)
})

# 5. full = TRUE, sliding = TRUE, seasonal = TRUE

results5 <- find_pseudo_maxima_block_length(data,
                                            block_length = a_block_length,
                                            full = TRUE, sliding = TRUE,
                                            seasonal = TRUE)
correct5 <- correct2
correct5[c(1:3, 5:9, 16:17), ] <- NA
test_that(
  "find_pseudo_maxima_block_length(): full, sliding, seasonal", {
  testthat::expect_equal(results5, correct5, ignore_attr = TRUE)
})

# 6. full = FALSE, sliding = FALSE (seasonal = FALSE)

results6 <- find_pseudo_maxima_block_length(data,
                                            block_length = a_block_length,
                                            full = TRUE, sliding = FALSE,
                                            seasonal = FALSE)
correct6 <- correct3a
correct6[c(1:3, 5), ] <- NA
test_that(
  "find_pseudo_maxima_block_length(): full, not sliding", {
  testthat::expect_equal(results6, correct6, ignore_attr = TRUE)
})
