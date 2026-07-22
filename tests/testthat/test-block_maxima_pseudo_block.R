# Check that block_maxima() works as intended in cases where pseudo = TRUE and
# sliding = FALSE or sliding = TRUE, and
# seasonal = FALSE or seasonal = TRUE

# Uses argument block

# Block sizes: 4, 5, 4, 5
data <- c(1, 2, NA, 3,
          4, 3, NA, NA, 7,
          6, 5, 8, 7,
          1, 2, 3, 4, 5)
block <- c(rep(1, 4), rep(2, 5), rep(3, 4), rep(4, 5))

# 1. !full, sliding, seasonal

col1 <- c(NA, 4, 4, 4, NA, NA, NA, NA, 8, 7, 7, 7, 7, 4, 5)
col2 <- c(4, 4, NA, NA, NA, 7, 7, 7, 7, 6, 5, 3, 3, 5, NA)
c1 <- cbind(col1, col2)
colnames(c1) <- 1:2
rownames(c1) <- 1:nrow(c1)
r1 <- find_pseudo_maxima_block(data = data, block = block,
                               full = FALSE, sliding = TRUE, seasonal = TRUE)
test_that("find_pseudo_block_maxima_block(): block, !full, sliding, seasonal", {
  testthat::expect_equal(r1, c1, ignore_attr = FALSE)
})
# Repeat using years 1990, 1991, 1992, 1993 as block indicators
r1 <- find_pseudo_maxima_block(data = data, block = block + 1990,
                               full = FALSE, sliding = TRUE, seasonal = TRUE)
colnames(c1) <- 1991:1992
test_that("find_pseudo_block_maxima_block(): block, !full, sliding, seasonal", {
  testthat::expect_equal(r1, c1, ignore_attr = FALSE)
})
# Also, check the output from block_maxima
b1 <- block_maxima(data = data, block = block + 1990, pseudo = TRUE,
                   full = FALSE, sliding = TRUE, seasonal = TRUE)
test_that("block_maxima() block, !full, sliding, seasonal", {
  testthat::expect_equal(b1$pseudo_maxima, rm_NA_rows(c1), ignore_attr = FALSE)
})

# 2. !full, sliding, !seasonal

col1 <- c(NA, NA, NA, NA, NA, NA, NA, NA, 8, 7, 8, 8, 7, 4, 5)
col2 <- c(4, NA, NA, NA, NA, NA, NA, NA, 7, 6, 8, 8, 7, 5, NA)
c2 <- cbind(col1, col2)
colnames(c2) <- 1:2
rownames(c2) <- 1:nrow(c2)
r2 <- find_pseudo_maxima_block(data = data, block = block,
                               full = FALSE, sliding = TRUE, seasonal = FALSE)
test_that("find_pseudo_block_maxima_block(): block, !full, sliding, !seasonal", {
  testthat::expect_equal(r2, c2, ignore_attr = FALSE)
})

# 3. full, sliding, seasonal

col1 <- c(NA, NA, NA, NA, NA, NA, NA, NA, 8, 7, 7, 7, 7, 4, 5)
col2 <- c(NA, NA, NA, NA, NA, NA, NA, NA, 7, 6, 5, 3, 3, 5, NA)
c3 <- cbind(col1, col2)
colnames(c3) <- 1:2
rownames(c3) <- 1:nrow(c3)
r3 <- find_pseudo_maxima_block(data = data, block = block,
                               full = TRUE, sliding = TRUE, seasonal = TRUE)
test_that("find_pseudo_block_maxima_block(): block, full, sliding, seasonal", {
  testthat::expect_equal(r3, c3, ignore_attr = FALSE)
})

# 4. full, sliding !seasonal

col1 <- c(NA, NA, NA, NA, NA, NA, NA, NA, 8, 7, 8, 8, 7, 4, 5)
col2 <- c(NA, NA, NA, NA, NA, NA, NA, NA, 7, 6, 8, 8, 7, 5, NA)
c4 <- cbind(col1, col2)
colnames(c4) <- 1:2
rownames(c4) <- 1:nrow(c4)
r4 <- find_pseudo_maxima_block(data = data, block = block,
                               full = TRUE, sliding = TRUE, seasonal = FALSE)
test_that("find_pseudo_block_maxima_block(): block, full, sliding, !seasonal", {
  testthat::expect_equal(r4, c4, ignore_attr = FALSE)
})

# 5. !full, !sliding, !seasonal

col1 <- c(NA, NA, 7, 4)
col2 <- c(4, NA, 6, 5)
c5 <- cbind(col1, col2)
colnames(c5) <- 1:2
rownames(c5) <- 1:nrow(c5)
r5 <- find_pseudo_maxima_block(data = data, block = block,
                               full = FALSE, sliding = FALSE, seasonal = FALSE)
test_that("find_pseudo_block_maxima_block(): block, !full, !sliding, !seasonal", {
  testthat::expect_equal(r5, c5, ignore_attr = FALSE)
})

# 6. full, !sliding, !seasonal

col1 <- c(NA, NA, 7, 4)
col2 <- c(NA, NA, 6, 5)
c6 <- cbind(col1, col2)
colnames(c6) <- 1:2
rownames(c6) <- 1:nrow(c6)
r6 <- find_pseudo_maxima_block(data = data, block = block,
                               full = TRUE, sliding = FALSE, seasonal = FALSE)
test_that("find_pseudo_block_maxima_block(): block, full, !sliding, !seasonal", {
  testthat::expect_equal(r6, c6, ignore_attr = FALSE)
})

