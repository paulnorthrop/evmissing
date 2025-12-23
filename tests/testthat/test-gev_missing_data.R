# Check that gev_mle() works as expected when there are missing data

# sdata are simulated in test/testthat/setup.R
# and block_length and block are set there
# Update: now I get these from sdata

# Make adjustment for the numbers of non-missing values per block
adj_block_length <- gev_mle(sdata$data_miss,
                            block_length = sdata$block_length)
adj_block <- gev_mle(sdata$data_miss,
                     block = sdata$block)

# Do not make the adjustment
no_adj_block_length <- gev_mle(sdata$data_miss,
                               block_length = sdata$block_length,
                               adjust = FALSE)
no_adj_block <- gev_mle(sdata$data_miss,
                        block = sdata$block, adjust = FALSE)
# block_length vs block
test_that("gev_mle(): missing data, block_length vs block, adjust = TRUE", {
  testthat::expect_equal(adj_block_length, adj_block, ignore_attr = TRUE)
})
test_that("gev_mle(): missing data, block_length vs block, adjust = FALSE", {
  testthat::expect_equal(no_adj_block_length, no_adj_block, ignore_attr = TRUE)
})

# Check that we get the same fit based on block maxima from block_maxima()
bm <- block_maxima(sdata$data_miss, block_length = sdata$block_length)
bm_fit <- gev_mle(bm)
test_that("gev_mle(): missing data, block_length vs block, adjust = FALSE", {
  testthat::expect_equal(bm_fit, adj_block, ignore_attr = TRUE)
})
