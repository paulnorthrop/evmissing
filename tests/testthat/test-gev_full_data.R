# Check that gev_mle() and gev_ts() work when there are no missing data

# sdata are simulated in test/testthat/setup.R
# and block_length and block are set there
# Update: now I get these from sdata

# Supply the same block information using block_length and block
fit_block_length <- gev_mle(sdata$data_full, block_length = sdata$block_length)
fit_block <- gev_mle(sdata$data_full, block = sdata$block)
fit_block_no_adjust <- gev_mle(sdata$data_full, block = sdata$block,
                               adjust = FALSE)

# block_length vs block
test_that("gev_mle(): full data, block_length vs block, adjust = TRUE", {
  testthat::expect_equal(fit_block_length, fit_block, ignore_attr = TRUE)
})

# adjust = TRUE vs adjust = FALSE
# Note: object$adjust must be different!
fit_block_no_adjust$adjust <- TRUE
test_that("gev_mle(): full data, block, adjust = TRUE vs adjust = FALSE", {
  testthat::expect_equal(fit_block, fit_block_no_adjust, ignore_attr = TRUE)
})

# gev_ts()

# gev_ts() with pseudo = TRUE vs gev_mle()
fit_block_length_ts <- gev_ts(sdata$data_full,
                              block_length = sdata$block_length,
                              pseudo = TRUE)
test_that("gev_ts() vs gev_mle, coef: full data, block length", {
  testthat::expect_equal(coef(fit_block_length), coef(fit_block_length_ts),
                         ignore_attr = FALSE)
})
test_that("gev_ts() vs gev_mle, vcov: full data, block length", {
  testthat::expect_equal(vcov(fit_block_length), vcov(fit_block_length_ts),
                         ignore_attr = FALSE)
})
