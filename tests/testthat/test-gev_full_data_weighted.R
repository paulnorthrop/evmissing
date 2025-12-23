# Check that gev_mle() and gev_weighted() give the same results when there are
# no missing data

# sdata are simulated in test/testthat/setup.R
# and block_length and block are set there
# Update: now I get these from sdata

# Fit a GEV distribution to block maxima from the full data
fit0 <- gev_mle(sdata$data_full, block_length = sdata$block_length)

# Fit to the full data using weighting scheme 1
fit1 <- gev_weighted(sdata$data_full, scheme = 1,
                     block_length = sdata$block_length)

# Fit to the full data using weighting scheme 2
fit2 <- gev_weighted(sdata$data_full, scheme = 2,
                     block_length = sdata$block_length)

# coef
test_that("gev_mle() vs gev_weighted scheme 1: coef", {
  testthat::expect_equal(coef(fit0), coef(fit1))
})
test_that("gev_mle() vs gev_weighted scheme 2: coef", {
  testthat::expect_equal(coef(fit0), coef(fit2))
})

# symmetric confint
conf0 <- confint(fit0)
test_that("gev_mle() vs gev_weighted scheme 1: symmetric confint", {
  testthat::expect_equal(conf0, confint(fit1))
})
test_that("gev_mle() vs gev_weighted scheme 2: symmetric confint", {
  testthat::expect_equal(conf0, confint(fit2))
})

# profile confint
conf0 <- confint(fit0, profile = TRUE)
test_that("gev_mle() vs gev_weighted scheme 1: symmetric confint", {
  testthat::expect_equal(conf0, confint(fit1, profile = TRUE))
})
test_that("gev_mle() vs gev_weighted scheme 2: symmetric confint", {
  testthat::expect_equal(conf0, confint(fit2, profile = TRUE))
})
