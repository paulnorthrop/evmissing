# Check that gev_bayes() works as expected when there are no missing data

# sdata are simulated in test/testthat/setup.R
# and block_length and block are set there
# Update: now I get these from sdata

# Supply the same block information using block_length and block
n <- 10
set.seed(10)
fit <- gev_bayes(sdata$data_full, block = sdata$block, n = n)
# Repeat with the same seed and with list input from block_maxima()
set.seed(10)
sdata_list <- block_maxima(sdata$data_full,  block = sdata$block)
fit_no_adjust <- gev_bayes(sdata_list, block = sdata$block,
                           adjust = FALSE, n = n)

# adjust = TRUE vs adjust = FALSE
test_that("gev_bayes(): full data, block, adjust = TRUE vs adjust = FALSE", {
  testthat::expect_equal(fit$sim_vals,
                         fit_no_adjust$sim_vals)
})

# Repeat with the same seed and with data frame input
set.seed(10)
fit_no_adjust_df <- gev_bayes(as.data.frame(sdata_list), block = sdata$block,
                              adjust = FALSE, discard = 1, n = n)

# list input vs data frame input (and discard > 0)
test_that("gev_bayes(): full data, list input vs df input", {
  testthat::expect_equal(fit_no_adjust$sim_vals,
                         fit_no_adjust_df$sim_vals)
})
