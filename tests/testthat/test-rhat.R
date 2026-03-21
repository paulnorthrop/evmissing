# Check that rhat() in evmissing-internal.R works on a simple example

## Simulate example data
set.seed(7032025)
data <- rexp(15)
# Create some missing values
data[c(5, 7:8)] <- NA
# 5 blocks (columns), each with 3 observations
# Supplying block
block <- rep(1:5, each = 3)

# Calculate the block maxima and missing value information
maxima_notNA <- block_maxima_ts(data, block = block)

# Set some default GEV parameters
# Use different values of mu to create estimates of p_i that are inside and
# outside the interval (n_i / n, 1)

# Case 1: both rhats are > 1
parameters <- c(mu = 0 , sigma = 1, xi = 0)
rhats <- evmissing:::rhat(parameters, maxima_notNA)
cond1 <- all(rhats == c(1, 1))
cond2 <- all(attr(rhats, "unconstrained") > c(1, 1))
test_that("rhat(): mu = 0, both rhats become 1", {
  testthat::expect_true(cond1 & cond2)
})

# Case 2: One rhat is < 1 (but > propn_notNA) and the other > 1
parameters <- c(mu = 1 , sigma = 1, xi = 0)
rhats <- evmissing:::rhat(parameters, maxima_notNA)
uncon <- attr(rhats, "unconstrained")
propn <- attr(rhats, "propn_notNA")
cond1 <- rhats[1] == 1 & rhats[2] < 1 & rhats[2] == uncon[2]
cond2 <- uncon[1] >= 1 & uncon[2] < 1 & rhats[2] > propn[2]
test_that("rhat(): mu = 1, one rhat < 1 (but > propn_notNA), the other > 1", {
  testthat::expect_true(cond1 & cond2)
})

# Case 3: One rhat is < 1 (but > propn_notNA) and the other is < propn_notNA
parameters <- c(mu = 1.5 , sigma = 1, xi = 0)
rhats <- evmissing:::rhat(parameters, maxima_notNA)
uncon <- attr(rhats, "unconstrained")
propn <- attr(rhats, "propn_notNA")
cond1 <- rhats[1] < 1 & rhats[1] == uncon[1] &
  rhats[2] < 1 & rhats[2] == propn[2]
cond2 <- uncon[1] < 1 & uncon[1] > propn[1] &
  uncon[2] < 1 & uncon[2] < propn[1]
test_that("rhat(): mu = 1.5, both rhat < 1, one > propn_notNA, the other not", {
  testthat::expect_true(cond1 & cond2)
})
