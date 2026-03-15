# Check that phat() in evmissing-internal.R works on a simple example

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

# Case 1: both phats are > 1
parameters <- c(mu = 0 , sigma = 1, xi = 0)
phats <- evmissing:::phat(parameters, maxima_notNA)
cond1 <- all(phats == c(1, 1))
cond2 <- all(attr(phats, "unconstrained") > c(1, 1))
test_that("phat(): mu = 0, both phats become 1", {
  testthat::expect_true(cond1 & cond2)
})

# Case 2: One phat is < 1 (but > propn_notNA) and the other > 1
parameters <- c(mu = 1 , sigma = 1, xi = 0)
phats <- evmissing:::phat(parameters, maxima_notNA)
uncon <- attr(phats, "unconstrained")
propn <- attr(phats, "propn_notNA")
cond1 <- phats[1] == 1 & phats[2] < 1 & phats[2] == uncon[2]
cond2 <- uncon[1] >= 1 & uncon[2] < 1 & phats[2] > propn[2]
test_that("phat(): mu = 1, one phat < 1 (but > propn_notNA), the other > 1", {
  testthat::expect_true(cond1 & cond2)
})

# Case 3: One phat is < 1 (but > propn_notNA) and the other is < propn_notNA
parameters <- c(mu = 1.5 , sigma = 1, xi = 0)
phats <- evmissing:::phat(parameters, maxima_notNA)
uncon <- attr(phats, "unconstrained")
propn <- attr(phats, "propn_notNA")
cond1 <- phats[1] < 1 & phats[1] == uncon[1] &
  phats[2] < 1 & phats[2] == propn[2]
cond2 <- uncon[1] < 1 & uncon[1] > propn[1] &
  uncon[2] < 1 & uncon[2] < propn[1]
test_that("phat(): mu = 1.5, both phat < 1, one > propn_notNA, the other not", {
  testthat::expect_true(cond1 & cond2)
})
