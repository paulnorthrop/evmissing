# Check that the default setting of sim_data() produces exactly the same
# simulated data as Emma's original code

# Note: sdata are simulated in test/testthat/setup.R using the seed 17032025
# and block_length, block and blocks are set there
# Update: now I get these from sdata

the_seed <- 17032025
set.seed(the_seed)
data <- rexp(sdata$block_length * sdata$blocks)
data_original <- data
## Artificially remove some observations at random
missing <- NULL
for(miss_iter in 1:sdata$blocks){
  prop_miss <- runif(1, 0, 0.5)
  num_miss  <- ceiling(prop_miss * sdata$block_length)
  year_miss <- sample(sdata$block_length, num_miss, replace = FALSE)
  missing <- c(missing, (year_miss + sdata$block_length * (miss_iter - 1)))
}
data[missing] <- NA

# Raw data, no missings
test_that("sim_data(): no missing values", {
  testthat::expect_equal(sdata$data_full, data_original)
})
# With missings
test_that("sim_data(): with missing values", {
  testthat::expect_equal(sdata$data_miss, data)
})

# Check mcar2

blocks <- 3
block_length <- 4
pmiss <- 0.25
sdata2 <- sim_data(blocks = blocks, block_length = block_length,
                   missing_fn = mcar2, missing_args = list(pmiss = pmiss))
nmiss <- sum(is.na(sdata2$data_miss))

test_that("sim_data(): mcar2()", {
  testthat::expect_equal(nmiss, blocks * block_length * pmiss)
})

# Tests for simulated (exponential) data for which the point estimate of
# the GEV shape parameter xi is negative

fit_sim <- gev_mle(sdata$data_miss, block_length = sdata$block_length)
conf_sim <- confint(fit_sim, profile = TRUE, faster = TRUE, mult = 16,
                    epsilon = 0)
conf_xi <- plot(conf_sim, parm = "xi", type = "b")

test_that("sim data: plot.conf_gev() returns the correct interval for xi", {
  testthat::expect_equal(conf_sim[3, ], conf_xi)
})
