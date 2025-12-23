# Check that gev_influence() behaves as expected

# Influence curve for the default mu = 0, sigma = 1, xi = 0 case
z <- seq(from = -3, to = 3, by = 1)
inf1 <- gev_influence(z = z)

# Influence curve for the default mu = 0, sigma = 1, xi = 0 case
z <- seq(from = -3, to = 3, by = 1)
inf2 <- gev_influence(z = z, mu = 2)

# Influence curve for the default mu = 0, sigma = 1, xi = 0 case
z <- seq(from = -3, to = 3, by = 1)
inf3 <- gev_influence(z = z, mu = 2, sigma = 2)

test_that("GEV influence for mu and sigma scale with sigma", {
  testthat::expect_equal(inf2[, 3:4] * 2, inf3[, 3:4])
})

# Check that the GEV influence function for mu is identical to that for
# the 1 / (1 - exp(-1)) return level
z <- seq(from = -3, to = 3, by = 1)
m <- c(1 / (1 - exp(-1)), 50, 100)
infmu <- gev_influence(z = z, mu = 0, sigma = 1, xi = 0.1)
infrl <- gev_influence_rl(z = z, mu = 0, sigma = 1, xi = 0.1, m = m)

test_that("GEV influence fn for mu = that for correspnding return level", {
  testthat::expect_equal(infmu[, 3], infrl[, 3])
})
