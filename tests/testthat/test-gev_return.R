# Check that the limits of profile-based CIs for return level increase with
# return period

fit <- gev_mle(PlymouthOzoneMaxima, method = "L-BFGS-B")
rl <- gev_return(fit, m = c(50, 100))

# Profile-based confidence intervals, faster = TRUE
prof <- confint(rl, profile = TRUE, mult = 32, faster = TRUE)
test_that("RL CI lower limits increase with return period, faster = TRUE", {
  testthat::expect_lt(prof[1, 1], prof[2, 1])
})
test_that("RL CI upper limits increase with return period, faster = TRUE", {
  testthat::expect_lt(prof[1, 2], prof[2, 2])
})

# Profile-based confidence intervals, faster = FALSE
prof <- confint(rl, profile = TRUE, mult = 32, faster = FALSE)
test_that("RL CI lower limits increase with return period, faster = FALSE", {
  testthat::expect_lt(prof[1, 1], prof[2, 1])
})
test_that("RL CI upper limits increase with return period, faster = FALSE", {
  testthat::expect_lt(prof[1, 2], prof[2, 2])
})

# Plot method
interval1 <- plot(prof, parm = 1)
interval2 <- plot(prof, parm = 2)
test_that("Symmetric intervals for return level 1", {
  expect_equal(interval1, prof[1, ])
})
test_that("Symmetric intervals for return level 2", {
  expect_equal(interval2, prof[2, ])
})

# Repeat for gev_weighted()

fit <- gev_weighted(PlymouthOzoneMaxima, method = "L-BFGS-B")
rl <- gev_return(fit, m = c(50, 100))

# Profile-based confidence intervals, faster = TRUE
prof <- confint(rl, profile = TRUE, mult = 32, faster = TRUE)
test_that("RL CI lower limits increase with return period, faster = TRUE", {
  testthat::expect_lt(prof[1, 1], prof[2, 1])
})
test_that("RL CI upper limits increase with return period, faster = TRUE", {
  testthat::expect_lt(prof[1, 2], prof[2, 2])
})
# Profile-based confidence intervals, faster = FALSE
prof <- confint(rl, profile = TRUE, mult = 32, faster = FALSE)
test_that("RL CI lower limits increase with return period, faster = FALSE", {
  testthat::expect_lt(prof[1, 1], prof[2, 1])
})
test_that("RL CI upper limits increase with return period, faster = FALSE", {
  testthat::expect_lt(prof[1, 2], prof[2, 2])
})
