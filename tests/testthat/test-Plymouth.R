# Check that that fitting to the Plymouth data gives the correct MLEs and SEs
# in various circumstances, including checks against rounded estimates from
# Table 3

# adjust = TRUE, discard = 0

fitPlymouth1 <- gev_mle(PlymouthOzoneMaxima)
mles <- c(mu = 128.77, sigma = 18.81, xi = 0.00)
ses <- c(mu = 4.40, sigma = 2.63, xi = 0.16)
test_that("adjusted GEV fit to Plymouth Ozone data gives the correct MLEs", {
  testthat::expect_equal(round(coef(fitPlymouth1), 2), mles)
})
test_that("adjusted GEV fit to Plymouth Ozone data gives the correct SEs", {
  testthat::expect_equal(round(sqrt(diag(vcov(fitPlymouth1))), 2), ses)
})

# Check nobs.evmissing()
test_that("nobs.evmissing, no discarging of data", {
  testthat::expect_equal(nobs(fitPlymouth1), nrow(PlymouthOzoneMaxima))
})

# adjust = FALSE, discard = 0

fitPlymouth2 <- gev_mle(PlymouthOzoneMaxima, adjust = FALSE)
mles2 <- c(mu = 126.52, sigma = 25.50, xi = -0.28)
ses2 <- c(mu = 5.53, sigma = 4.00, xi = 0.15)
test_that("unadjusted GEV fit to Plymouth Ozone data gives the correct MLEs", {
  testthat::expect_equal(round(coef(fitPlymouth2), 2), mles2)
})
test_that("unadjusted GEV fit to Plymouth Ozone data gives the correct SEs", {
  testthat::expect_equal(round(sqrt(diag(vcov(fitPlymouth2))), 2), ses2)
})

## Setting discard = 50, say, removes the maxima for years 2001 and 2006

# adjust = TRUE, discard = 50

fitPlymouth3 <- gev_mle(PlymouthOzoneMaxima, adjust = TRUE, discard = 50)
mles3 <- c(mu = 129.55, sigma = 17.65, xi = 0.04)
ses3 <- c(mu = 4.60, sigma = 3.41, xi = 0.26)
test_that("Discard 2001 and 2006, adjust = TRUE, Plymouth Ozone data: MLEs", {
  testthat::expect_equal(round(coef(fitPlymouth3), 2), mles3)
})
test_that("Discard 2001 and 2006, adjust = TRUE, Plymouth Ozone data: SEs", {
  testthat::expect_equal(round(sqrt(diag(vcov(fitPlymouth3))), 2), ses3)
})

# adjust = FALSE, discard = 50

fitPlymouth4 <- gev_mle(PlymouthOzoneMaxima, adjust = FALSE, discard = 50)
mles4 <- c(mu = 128.59, sigma = 17.78, xi = 0.04)
ses4 <- c(mu = 4.55, sigma = 3.61, xi = 0.27)
test_that("Discard 2001 and 2006, adjust = FALSE, Plymouth Ozone data: MLEs", {
  testthat::expect_equal(round(coef(fitPlymouth4), 2), mles4)
})
test_that("Discard 2001 and 2006, adjust = FALSE, Plymouth Ozone data: SEs", {
  testthat::expect_equal(round(sqrt(diag(vcov(fitPlymouth4))), 2), ses4)
})

## Repeat by removing the maxima for years 2001 and 2006 manually

rmYears <- is.element(rownames(PlymouthOzoneMaxima), c(2001, 2006))
newData <- PlymouthOzoneMaxima[!rmYears, ]

# Check nobs.evmissing()
test_that("nobs.evmissing after discarding data", {
  testthat::expect_equal(nobs(fitPlymouth4), nrow(newData))
})

# adjust = TRUE

fitPlymouth5 <- gev_mle(newData, adjust = TRUE)
test_that("Removed 2001 and 2006, adjust = TRUE, Plymouth Ozone data: MLEs", {
  testthat::expect_equal(round(coef(fitPlymouth5), 2), mles3)
})
test_that("Removed 2001 and 2006, adjust = TRUE, Plymouth Ozone data: SEs", {
  testthat::expect_equal(round(sqrt(diag(vcov(fitPlymouth5))), 2), ses3)
})

# adjust = FALSE

fitPlymouth6 <- gev_mle(newData, adjust = FALSE)
test_that("Removed 2001 and 2006, adjust = TRUE, Plymouth Ozone data: MLEs", {
  testthat::expect_equal(round(coef(fitPlymouth6), 2), mles4)
})
test_that("Removed 2001 and 2006, adjust = FALSE, Plymouth Ozone data: SEs", {
  testthat::expect_equal(round(sqrt(diag(vcov(fitPlymouth6))), 2), ses4)
})

## Also check that discard and manual removal of data give the same results

test_that("Plymouth, discard = 50 vs manual removal, adjust = TRUE, MLEs", {
  testthat::expect_equal(coef(fitPlymouth3), coef(fitPlymouth5))
})
test_that("Plymouth, discard = 50 vs manual removal, adjust = TRUE, SEs", {
  testthat::expect_equal(sqrt(diag(vcov(fitPlymouth3))),
                              sqrt(diag(vcov(fitPlymouth5))))
})

test_that("Plymouth, discard = 50 vs manual removal, adjust = FALSE, MLEs", {
  testthat::expect_equal(coef(fitPlymouth4), coef(fitPlymouth6))
})
test_that("Plymouth, discard = 50 vs manual removal, adjust = FALSE, SEs", {
  testthat::expect_equal(sqrt(diag(vcov(fitPlymouth4))),
                         sqrt(diag(vcov(fitPlymouth6))))
})

# Check that the return level estimates and confidence intervals returned
# from plot.evmissing for profile = FALSE and profile = TRUE are as expected

profFALSE <- plot(fitPlymouth1, profile = FALSE, num = 2, m = c(2, 100))
profTRUE <- plot(fitPlymouth1, profile = TRUE, num = 2, m = c(2, 100))

test_that("plot.evmissing, RLs equal for prof = TRUE and prof = FALSE", {
  testthat::expect_equal(profTRUE[, 2], profFALSE[, 2])
})
test_that("plot.evmissing, RL CI upper prof = TRUE > prof = FALSE", {
  testthat::expect_gt(profTRUE[2, 3], profFALSE[2, 3])
})
test_that("plot.evmissing, RL CI lower prof = TRUE > prof = FALSE", {
  testthat::expect_gt(profTRUE[2, 1], profFALSE[2, 1])
})

# Check print.conf_gev()

confPlymouth <- confint(fitPlymouth1, profile = TRUE)

test_that("print.conf_gev() returns its first argument", {
  testthat::expect_equal(print(confPlymouth), confPlymouth)
})

# Check plot.conf_gev()

conf_mu <- plot(confPlymouth)
conf_sigma <- plot(confPlymouth, parm = "sigma")
conf_xi <- plot(confPlymouth, parm = "xi")

test_that("plot.conf_gev() returns the correct interval for mu", {
  testthat::expect_equal(confPlymouth[1, ], conf_mu)
})

test_that("plot.conf_gev() returns the correct interval for sigma", {
  testthat::expect_equal(confPlymouth[2, ], conf_sigma)
})

test_that("plot.conf_gev() returns the correct interval for xi", {
  testthat::expect_equal(confPlymouth[3, ], conf_xi)
})
