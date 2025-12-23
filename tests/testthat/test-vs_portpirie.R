# Check GEV inferences vs the portpirie annual maxima in ismev and Coles (2001)

# The Portpirie sea levels
p <- c(4.03, 3.83, 3.65, 3.88, 4.01, 4.08, 4.18, 3.80, 4.36, 3.96, 3.98, 4.69,
       3.85, 3.96, 3.85, 3.93, 3.75, 3.63, 3.57, 4.25, 3.97, 4.05, 4.24, 4.22,
       3.73, 4.37, 4.06, 3.71, 3.96, 4.06, 4.55, 3.79, 3.89, 4.11, 3.85, 3.86,
       3.86, 4.21, 4.01, 4.11, 4.24, 3.96, 4.21, 3.74, 3.85, 3.88, 3.66, 4.11,
       3.71, 4.18, 3.90, 3.78, 3.91, 3.72, 4.00, 3.66, 3.62, 4.33, 4.55, 3.75,
       4.08, 3.90, 3.88, 3.94, 4.33)

# Fit a GEV model
x <- gev_mle(p, block_length = 1)
# From
# ifit <- ismev::gev.fit(p)
# ifit$mle
ismev_mle <- c(3.87474692, 0.19804120, -0.05008773)
names(ismev_mle) <- names(coef(x))

# MLEs
test_that("portpirie: MLEs correct", {
  testthat::expect_equal(coef(x), ismev_mle, tolerance = 1e-4)
})
# Number of observations
# length(ifit$data)
test_that("portpirie: number of observations correct", {
  testthat::expect_equal(nobs(x), 65, tolerance = 1e-4)
})
# Log-likelihood
# -ifit$nllh
test_that("portpirie: log-likelihood correct", {
  testthat::expect_equal(logLik(x), 4.339058, tolerance = 1e-4,
                         ignore_attr = TRUE)
})

# VC matrices
# From page 59 of Coles (2001)
ismev_vc <- matrix(c(0.000780, 0.000197, -0.00107,
                     0.000197, 0.000410, -0.000778,
                     -0.00107, -0.000778, 0.00965), 3, 3)
dimnames(ismev_vc) <- dimnames(vcov(x))
test_that("portpirie: VC matrix correct", {
  testthat::expect_equal(vcov(x), ismev_vc, tolerance = 1e-3)
})

# 10-year and 100-year return levels
rl <- gev_return(x, m = c(10, 100))
# From page 60 of Coles (2001)
ismev_rl <- c(4.30, 4.69)
names(ismev_rl) <- names(rl)
test_that("portpirie: return level estimates correct", {
  testthat::expect_equal(rl, ismev_rl, tolerance = 1e-3, ignore_attr = TRUE)
})

# Variance of the 10-year return level variance
# From page 60 of Coles (2001)
rl10 <- gev_return(x, m = 10)
rl_variance <- vcov(rl10)[1, 1]
ismev_rl_variance <- 0.00303
dimnames(ismev_rl_variance) <- dimnames(rl_variance)
test_that("portpirie: 10-year return level variance estimate correct", {
  testthat::expect_equal(rl_variance, ismev_rl_variance, tolerance = 1e-2,
                         ignore_attr = TRUE)
})

# Calculating 10-year and 100-year return level estimates together and
# separately gives the same results
rl100 <- gev_return(x, m = 100)
rl_together <- coef(rl)
rl_separately <- c(coef(rl10), coef(rl100))
test_that("portpirie: return levels together vs separately", {
  testthat::expect_equal(rl_together, rl_separately)
})

# Similarly for the estimates VC matrices
vc_together <- diag(vcov(rl))
vc_separately <- c(vcov(rl10), vcov(rl100))
names(vc_separately) <- names(vc_together)
test_that("portpirie: variances of return levels together vs separately", {
  testthat::expect_equal(vc_together, vc_separately)
})
