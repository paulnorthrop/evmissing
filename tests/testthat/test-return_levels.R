# Check that two different ways to calculate return level estimates agree
# Likewise for the estimated variances of the return level estimators

# Fit a GEV distribution to block maxima from the full data
fit2 <- gev_mle(sdata$data_miss, block_length = sdata$block_length,
                init = "moments")
summary(fit2)

# Extract the MLEs of the GEV parameters
mles <- fit2$par
mu <- mles["mu"]
sigma <- mles["sigma"]
xi <- mles["xi"]
# Set the return period m and number of maxima per year
m <- c(10, 100)
npy <- 1
# Set the annual probability of exceedance based on m and npy
p <- 1 - (1 - 1 / m) ^ (1 / npy)

## Calculation of return level estimates

# Using nieve
rl1 <- nieve::qGEV(p, loc = mu, scale = sigma, shape = xi,
                   lower.tail = FALSE)
# Using my code
rl2 <- mu - sigma * box_cox_vec(x = -log(1 - p), lambda = -xi)

# Using gev_return()

rl3 <- gev_return(fit2, m = c(10, 100))

test_that("Quantiles for return levels: nieve and I agree", {
  testthat::expect_equal(rl1, rl2, ignore_attr = TRUE)
})

test_that("Quantiles for return levels: nieve and gev_return() agree", {
  testthat::expect_equal(rl1, rl3, ignore_attr = TRUE)
})

## Calculation of VC matrix of return level estimators

# Calculate the estimated VC matrix for all return levels in object
yp <- -log(1 - p)

# Using my code
delta <- matrix(0, 3, length(m))
delta[1, ] <- 1
delta[2, ] <- revdbayes::qgev(p, loc = 0, scale = 1, shape = xi,
                              lower.tail = FALSE)
delta[3, ] <- sigma * box_cox_deriv(yp, lambda = -xi)

# Using nieve
temp <- nieve::qGEV(p, loc = mu, scale = sigma, shape = xi, lower.tail = FALSE,
                    deriv = TRUE)
delta_nieve <- t(attr(temp, "gradient"))

test_that("Derivatives of quantiles for return levels: nieve and I agree", {
  testthat::expect_equal(delta, delta_nieve, ignore_attr = TRUE)
})

# Using vcov.return_level()

test_that("vcov for return levels: nieve and vcov.return_level() agree", {
  testthat::expect_equal(t(delta_nieve) %*% vcov(fit2) %*% delta_nieve,
                         vcov(rl3), ignore_attr = TRUE)
})
