# Check that gev_profile_init() works as expected

gev_profile_init <- getFromNamespace("gev_profile_init", "evmissing")

constraint_fn <- function(pars, data) {
  return(1 + pars[3] * (data - pars[1]) / pars[2])
}

# Exponential raw data
# Make adjustment for the numbers of non-missing values per block
fit <- gev_mle(sdata$data_miss, block_length = sdata$block_length)
cis <- confint(fit)
pars <- cbind(cis[, 1], coef(fit), cis[, 2])
colnames(pars) <- c("lower", "mle", "upper")

# Calculate initial estimates for all cases
data <- fit$maxima
mu_given <- mapply(FUN = gev_profile_init, mu = pars[1, ],
                   MoreArgs = list(data = data))
sigma_given <- mapply(FUN = gev_profile_init, sigma = pars[2, ],
                      MoreArgs = list(data = data))
xi_given <- mapply(FUN = gev_profile_init, xi = pars[3, ],
                   MoreArgs = list(data = data))

# Check that all values of sigma are > 0
test_that("gev_profile_init. mu_given: sigmas > 0", {
  testthat::expect_gt(min(mu_given["sigma", ]), 0)
})
test_that("gev_profile_init. sigma_given: sigmas > 0", {
  testthat::expect_gt(min(sigma_given["sigma", ]), 0)
})
test_that("gev_profile_init. xi_given: sigmas > 0", {
  testthat::expect_gt(min(xi_given["sigma", ]), 0)
})

constraints_mu <- apply(mu_given, 2, constraint_fn, data = data)
constraints_sigma <- apply(sigma_given, 2, constraint_fn, data = data)
constraints_xi <- apply(xi_given, 2, constraint_fn, data = data)

# Check that all observations are within bounds
test_that("gev_profile_init. mu_given: constraints OK", {
  testthat::expect_gt(min(constraints_mu), 0)
})
# Check that all observations are within bounds
test_that("gev_profile_init. sigma_given: constraints OK", {
  testthat::expect_gt(min(constraints_sigma), 0)
})
# Check that all observations are within bounds
test_that("gev_profile_init. xi_given: constraints OK", {
  testthat::expect_gt(min(constraints_xi), 0)
})

# Repeat for t(2) raw data
# Make adjustment for the numbers of non-missing values per block
fit <- gev_mle(sdata_t$data_miss, block_length = sdata_t$block_length)
cis <- confint(fit)
pars <- cbind(cis[, 1], coef(fit), cis[, 2])
colnames(pars) <- c("lower", "mle", "upper")

# Calculate initial estimates for all cases
data <- fit$maxima
mu_given <- mapply(FUN = gev_profile_init, mu = pars[1, ],
                   MoreArgs = list(data = data))
sigma_given <- mapply(FUN = gev_profile_init, sigma = pars[2, ],
                      MoreArgs = list(data = data))
xi_given <- mapply(FUN = gev_profile_init, xi = pars[3, ],
                   MoreArgs = list(data = data))

# Check that all values of sigma are > 0
test_that("gev_profile_init. t-distn, mu_given: sigmas > 0", {
  testthat::expect_gt(min(mu_given["sigma", ]), 0)
})
test_that("gev_profile_init. t-distn, sigma_given: sigmas > 0", {
  testthat::expect_gt(min(sigma_given["sigma", ]), 0)
})
test_that("gev_profile_init. t-distn, xi_given: sigmas > 0", {
  testthat::expect_gt(min(xi_given["sigma", ]), 0)
})

constraints_mu <- apply(mu_given, 2, constraint_fn, data = data)
constraints_sigma <- apply(sigma_given, 2, constraint_fn, data = data)
constraints_xi <- apply(xi_given, 2, constraint_fn, data = data)

# Check that all observations are within bounds
test_that("gev_profile_init. t-distn, mu_given: constraints OK", {
  testthat::expect_gt(min(constraints_mu), 0)
})
# Check that all observations are within bounds
test_that("gev_profile_init. t-distn, sigma_given: constraints OK", {
  testthat::expect_gt(min(constraints_sigma), 0)
})
# Check that all observations are within bounds
test_that("gev_profile_init. t-distn, xi_given: constraints OK", {
  testthat::expect_gt(min(constraints_xi), 0)
})
