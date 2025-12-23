# The following is based on a dataset simulated by the following code
# set.seed(12345)
#### Exponential data
# res <- sim_study(nsim = 1000, adjust = FALSE, discard = 25)

# For at least one of the simulated datasets, when we discard block maxima
# for blocks with more than 25% missing values there is a problem with the GEV
# fit: optim throws an error. Further investigation suggests that the MLE is not
# away from the parameter boundary and although ismev and evd return results,
# their estimates of xi are close to -1 and issues the are convergence issues.

# This test checks that for these data NA values are returned..
# This is the dataset that causes the problem

maxima <- c(5.037332, 6.472600, 7.631087, 7.349181, 5.459609, 7.488644,
            7.478559, 6.569493, 7.796953, 5.158837, 4.864538, 7.674861,
            4.641970, 7.365161, 7.142155, 6.287263, 6.445659, 7.112482,
            5.389005, 4.925539)
notNA <- c(305, 348, 293, 333, 335, 340, 286, 300, 298, 306,
           328, 347, 296, 348, 349, 345, 354, 356, 331, 278)
n <- rep(365, length(maxima))
init <- c(6.040442, 1.311290, 0.000000)
maxima_notNA <- list(maxima = maxima, notNA = notNA, n = n)
res <- suppressWarnings(gev_mle(maxima_notNA))

val <- rep(NA, 3)
names(val) <- names(coef(res))
test_that("gev_mle: optim throws and error", {
  testthat::expect_equal(coef(res), val)
})
