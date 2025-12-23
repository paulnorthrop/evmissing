# Check that BrestSurgeMaxima and BrestSurgeMissing imply the same number of
# missing values

missing1 <- sum(BrestSurgeMissing)
missing2 <- sum(BrestSurgeMaxima$n - BrestSurgeMaxima$notNA)

test_that("Brest datasets are consistent re the numbers of missing values", {
  testthat::expect_equal(missing1, missing2)
})
