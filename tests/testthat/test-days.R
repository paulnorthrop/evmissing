# Check days_in_year() and days_in_month()

test_that("Days in year", {
  testthat::expect_equal(days_in_year(2020:2025),
                         c(366, 365, 365, 365, 366, 365))
})

test_that("Days in month, one year", {
  testthat::expect_equal(days_in_month(2024, 1:12),
                         c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
})

test_that("Days in February, 3 years", {
  testthat::expect_equal(days_in_month(2023:2025, 2),
                         c(28, 29, 28))
})

test_that("Days in month, 3 years", {
  testthat::expect_equal(days_in_month(2023:2025, 1:3),
                         c(31, 28, 31, 31, 29, 31, 31, 28, 31))
})
