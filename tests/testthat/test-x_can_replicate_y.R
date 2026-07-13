# Check the internal function x_can_replicate_y()

x1 <- c(NA, 20, NA, 20, 20, 20)
x2 <- c(NA, NA, NA, 20, NA, 20)
x3 <- c(NA, NA, 20, 20, 20, 20)

# x can replicate y if there are no elements for which x == NA and y != NA

# x1 and x2
test_that("x1 can replicate x2", {
  testthat::expect_true(x_can_replicate_y(x = x1, y = x2))
})
test_that("x2 cannot replicate x1", {
  testthat::expect_false(x_can_replicate_y(x = x2, y = x1))
})

# x1 and x3
test_that("x1 cannot replicate x3", {
  testthat::expect_false(x_can_replicate_y(x = x1, y = x3))
})
test_that("x3 cannot replicate x1", {
  testthat::expect_false(x_can_replicate_y(x = x3, y = x1))
})

# x2 and x3
test_that("x2 cannot replicate x3", {
  testthat::expect_false(x_can_replicate_y(x = x2, y = x3))
})
test_that("x3 can replicate x2", {
  testthat::expect_true(x_can_replicate_y(x = x3, y = x2))
})

# x1, x2 and x3 and themselves
test_that("x1 can replicate x1", {
  testthat::expect_true(x_can_replicate_y(x = x1, y = x1))
})
test_that("x2 can replicate x2", {
  testthat::expect_true(x_can_replicate_y(x = x2, y = x2))
})
test_that("x3 can replicate x3", {
  testthat::expect_true(x_can_replicate_y(x = x3, y = x3))
})
