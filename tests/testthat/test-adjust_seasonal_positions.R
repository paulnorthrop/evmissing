# Does adjust_seasonal_positions() changes the seasonal indicators correctly?

# (a) Donor block has more seasonal positions than the receiver block

# Receiver block
y <- c(1, 2, 3, 4)
# Donor block x
x <- c(1, 2, 3, 4)

# Case 1
names(x) <- 1:4
newx <- adjust_seasonal_positions(x = x, y = y)
test_that("adjust_seasonal_positions, case 1", {
  expect_equal(names(newx), as.character(c(1, 2, 3, 4)))
})
# Case 2
names(x) <- c(2, 3, 4, 5)
newx <- adjust_seasonal_positions(x = x, y = y)
test_that("adjust_seasonal_positions, case 2", {
  expect_equal(names(newx), as.character(c(2, 3, 4, 1)))
})
# Case 3
names(x) <- c(3, 4, 5, 1)
newx <- adjust_seasonal_positions(x = x, y = y)
test_that("adjust_seasonal_positions, case 3", {
  expect_equal(names(newx), as.character(c(3, 4, 1, 2)))
})
# Case 4
names(x) <- c(4, 5, 1, 2)
newx <- adjust_seasonal_positions(x = x, y = y)
test_that("adjust_seasonal_positions, case 4", {
  expect_equal(names(newx), as.character(c(4, 1, 2, 3)))
})

# (b) Donor block has fewer seasonal positions than the receiver block

# Receiver block
y <- c(1, 2, 3, 4, 5)
# Donor block x
x <- c(1, 2, 3, 4, 5)

# Case 1
names(x) <- c(1, 2, 3, 4, 1)
newx <- adjust_seasonal_positions(x = x, y = y)
test_that("adjust_seasonal_positions, case 1", {
  expect_equal(names(newx), as.character(c(1, 2, 3, 4, 5)))
})
# Case 2
names(x) <- c(2, 3, 4, 1, 2)
newx <- adjust_seasonal_positions(x = x, y = y)
test_that("adjust_seasonal_positions, case 2", {
  expect_equal(names(newx), as.character(c(2, 3, 4, 5, 1)))
})
# Case 3
names(x) <- c(3, 4, 1, 2, 3)
newx <- adjust_seasonal_positions(x = x, y = y)
test_that("adjust_seasonal_positions, case 3", {
  expect_equal(names(newx), as.character(c(3, 4, 5, 1, 2)))
})
# Case 4
names(x) <- c(4, 1, 2, 3, 4)
newx <- adjust_seasonal_positions(x = x, y = y)
test_that("adjust_seasonal_positions, case 4", {
  expect_equal(names(newx), as.character(c(4, 5, 1, 2, 3)))
})
