# Check the output from block_maxima() for the LaPlagne data

# Precipitation block maxima (wlog use LaPlagnePrecipMaximaSeasonal)

m <- tapply(LaPlagne$rr24, INDEX = list(LaPlagne$winter), FUN = max, na.rm = TRUE)
p <- tapply(LaPlagne$rr24, INDEX = list(LaPlagne$winter),
            FUN = function(x) sum(!is.na(x)) / length(x))
full_maxima <- m[p == 1]
partial_maxima <- m[p != 1]

test_that("LaPlagne: full maxima", {
  testthat::expect_equal(LaPlagnePrecipMaximaSeasonal$full_maxima,
                         full_maxima)
})
test_that("LaPlagne: partial maxima", {
  testthat::expect_equal(LaPlagnePrecipMaximaSeasonal$partial_maxima,
                         partial_maxima, ignore_attr = TRUE)
})

# Precipitation pseudo maxima

# 1998 (rows 303-453). The only full block: can donate to any partial block
w1998 <- min(which(LaPlagne$winter == 1998))
winter1998 <- which(rownames(LaPlagnePrecipMaximaSeasonal$pseudo_maxima)
                    == w1998)
pm2013 <- LaPlagnePrecipMaximaSeasonal$pseudo_maxima[winter1998, ]
# Mostly 41.6 (1998's maximum)
# 2013 is 22.8. Let's check it
rr2013 <- LaPlagne$rr24[LaPlagne$winter == 2013]
rr1998 <- LaPlagne$rr24[LaPlagne$winter == 1998]
rr1998[is.na(rr2013)] <- NA
pm2013check <- max(rr1998, na.rm = TRUE)
c2013 <- which(colnames(LaPlagnePrecipMaximaSeasonal$pseudo_maxima) == "2013")
test_that("LaPlagne: pseudo maxima, 1998 donating to 2013", {
  testthat::expect_equal(pm2013[c2013], pm2013check, ignore_attr = TRUE)
})

# Now use the snow data and a donor block not at the start of a winter

# Full maxima in 1998 and 2006
# LaPlagneSnowMaximaSeasonal$full_maxima

check_fn <- function(receiver_year, pick_block, seasonal = TRUE) {
  # Pick the start of a donor block
  block_start <- block_starts[pick_block]
  donor_winter <- LaPlagne[block_start, "winter"]

  # End of previous winter
  if (donor_winter > 1996) {
    end_previous_block <- max(which(LaPlagne$winter == donor_winter - 1))
  } else {
    end_previous_block <- 0
  }
  start_label <- as.numeric(rownames(LaPlagne))[block_start] -
    as.numeric(rownames(LaPlagne)[end_previous_block])

  donor_block <- LaPlagne[block_start:(block_start + 151 - 1), "ht_neige"]
  if (start_label == 1 || !seasonal) {
    the_names <- 1:151
  } else {
    the_names <- c(start_label:151, 1:(start_label - 1))
  }
  names(donor_block) <- the_names
  receiver_block <- LaPlagne[LaPlagne$winter == receiver_year, "ht_neige"]
  names(receiver_block) <- 1:length(receiver_block)

  where_receiver_NA <- names(receiver_block)[is.na(receiver_block)]
  where_make_donor_NA <- which(is.element(names(donor_block), where_receiver_NA))
  donor_block[where_make_donor_NA] <- NA

  res1 <- max(donor_block, na.rm = TRUE)
  names(res1) <- block_start
  res2 <- pseudo_for_receiver_year[pick_block]
  return(c(res1, res2))
}

# 1996
receiver_year <- 1996

# seasonal
receiver_column <- which(colnames(LaPlagneSnowMaximaSeasonal$pseudo_maxima)
                         == receiver_year)
sliding_block_nos <- as.numeric(
  which(!is.na(LaPlagneSnowMaximaSeasonal$pseudo_maxima[, receiver_column]))
)
pseudo_for_receiver_year <-
  LaPlagneSnowMaximaSeasonal$pseudo_maxima[sliding_block_nos, receiver_column]
block_starts <- as.numeric(names(pseudo_for_receiver_year))
# Look at block_starts to see how many values could be set for pick_block
res <- check_fn(receiver_year = receiver_year, pick_block = 1,
                seasonal = TRUE)
test_that("LaPlagne: sliding pseudo maxima, sliding block donating to 1996", {
  testthat::expect_equal(res[1], res[2])
})

# !seasonal
receiver_column <- which(colnames(LaPlagneSnowMaxima$pseudo_maxima)
                         == receiver_year)
sliding_block_nos <- as.numeric(
  which(!is.na(LaPlagneSnowMaxima$pseudo_maxima[, receiver_column]))
)
pseudo_for_receiver_year <-
  LaPlagneSnowMaxima$pseudo_maxima[sliding_block_nos, receiver_column]
block_starts <- as.numeric(names(pseudo_for_receiver_year))
# Look at block_starts to see how many values could be set for pick_block
res <- check_fn(receiver_year = receiver_year, pick_block = 1,
                seasonal = FALSE)
test_that("LaPlagne: sliding pseudo maxima, sliding block donating to 1996", {
  testthat::expect_equal(res[1], res[2])
})

# 2020
receiver_year <- 2020

# seasonal
receiver_column <- which(colnames(LaPlagneSnowMaximaSeasonal$pseudo_maxima)
                         == receiver_year)
sliding_block_nos <- as.numeric(
  which(!is.na(LaPlagneSnowMaximaSeasonal$pseudo_maxima[, receiver_column]))
)
pseudo_for_receiver_year <-
  LaPlagneSnowMaximaSeasonal$pseudo_maxima[sliding_block_nos, receiver_column]
block_starts <- as.numeric(names(pseudo_for_receiver_year))
# Look at block_starts to see how many values could be set for pick_block
res <- check_fn(receiver_year = receiver_year, pick_block = 17,
                seasonal = TRUE)
test_that("LaPlagne: sliding pseudo maxima, sliding block donating to 2020", {
  testthat::expect_equal(res[1], res[2])
})

# !seasonal
receiver_column <- which(colnames(LaPlagneSnowMaxima$pseudo_maxima)
                         == receiver_year)
sliding_block_nos <- as.numeric(
  which(!is.na(LaPlagneSnowMaxima$pseudo_maxima[, receiver_column]))
)
pseudo_for_receiver_year <-
  LaPlagneSnowMaxima$pseudo_maxima[sliding_block_nos, receiver_column]
block_starts <- as.numeric(names(pseudo_for_receiver_year))
# Look at block_starts to see how many values could be set for pick_block
res <- check_fn(receiver_year = receiver_year, pick_block = 17,
                seasonal = FALSE)
test_that("LaPlagne: sliding pseudo maxima, sliding block donating to 2020", {
  testthat::expect_equal(res[1], res[2])
})
