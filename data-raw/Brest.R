# We modify the `Brest` dataset provided in the `Renext` package.
# The aim is to provide a dataset of hourly surge heights with explicit `NA`
# values for hours in which the underlying data value was missing.
# We also provide a smaller dataset containing the annual maxima and the
# numbers of non-missing values in each block, consistent with an object
# returned from `block_maxima()`.

# Validated hourly data for surge
# Global REFMAR DOI : http://dx.doi.org/10.17183/REFMAR

library(Renext)

## BrestSurgeMaxima

# 1. The annual maxima

# Extract the dates of the surges and then their years
dates <- Brest$OTdata$date
year <- format(dates, format="%Y")
surge <- Brest$OTdata$Surge
# Create a table that contains a 0 for any year with no exceedances of 30cm.
# We do not include 2008. The end date in Brest$OTinfo is "2009-01-01 GMT", but
# the last row in Brest$OTmissing is 2008-01-10 to 2009-01-01 and there are no
# entries for Surge in Brest$OTdata for 2008 for the 9 days of 2008 for which
# there could have been data. While we could have added
#   2008 NA 9 266 163
# to the end of these data it seems clearer not to do this
year_range <- 1846:2007
x <- table(c(year, year_range)) - 1
# Where are the 0s? That is, where are the years that have no surge data
which(x == 0)
# Are these owing to whole years of missing data?
Brest$OTmissing
# 1858: yes
# 1938: yes
# 1945-1951: yes
# For the sake of completeness, we will include these missing years in the data
surge <- c(surge, rep(NA, length(which(x == 0))))
year <- c(year, names(which(x == 0)))
# Calculate the surge height annual maxima
surgeYear <- tapply(X = surge, INDEX = year, FUN = max)

# 2. The number of non-missing values in each year

days_between_dates <- function(start_date, end_date) {
  days <- as.numeric(difftime(end_date, start_date, units = "days")) + 1
  return(days)
}

start_year <- format(Brest$OTmissing$start, format="%Y")
end_year <- format(Brest$OTmissing$end, format="%Y")

number_of_missing_days <- function(year) {
  # Set up a counter of the number of non-missing days
  notNA <- 0
  # There are 4 cases
  # 1. start_year = end_year = year
  case1 <- start_year == year & end_year == year
  if (any(case1)) {
    start_date <- Brest$OTmissing[which(case1), 1]
    end_date <- Brest$OTmissing[which(case1), 2]
    notNA <- notNA + sum(days_between_dates(start_date, end_date))
  }
  # 2. start_year = year, end_year > year
  case2 <- start_year == year & end_year > year
  if (any(case2)) {
    start_date <- Brest$OTmissing[which(case2), 1]
    end_date <- as.Date(paste0(year, "-12-31"))
    notNA <- notNA + sum(days_between_dates(start_date, end_date))
  }
  # 3. start_year < year, end_year = year
  case3 <- start_year < year & end_year == year
  if (any(case3)) {
    start_date <- as.Date(paste0(year, "-01-01"))
    end_date <- Brest$OTmissing[which(case3), 2]
    notNA <- notNA + sum(days_between_dates(start_date, end_date))
  }
  # 4. start_year < year, end_year > year
  case4 <- start_year < year & end_year > year
  if (any(case4)) {
    start_date <- as.Date(paste0(year, "-01-01"))
    end_date <- as.Date(paste0(year, "-12-31"))
    notNA <- notNA + sum(days_between_dates(start_date, end_date))
  }
  return(notNA)
}

numbersOfNAs <- vapply(year_range, number_of_missing_days, 0)

# 3. The largest number of values that could have been observed in each year

days_in_year <- function(year) {
  start_date <- as.Date(paste0(year, "-01-01"))
  end_date <- as.Date(paste0(year, "-12-31"))
  days_in_year <- as.numeric(difftime(end_date, start_date, units = "days")) + 1
  return(days_in_year)
}

n <- days_in_year(names(surgeYear))

# Infer the number of non-missing days

notNA <- n - numbersOfNAs

# 4. Create the data frame

BrestSurgeMaxima <- data.frame(maxima = surgeYear,
                               notNA = notNA,
                               n = n,
                               block = 1:length(surgeYear))
rownames(BrestSurgeMaxima) <- names(surgeYear)

## Plots

# Time series plot of annual maxima surges
plot(rownames(BrestSurgeMaxima), BrestSurgeMaxima$maxima,
     ylab = "surge (cm)", xlab = "year", pch = 16)

# Time series plot of proportion of non-missing days
plot(rownames(BrestSurgeMaxima), BrestSurgeMaxima$notNA / BrestSurgeMaxima$n,
     ylab = "surge (cm)", xlab = "year", pch = 16)

# Plot surges against the proportion of non-missing days
plot(BrestSurgeMaxima$notNA / BrestSurgeMaxima$n, BrestSurgeMaxima$maxima,
     ylab = "surge (cm)", xlab = "proportion of non-missing days", pch = 16)

# Create the package data
usethis::use_data(BrestSurgeMaxima, overwrite = TRUE)

## BrestSurgeMissing

start_month <- format(Brest$OTmissing$start, format="%m")
end_month <- format(Brest$OTmissing$end, format="%m")

days_in_month <- function(year, month) {
  month <- as.numeric(month)
  next_month <- ifelse(month != 12,
                       paste(year, month + 1, "01", sep = "-"),
                       paste(year + 1, "01", "01", sep = "-"))
  next_month <- as.Date(next_month)
  last_day_of_month <- next_month - 1
  return(as.integer(format(last_day_of_month, "%d")))
}

number_of_missing_days_year_month <- function(year, month) {

  # Set up a counter of the number of non-missing days
  notNA <- 0
  # There are 4 cases
  # 1. Start and end within a single month
  #    start_year = end_year = year and start_month = end_month = month
  case1 <- start_year == year & end_year == year &
    start_month == month & end_month == month
  if (any(case1)) {
    start_date <- Brest$OTmissing[which(case1), 1]
    end_date <- Brest$OTmissing[which(case1), 2]
    notNA <- notNA + sum(days_between_dates(start_date, end_date))
  }
  # 2. Starts within a month and continues beyond the end of that month
  #    start_year = year, start_month = month and
  #    end_year > year or end_month != month
  case2 <- start_year == year & start_month == month &
    (end_year > year | end_month != month)
  if (any(case2)) {
    start_date <- Brest$OTmissing[which(case2), 1]
    end_date <- as.Date(paste(year, month, days_in_month(year, month),
                              sep = "-"))
    notNA <- notNA + sum(days_between_dates(start_date, end_date))
  }
  # 3. Starts before the start of a month stops before the end of the month
  #    (start_year = year and start_month < month) or (start_year < year) and
  #    (end_year = year and end_month > month) or (end_year > year)
  case3a <- (start_year == year & start_month < month) | (start_year < year)
  case3b <- end_year == year & end_month == month
  case3 <- case3a & case3b
  if (any(case3)) {
    start_date <- as.Date(paste(year, month, "01", sep = "-"))
    end_date <- Brest$OTmissing[which(case3), 2]
    notNA <- notNA + sum(days_between_dates(start_date, end_date))
  }
  # 4. Start before the start of a month and ends after the end of the month
  #    (start_year = year and start_month < month) or (start_year < year) and
  #    (end_year = year and end_month > month) or (end_year > year)
  case4a <- (start_year == year & start_month < month) | (start_year < year)
  case4b <- (end_year == year & end_month > month) | (end_year > year)
  case4 <- case4a & case4b
  if (any(case4)) {
    notNA <- notNA + sum(days_in_month(year, month))
  }
  return(notNA)
}

year_month_vec <- Vectorize(number_of_missing_days_year_month)
month_range <- c(paste0(0, 1:9), 10:12)
BrestSurgeMissing <- outer(year_range, month_range, FUN = year_month_vec)
rownames(BrestSurgeMissing) <- year_range
colnames(BrestSurgeMissing) <- month.abb
BrestSurgeMissing <- as.data.frame(BrestSurgeMissing)

# Create a data frame containing the number of days for each month in the data

BrestSurgeDays <- outer(year_range, month_range, FUN = days_in_month)
rownames(BrestSurgeDays) <- year_range
colnames(BrestSurgeDays) <- month.abb
BrestSurgeDays <- as.data.frame(BrestSurgeDays)

# Create the package data
usethis::use_data(BrestSurgeDays, overwrite = TRUE)

# Plots

# Proportion of missing values by year
propn_year <- rowSums(BrestSurgeMissing) /
  days_in_year(rownames(BrestSurgeMissing))
plot(rownames(BrestSurgeMissing), propn_year,
     ylab = "proportion of missing values", xlab = "year", pch = 16)

# Proportion of missing values by year and month
propn_year_month <- BrestSurgeMissing / BrestSurgeDays

# Proportion of missing values by month
plot(1:12, colMeans(propn_year_month), axes = FALSE,
     ylab = "proportion of missing values", xlab = "month", pch = 16)
axis(1, at = 1:12, labels = 1:12)
axis(2)
box()

# Create the package data
usethis::use_data(BrestSurgeMissing, overwrite = TRUE)
