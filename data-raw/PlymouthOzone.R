# Provide the Plymouth Ozone annual maxima data in the same format as
# BrestSurgeMaxima
# (Also provide the raw data)

# Source: [Data Selector - DEFRA UK Air - GOV.UK](https://uk-air.defra.gov.uk/data/data_selector)

## Read in the raw daily ozone values

PlymouthOzone <- read.csv("data-raw/OzonePlymouth.csv", skip = 10,
                          na.strings = "No data")
# Make PlymouthOzone$Date a proper date object
PlymouthOzone$Date <- as.Date(PlymouthOzone$Date, "%d/%m/%Y")
head(PlymouthOzone)
tail(PlymouthOzone)

# Add a column for year
year <- as.numeric(format(PlymouthOzone$Date, format="%Y"))

PlymouthOzone <- data.frame(Date = PlymouthOzone$Date, Year = year,
                            Ozone = PlymouthOzone$Ozone)
head(PlymouthOzone)
tail(PlymouthOzone)

# 1997 and 2025 are incomplete
# We could retain them, but 1997 has only 2 months and 2025 almost 5 months
# Let's remove 1997 and 2025

retain <- PlymouthOzone$Year > 1997 & PlymouthOzone$Year < 2025
PlymouthOzone <- PlymouthOzone[retain, ]
head(PlymouthOzone)
tail(PlymouthOzone)

# Calculate the annual maxima, notNA, block size n and block indicator

maxima <- tapply(PlymouthOzone$Ozone, INDEX = PlymouthOzone$Year,
                 FUN = max, na.rm = TRUE)
notNA <- tapply(PlymouthOzone$Ozone, INDEX = PlymouthOzone$Year,
                FUN = function(x) sum(!is.na(x)))
n <- evmissing::days_in_year(names(maxima))
block <- 1:length(maxima)

# Create the data frame

PlymouthOzoneMaxima <- data.frame(maxima = maxima,notNA = notNA, n = n,
                                  block = block)
rownames(PlymouthOzoneMaxima) <- names(maxima)

# Fit using the raw data
fit1 <- evmissing::gev_mle(PlymouthOzone$Ozone, block = PlymouthOzone$Year)
fit2 <- evmissing::gev_mle(PlymouthOzone$Ozone, block = PlymouthOzone$Year,
                        adjust = FALSE)
coef(fit1)
coef(fit2)

# Fit using the data frame of maxima
fit1 <- evmissing::gev_mle(PlymouthOzoneMaxima)
fit2 <- evmissing::gev_mle(PlymouthOzoneMaxima, adjust = FALSE)
coef(fit1)
coef(fit2)

# Create the package data
usethis::use_data(PlymouthOzoneMaxima, overwrite = TRUE)
usethis::use_data(PlymouthOzone, overwrite = TRUE)

# Explore the distribution of lengths of missing and non-missing episodes

x <- rle(is.na(PlymouthOzone$Ozone))
missing_periods <- x$lengths[x$values]
non_missing_periods <- x$lengths[!x$values]
tab_missing <- table(c(missing_periods, 1:max(missing_periods))) - 1
barplot(tab_missing)
tab_non_missing <- table(c(non_missing_periods, 1:max(non_missing_periods))) - 1
barplot(tab_non_missing)
plot(non_missing_periods, c(missing_periods, NA),
     xlab = "length of prior presence period",
     ylab = "length of following missing period")
plot(non_missing_periods, c(missing_periods, NA),
     xlab = "length of prior presence period",
     ylab = "length of following missing period", log = "xy")
plot(missing_periods, non_missing_periods[-1],
     xlab = "length of prior missing period",
     ylab = "length of following presence period", log = "xy")

# It seems reasonable to model the missing/non-missing process as an
# alternating renewal process using either
# (a) empirical distributions from these data, or
# (b) fits of suitable probability distributions to these data

# Alternatively, I could base the patterns of missingness exactly on
# those observed in PlymouthOzone

# I should also explore whether the onset of a period of missing values seems
# related at all to the most recent non-missing value of ozone.

beforeNA <- function(x) {
  # Find all NAs
  whereNA <- which(is.na(x))
  # Remove whereNA = 1, if this occurs
  whereNA <- whereNA[whereNA > 1]
  beforeNA <- whereNA - 1
  # Extract values immediately before each NA
  before <- x[beforeNA]
  before <- before[!is.na(before)]
  elsewhere <- x[-beforeNA]
  elsewhere <- elsewhere[!is.na(elsewhere)]
  # Return only non-missing such values
  return(list(before = before, elsewhere = elsewhere))
}

# It seems not. The sample quartiles are very similar. The range of the
# "elsewhere" values is larger but this is because there are 8770 of these and
# only 82 "before" values.
b <- beforeNA(PlymouthOzone$Ozone)
summary(b$before)
summary(b$elsewhere)
