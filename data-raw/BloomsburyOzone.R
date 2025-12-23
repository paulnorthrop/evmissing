# Provide the Bloomsbury Ozone annual maxima data in the same format as
# BrestSurgeMaxima
# Also provide the raw data

# Source: [Data Selector - DEFRA UK Air - GOV.UK](https://uk-air.defra.gov.uk/data/data_selector)

# Read in the raw daily ozone values

BloomsburyOzone <- read.csv("data-raw/OzoneBloomsbury.csv")
# Make BloomsburyOzone$Date a proper date object
BloomsburyOzone$Date <- as.Date(BloomsburyOzone$Date, "%d/%m/%Y")
head(BloomsburyOzone)
tail(BloomsburyOzone)

# Calculate the annual maxima, notNA, block size n and block indicator

maxima <- tapply(BloomsburyOzone$Ozone, INDEX = BloomsburyOzone$Year,
                 FUN = max, na.rm = TRUE)
notNA <- tapply(BloomsburyOzone$Ozone, INDEX = BloomsburyOzone$Year,
                FUN = function(x) sum(!is.na(x)))
n <- evmissing::days_in_year(names(maxima))
block <- 1:length(maxima)

# Create the data frame

BloomsburyOzoneMaxima <- data.frame(maxima = maxima,notNA = notNA, n = n,
                                    block = block)
rownames(BloomsburyOzoneMaxima) <- names(maxima)

# Fit using the raw data
fit1 <- evmissing::gev_mle(BloomsburyOzone$Ozone, block = BloomsburyOzone$Year)
fit2 <- evmissing::gev_mle(BloomsburyOzone$Ozone, block = BloomsburyOzone$Year,
                        adjust = FALSE)
coef(fit1)
coef(fit2)

# Fit using the data frame of maxima
fit1 <- evmissing::gev_mle(BloomsburyOzoneMaxima)
fit2 <- evmissing::gev_mle(BloomsburyOzoneMaxima, adjust = FALSE)
coef(fit1)
coef(fit2)

# Create the package data
usethis::use_data(BloomsburyOzoneMaxima, overwrite = TRUE)
usethis::use_data(BloomsburyOzone, overwrite = TRUE)
