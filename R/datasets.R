#' Annual maxima sea surge heights at Brest, France
#'
#' Annual maxima of sea surge heights near high tide at Brest tide gauge
#' station (France) for the years 1846-2007 inclusive.
#'
#' @format `BrestSurgeMaxima` is a data frame with 162 rows (years 1846 to
#'   2007) and the 4 variables:
#'
#' * `maxima`: annual maximum surge height at high tide in cm.
#' * `notNA` : the number of days of the year for which raw data were available.
#' * `n` : the number of days in the year (365 or 366).
#' * `block` : a block number of 1 for year 1846 through to 162 for year 2007.
#'
#' The row names of `BrestSurgeMaxima` are the years `1946:2007`.
#'
#' @note The raw data are missing for approximately \eqn{9\%} of the days.
#' The data were declustered by the original providers in order to provide a
#' series of independent surge heights at high tide. Specifically, these
#' surge heights are separated by at least two days. A correction was applied
#' to account for trend in the sea-level over the observation period.
#' Although the declustering of the data means that the effective block size is
#' smaller than `n`, it may be reasonable to suppose that the proportion
#' `notNA/n` of non-missing values provides a useful measure of the extent to
#' which the size of an annual maximum is likely to be affected by missingness.
#'
#' @source The dataset `Brest` in the `Renext` R package, specifically
#'   `Brest$OTdata` and `Brest$OTmissing`. Originally, the source was
#'   \url{https://data.shom.fr/}.
#' @references Deville Y. and Bardet L. (2023). Renext: Renewal Method for
#'   Extreme Values Extrapolation. R package version 3.1-4.
#'   \doi{https://doi.org/10.32614/CRAN.package.Renext}
#' @seealso
#'
#'  * [`BrestSurgeMissing`]: numbers of missing values in each month.
#'  * [`BrestSurgeDays`]: Number of days per month in 1846-2007.
#'
#' @examples
#' head(BrestSurgeMaxima)
#'
#' # Time series plot of annual maxima surges
#' plot(rownames(BrestSurgeMaxima), BrestSurgeMaxima$maxima,
#'      ylab = "surge (cm)", xlab = "year", pch = 16)
#'
#' # Time series plot of proportion of non-missing days
#' plot(rownames(BrestSurgeMaxima), BrestSurgeMaxima$notNA / BrestSurgeMaxima$n,
#'      ylab = "proportion of non-missing days", xlab = "year", pch = 16)
#'
#' # Plot surges against the proportion of non-missing days
#' plot(BrestSurgeMaxima$notNA / BrestSurgeMaxima$n, BrestSurgeMaxima$maxima,
#'      ylab = "surge (cm)", xlab = "proportion of non-missing days", pch = 16)
"BrestSurgeMaxima"

#' Missing values in sea surge heights at Brest, France
#'
#' Numbers of missing values in each month of the Brest sea surge heights data
#' [`BrestSurgeMaxima`].
#'
#' @format `BrestSurgeMissing` is a data frame with 162 rows (years 1846 to
#'   2007) and the 12 variables (one for each month of the year). Each value
#'   in the data frame gives the number of days for which the surge height data
#'   were missing in the month in question.
#'
#' The row names of `BrestSurgeMaxima` are the years `1946:2007` and the column
#' names are the abbreviated names of the months.
#'
#' @source The dataset `Brest` in the `Renext` R package, specifically
#'   `Brest$OTmissing`. Originally, the source was
#'   \url{https://data.shom.fr/}.
#' @references Deville Y. and Bardet L. (2023). Renext: Renewal Method for
#'   Extreme Values Extrapolation. R package version 3.1-4.
#'   \doi{https://doi.org/10.32614/CRAN.package.Renext}
#' @seealso
#'
#'  * [`BrestSurgeMaxima`]: Annual maxima surge heights at Brest, France.
#'  * [`BrestSurgeDays`]: Number of days per month in 1846-2007.
#'
#' @examples
#' head(BrestSurgeMissing)
#'
#' # Proportion of missing values by year
#' propn_year <- rowSums(BrestSurgeMissing) /
#'   days_in_year(rownames(BrestSurgeMissing))
#' plot(rownames(BrestSurgeMissing), propn_year,
#'      ylab = "proportion of missing values", xlab = "year", pch = 16)
#'
#' # Proportion of missing values by year and month
#' propn_year_month <- BrestSurgeMissing / BrestSurgeDays
#'
#' # Proportion of missing values by month
#' plot(1:12, colMeans(propn_year_month), axes = FALSE,
#'         ylab = "proportion of missing values", xlab = "month", pch = 16)
#' axis(1, at = 1:12, labels = 1:12)
#' axis(2)
#' box()
"BrestSurgeMissing"

#' Number of days per month in 1846-2007
#'
#' Number of days in each month relevant to the Brest sea surge heights data
#' [`BrestSurgeMaxima`].
#'
#' @format `BrestSurgeDays` is a data frame with 162 rows (years 1846 to
#'   2007) and the 12 variables (one for each month of the year). Each value
#'   in the data frame gives the number of days in the month in question.
#'
#' The row names of `BrestSurgeMaxima` are the years `1946:2007` and the column
#' names are the abbreviated names of the months.
#'
#' @seealso
#'
#'  * [`BrestSurgeMaxima`]: Annual maxima surge heights at Brest, France.
#'  * [`BrestSurgeMissing`]: numbers of missing values in each month.
#'
#' @examples
#' head(BrestSurgeDays)
"BrestSurgeDays"

#' Annual maxima ozone levels at Bloomsbury, UK
#'
#' Annual maxima of daily maximum ozone levels at Bloomsbury in London (UK)
#' for the years 1992-2024 inclusive.
#'
#' @format `BloomsburyOzoneMaxima` is a data frame with 33 rows (years 1992 to
#'   2024) and the 4 variables:
#'
#' * `maxima`: annual maximum ozone level in \eqn{\mu}g/m\eqn{^3}.
#' * `notNA` : the number of days of the year for which raw data were available.
#' * `n` : the number of days in the year (365 or 366).
#' * `block` : a block number of 1 for year 1992 through to 33 for year 2024.
#'
#' The row names of `BloomsburyOzoneMaxima` are the years `1992:2024`.
#' The raw data are missing for approximately \eqn{5\%} of the days.
#'
#' @source The Department for Environment Food and Rural Affair (DEFRA).
#'   The London Bloomsbury monitoring site at
#'   the [UK-AIR](https://uk-air.defra.gov.uk/) database
#'   [Data Selector](https://uk-air.defra.gov.uk/data/data_selector).
#' @seealso [`BloomsburyOzone`] for the raw time series.
#' @examples
#' head(BloomsburyOzoneMaxima)
#'
#' # Time series plot of annual maxima ozone levels
#' plot(rownames(BloomsburyOzoneMaxima), BloomsburyOzoneMaxima$maxima,
#'      ylab = "ozone (micrograms / metre cubed)", xlab = "year", pch = 16)
#'
#' # Time series plot of proportion of non-missing days
#' plot(rownames(BloomsburyOzoneMaxima),
#'      BloomsburyOzoneMaxima$notNA / BloomsburyOzoneMaxima$n,
#'      ylab = "proportion of non-missing days", xlab = "year", pch = 16)
#'
#' # Plot ozone levels against the proportion of non-missing days
#' plot(BloomsburyOzoneMaxima$notNA / BloomsburyOzoneMaxima$n,
#'      BloomsburyOzoneMaxima$maxima,
#'      ylab = "ozone (micrograms / metre cubed)",
#'      xlab = "proportion of non-missing days", pch = 16)
"BloomsburyOzoneMaxima"

#' Ozone levels at Bloomsbury, UK
#'
#' Daily maximum ozone levels at Bloomsbury in London (UK) for the years
#' 1992-2024 inclusive.
#'
#' @format `BloomsburyOzone` is a data frame with 12054 rows and the 3 variables:
#'
#' * `Date`: with class `"Date"` in the format `YYYY-MM-DD`.
#' * `Year`: Values in 1992-2024.
#' * `Ozone`: daily maximum ozone level in \eqn{\mu}g/m\eqn{^3}.
#'
#' @source The Department for Environment Food and Rural Affair (DEFRA).
#'   The London Bloomsbury monitoring site at
#'   the [UK-AIR](https://uk-air.defra.gov.uk/) database
#'   [Data Selector](https://uk-air.defra.gov.uk/data/data_selector).
#' @seealso [`BloomsburyOzoneMaxima`] for the annual maxima and numbers of
#'   missing values per year.
#' @examples
#' head(BloomsburyOzone)
#'
#' # Time series plot of annual maxima ozone levels
#' plot(BloomsburyOzone$Date, BloomsburyOzone$Ozone, xlab = "year",
#'      ylab = "ozone (micrograms / metre cubed)", pch = 16)
"BloomsburyOzone"

#' Annual maxima ozone levels at Plymouth, UK
#'
#' Annual maxima of daily maximum ozone levels at Plymouth in Devon (UK)
#' for the years 1998-2024 inclusive.
#'
#' @format `PlymouthOzoneMaxima` is a data frame with 27 rows (years 1998 to
#'   2024) and the 4 variables:
#'
#' * `maxima`: annual maximum ozone level in \eqn{\mu}g/m\eqn{^3}.
#' * `notNA` : the number of days of the year for which raw data were available.
#' * `n` : the number of days in the year (365 or 366).
#' * `block` : a block number of 1 for year 1998 through to 27 for year 2024.
#'
#' The row names of `PlymouthOzoneMaxima` are the years `1998:2024`.
#' The raw data are missing for approximately \eqn{10\%} of the days.
#'
#' @source The Department for Environment Food and Rural Affair (DEFRA).
#'   The Plymouth Centre monitoring site at
#'   the [UK-AIR](https://uk-air.defra.gov.uk/) database
#'   [Data Selector](https://uk-air.defra.gov.uk/data/data_selector).
#' @seealso [`PlymouthOzone`] for the raw time series.
#' @examples
#' head(PlymouthOzoneMaxima)
#'
#' # Time series plot of annual maxima ozone levels
#' plot(rownames(PlymouthOzoneMaxima), PlymouthOzoneMaxima$maxima,
#'      ylab = "ozone (micrograms / metre cubed)", xlab = "year", pch = 16)
#'
#' # Time series plot of proportion of non-missing days
#' plot(rownames(PlymouthOzoneMaxima),
#'      PlymouthOzoneMaxima$notNA / PlymouthOzoneMaxima$n,
#'      ylab = "proportion of non-missing days", xlab = "year", pch = 16)
#'
#' # Plot ozone levels against the proportion of non-missing days
#' plot(PlymouthOzoneMaxima$notNA / PlymouthOzoneMaxima$n,
#'      PlymouthOzoneMaxima$maxima,
#'      ylab = "ozone (micrograms / metre cubed)",
#'      xlab = "proportion of non-missing days", pch = 16)
"PlymouthOzoneMaxima"

#' Ozone levels at Plymouth, UK
#'
#' Daily maximum ozone levels at Plymouth in London (UK) for the years
#' 1998-2024 inclusive.
#'
#' @format `PlymouthOzone` is a data frame with 9862 rows and the 3 variables:
#'
#' * `Date`: with class `"Date"` in the format `YYYY-MM-DD`.
#' * `Year`: Values in 1998-2024.
#' * `Ozone`: daily maximum ozone level in \eqn{\mu}g/m\eqn{^3}.
#'
#' @source The Department for Environment Food and Rural Affair (DEFRA).
#'   The Plymouth Centre monitoring site at
#'   the [UK-AIR](https://uk-air.defra.gov.uk/) database
#'   [Data Selector](https://uk-air.defra.gov.uk/data/data_selector).
#' @seealso [`PlymouthOzoneMaxima`] for the annual maxima and numbers of
#'   missing values per year.
#' @examples
#' head(PlymouthOzone)
#'
#' # Time series plot of annual maxima ozone levels
#' plot(PlymouthOzone$Date, PlymouthOzone$Ozone, xlab = "year",
#'      ylab = "ozone (micrograms / metre cubed)", pch = 16)
"PlymouthOzone"
