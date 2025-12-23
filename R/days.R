#' Days in a year or in a month
#'
#' Returns the number of days in each of a vector of years or months.
#'
#' @param year An integer vector. The years of interest.
#' @param month An integer vector. A subset of `1:12`.The months of interest.
#' @details The length of the output vector is equal to the length of `month`.
#'   The argument `year` is recycled to the length of the output vector if
#'   necessary.
#'
#' @return A numeric vector of the numbers of days in each of the years in
#'   `year` or the months specified by `year` and `month`.
#' @examples
#' days_in_year(1999:2025)
#'
#' days_in_month(2024, 1:12)
#' days_in_month(2025, 1:12)
#' days_in_month(2024:2025, 1:3)
#' @name days
NULL


#' @rdname days
#' @export
days_in_year <- function(year) {
  start_date <- as.Date(paste0(year, "-01-01"))
  end_date <- as.Date(paste0(year, "-12-31"))
  return(as.numeric(difftime(end_date, start_date, units = "days")) + 1)
}

#' @rdname days
#' @export
days_in_month <- function(year, month) {
  month <- as.numeric(month)
  by_year_fun <- function(year) {
    next_month <- ifelse(month != 12,
                         paste(year, month + 1, "01", sep = "-"),
                         paste(year + 1, "01", "01", sep = "-"))
    next_month <- as.Date(next_month)
    last_day_of_month <- next_month - 1
    return(last_day_of_month)
  }
  last_day_of_month <- lapply(year, FUN = by_year_fun)
  last_day_of_month <- do.call(c, last_day_of_month)
  return(as.integer(format(last_day_of_month, "%d")))
}
