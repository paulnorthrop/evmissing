#' Plot method for objects of class `"sliding_block_maxima"`
#'
#' Plot method for objects of class `"sliding_block_maxima"` returned from
#'   [`sliding_block_maxima_ts`], for example [`PlymouthOzoneSlidingMaxima`].
#' @param x An object inheriting from class `"sliding_block_maxima"`, an
#'   object from a call to [`sliding_block_maxima_ts`].
#' @param which If `which = 1` then the sliding block maxima are plotted
#'   against block number. If `which = 2` then the sliding block maxima are
#'   plotted against the proportion of non-missing raw values.
#' @param ... Further arguments to [`plot`][graphics::plot].
#' @details When `which = 1` we obtain a time series plot in which there are
#'   periods where the value does not change. When `which = 2` we expect to
#'   see that the sliding block maximum tends to be smaller for blocks with a
#'   larger proportionof missing values.
#' @return Nothing is returned.
#' @examples
#' # Time series plots of sliding block maxima
#' plot(PlymouthOzoneSlidingMaxima)
#'
#' # Plot maxima against the proportion of non-missing daily values
#' plot(PlymouthOzoneSlidingMaxima, which = 2)
#' @export
plot.sliding_block_maxima_ts <- function(x, which = 1, ...) {
  # Create the data to plot
  yvals <- x$maxima
  the_ylab <- "sliding block maximum"
  if (which == 1) {
    xvals <- 1:length(x$maxima)
    the_type <- "l"
    the_xlab <- "block"
    the_pch <- 1
  } else if (which == 2) {
    xvals <- x$notNA / x$n
    the_type <- "p"
    the_xlab <- "proportion of non-missing raw values"
    the_pch <- 16
  } else {
    stop("''which' must be equal to 1 or 2.")
  }
  plot_fn <- function(x, y, ..., xlab = the_xlab, ylab = the_ylab,
                      lwd = 2, type = the_type, pch = the_pch) {
    graphics::plot(x, y, ..., xlab = xlab, ylab = ylab,
                   lwd = lwd, type = type, pch = pch)
  }
  plot_fn(x = xvals, y = yvals, ...)
  return(invisible())
}
