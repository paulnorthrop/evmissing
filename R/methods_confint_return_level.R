#' Methods for objects of class `"confint_return_level"`
#'
#' Methods for objects of class `"confint_return_level"` returned from
#' [`confint.return_level`].
#' @param x An object inheriting from class `"confint_return_level"`, a result
#'   of a call to [`confint.return_level`].
#' @param ... Further arguments. For `print.confint_return_level` to pass
#'   arguments to [`print`]). For `plot.confint_return_level` to pass graphical
#'   parameters to [`plot`][graphics::plot] to create the initial plot of the
#'   profile log-likelihood.
#' @details `print.confint_return_level`. A numeric matrix with 2 columns
#'   giving the lower and upper confidence limits for the parameters specified
#'   by the argument `parm` in [`confint.return_level`]. These columns are
#'   labelled as `(1-level)/2` and `1-(1-level)/2`, expressed as a percentage,
#'   by default `2.5%` and `97.5%`.
#'
#'   `plot.confint.return_level`. A plot is produced of the profile log-likelihood for
#'   the parameter chosen by `parm`.
#' @return `print.confint_return_level`: the argument `x` is returned, invisibly.
#'
#'   `plot.confint_return_level`: a numeric vector containing the confidence
#'   interval for the return level chosen for the plot.
#' @section Examples: See [`return_level_methods`].
#' @seealso [`gev_mle`], [`gev_return`] and [`return_level_methods`].
#' @name confint_return_level_methods
NULL
## NULL

# ======================== print.confint_return_level ======================= #

#' Print method for objects of class `"confint_return_level"`
#'
#' @rdname confint_return_level_methods
#' @export
print.confint_return_level <- function(x, ...) {
  print(x[1:nrow(x), , drop = FALSE], ...)
  return(invisible(x))
}

# ========================= plot.confint_return_level ======================= #

#' Plot method for objects of class `"confint_return_level"`
#'
#' @param parm An integer scalar. For which component, that is, which return
#'   level, in `x` we require a confidence interval.
#' @param add A logical scalar. If `add = TRUE` then the plot is annotated with
#'   a horizontal line indicating the critical value for the profile
#'   log-likelihood used to calculate the confidence limits, vertical lines
#'   indicating the values of these limits and a legend stating the
#'   confidence interval.
#' @param digits An integer. Passed to [`signif`] to round the confidence
#'   limits in the legend, if `add = TRUE`.
#' @rdname confint_return_level_methods
#' @export
plot.confint_return_level <- function(x, parm = 1, add = TRUE, digits = 2,
                                      ...) {
  # Check that confint.return_level() was called with profile = TRUE
  if (is.null(attr(x, "for_plot"))) {
    stop("Call confint.return_level() with profile = TRUE")
  }

  # Check parm
  if (length(parm) > 1) {
    stop("parm must be an integer scalar")
  }
  p_message <- 1:nrow(x)
  parm_values <- 1:nrow(x)
  if (!all(is.element(parm, parm_values))) {
    stop("''parm'' must be one of {", paste(parm_values, collapse = ","),
         "}")
  }
  parm <- parm[order(parm)]
  limits <- x[parm, ]
  my_xlab <- rownames(x)[parm]
  to_plot <- attr(x, "for_plot")[[parm]]
  crit <- attr(x, "crit")
  plot_fn <- function(x, ..., xlab = my_xlab, ylab = "profile log-likelihood",
                      lwd = 2, type = "l") {
    graphics::plot(x, ..., xlab = xlab, ylab = ylab, lwd = lwd, type = type)
  }
  plot_fn(to_plot, ...)
  if (add) {
    user_args <- list(...)
    graphics::abline(h = crit, lty = 2)
    graphics::abline(v = limits, lty = 2)
    level <- attr(x, "level") * 100
    rlimits <- signif(limits, digits)
    legend_text <- paste0(level, "% CI: (", rlimits[1], ",", rlimits[2], ")")
    graphics::legend("bottom", legend = legend_text)
  }
  return(invisible(limits))
}
