#' Methods for objects of class `"confint_gev"`
#'
#' Methods for objects of class `"confint_gev"` returned from
#'   [`confint.evmissing`].
#' @param x An object inheriting from class `"confint_gev"`, an object returned
#'   after a call to [`confint.evmissing`].
#' @param ... Further arguments. For `print.confint_gev` to pass arguments to
#'   [`print`]. For `plot.confint_gev` to pass graphical parameters to
#'   [`plot`][graphics::plot] to create the initial plot of the profile
#'   log-likelihood.
#' @details `print.confint_gev`. A numeric matrix with 2 columns giving the
#'   lower and upper confidence limits for the parameters specified by the
#'   argument `parm` in [`confint.evmissing`]. These columns are labelled as
#'   `(1-level)/2` and `1-(1-level)/2`, expressed as a percentage, by default
#'   `2.5%` and `97.5%`.
#'
#'   `plot.confint_gev`. A plot is produced of the profile log-likelihood for
#'   the parameter chosen by `parm`. Only the parameter values used to profile
#'   the log-likelihood in the call to [`confint.evmissing`] are included, so
#'   if `faster = TRUE` was used then the plot will not be of a smooth curve
#'   but will be triangular in the middle.
#' @return `print.confint_gev`: the argument `x` is returned, invisibly.
#'
#'   `plot.confint_gev`: a numeric vector containing the confidence limits for
#'   the parameter requested in `parm` is returned invisibly.
#' @section Examples: See [`evmissing_methods`].
#' @seealso [`gev_mle`] and [`evmissing_methods`].
#' @name confint_gev_methods
NULL
## NULL

# ============================ print.confint_gev ============================ #

#' Print method for objects of class `"confint_gev"`
#'
#' @rdname confint_gev_methods
#' @export
print.confint_gev <- function(x, ...) {
  print(x[1:nrow(x), , drop = FALSE], ...)
  return(invisible(x))
}

# ============================= plot.confint_gev ============================ #

#' Plot method for objects of class `"confint_gev"`
#'
#' @param parm A character scalar specifying the parameter for which
#'   a profile log-likelihood is plotted. Must be a single component of
#'   `c("mu", "sigma", "xi")`.
#' @param add A logical scalar. If `add = TRUE` then the plot is annotated with
#'   a horizontal line indicating the critical value for the profile
#'   log-likelihood used to calculate the confidence limits, vertical lines
#'   indicating the values of these limits and a legend stating the
#'   confidence interval.
#' @param digits An integer. Passed to [`signif`] to round the confidence
#'   limits in the legend, if `add = TRUE`. The confidence level is hard-coded
#'   to be expressed to 3 significant figures.
#' @rdname confint_gev_methods
#' @export
plot.confint_gev <- function(x, parm = c("mu", "sigma", "xi"), add = TRUE,
                             digits = 2, ...) {
  # Check that confint.evmissing() was called with profile = TRUE
  if (is.null(attr(x, "for_plot"))) {
    stop("Call confint.evmissing() with profile = TRUE")
  }
  # Check that parm was included in parm in the call to confint.evmissing()
  parameters_in_x <- attr(x, "dimnames")[[1]]
  if (any(!is.element(parm, parameters_in_x))) {
    stop(parm, " was not included in ''parm'' in the call to confint.evmissing.")
  }
  parm <- match.arg(parm)
  if (parm == "mu") {
    to_plot <- attr(x, "for_plot")$mu
    my_xlab <- "mu"
    limits <- x[rownames(x) == "mu", ]
  } else if (parm == "sigma") {
    to_plot <- attr(x, "for_plot")$sigma
    my_xlab <- "sigma"
    limits <- x[rownames(x) == "sigma", ]
  } else if (parm == "xi") {
    to_plot <- attr(x, "for_plot")$xi
    my_xlab <- "xi"
    limits <- x[rownames(x) == "xi", ]
  }
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
    legend_text <- paste0(signif(level, 3), "% CI: (", rlimits[1], ",",
                          rlimits[2], ")")
    graphics::legend("bottom", legend = legend_text)
  }
  return(invisible(limits))
}
