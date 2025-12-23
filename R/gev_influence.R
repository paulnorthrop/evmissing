#' GEV influence curves
#'
#' Calculates influence function curves for maximum likelihood estimators of
#' Generalised Extreme Value (GEV) parameters.
#'
#' @param z A numeric vector. Values of normal quantiles \eqn{z} at which to
#'  calculate the GEV influence function. See **Details**.
#' @param mu,sigma,xi Numeric scalars supplying the values of the GEV
#'   parameters \eqn{\mu}, \eqn{\sigma} and \eqn{\xi}.
#' @details An influence function measures the effect on a parameter estimator
#'   of changing one observation in a sample. See Hampel (2005) and, in the
#'   context of extreme value analyses of threshold exceedances,
#'   Davison and Smith (1990).
#'   Let \eqn{\theta = (\mu, \sigma, \xi)}. The GEV influence curve is defined,
#'   as a function of an observation \eqn{y}, by
#'   \eqn{IC(y) = \{IC_{\mu}(y), IC_{\sigma}(y), IC_{\xi}(y)\} =
#'   i_\theta^{-1} {\rm d}\ell(y; \theta)/{\rm d} \theta},
#'   where \eqn{\ell(y; \theta)} is the GEV log-likelihood function and
#'   \eqn{i_\theta^{-1}} is the GEV expected information matrix.
#'   The value of \eqn{y} is related to \eqn{z} by
#'   \eqn{y = G^{-1}\{\Phi(z)\}}, where \eqn{\Phi} and \eqn{G} are the
#'   distribution functions of a standard normal and
#'   GEV(\eqn{\mu, \sigma, \xi}) distribution, respectively. Plotting influence
#'   curves on the standard normal `z` scale can help to aid interpretation.
#'
#'   The value of \eqn{IC(y)} is invariant to the value of \eqn{\mu}.
#'   For a given value of \eqn{\xi}, the values of \eqn{IC_{\mu}(y)} and
#'   \eqn{IC_{\sigma}(y)} scale with the scale parameter \eqn{\sigma}, that is,
#'   doubling \eqn{\sigma} doubles \eqn{IC_{\mu}(y)} and \eqn{IC_{\sigma}(y).}
#'   Typically, interest focuses on the shape parameter \eqn{\xi}, but if the
#'   input scale parameter \eqn{\sigma} is large then this may hide the
#'   influence of \eqn{y} on \eqn{\xi}. Therefore, the default setting of
#'   `plot.gev_influence`, `sep_xi = TRUE`, separates the plotting of the
#'   influence curve for \eqn{\xi} from those of \eqn{\mu} and \eqn{\sigma}.
#'
#'   The example in **Examples** shows that extremely large block maxima have a
#'   strong positive influence on the MLE \eqn{\hat{\xi}} and also that
#'   extremely small block maxima have a strong negative influence on
#'   \eqn{\hat{\xi}},
#'
#' @return `gev_influence`: an object with class
#'   `c("gev_influence", "matrix", "array")`, a `length(z)` by `5` numeric
#'   matrix. The first two columns contain the input values in `z` and the
#'   corresponding values of `y`. Columns 3-5 contain the values of the GEV
#'   influence function for \eqn{\mu}, \eqn{\sigma} and \eqn{\xi} respectively
#'   at the values of `z`.
#'
#'   `plot.gev_influence`: nothing, only the plot is produced.
#' @references Hampel, F. R., Ronchetti, E. M., Rousseeuw, P. J., and
#'   Stahel, W. A. (2005). Robust Statistics. Wiley-Interscience, New York.
#'   \doi{10.1002/9781118186435}
#' @references Davison, A. C. and Smith, R. L. (1990). Models for exceedances
#'   over high thresholds. Journal of the Royal Statistical Society: Series B
#'   (Methodological), 52(3):393–425. \doi{10.1111/j.2517-6161.1990.tb01796.x}
#' @seealso [`gev_influence_rl`]
#' @examples
#' # Influence curve for the default mu = 0, sigma = 1, xi = 0 case
#' z <- seq(from = -3, to = 3, by = 0.01)
#' inf <- gev_influence(z = z)
#' plot(inf)
#'
#' # Influence curves based on the adjusted fit to the Plymouth ozone data
#' fit <- gev_mle(PlymouthOzoneMaxima)
#' pars <- coef(fit)
#' infp <- gev_influence(z = z, mu = pars[1], sigma = pars[2], xi = pars[3])
#' plot(infp)
#' @name gev_influence
NULL

#' @rdname gev_influence
#' @param x An object inheriting from class `"gev_influence"`, a result of a
#'   call to [`gev_influence`].
#' @param xvar A logical scalar.
#'   If `xvar = "z"` then the influence curves are plotted against the standard
#'   normal quantiles in `x[, "z"]`.
#'   If `xvar = "y"` then the influence curves are plotted against the
#'   corresponding GEV quantiles in `x[, "y"]`.
#' @param sep_xi A logical scalar. If `sep_xi = TRUE` then separate vertical
#'   scales are used for \eqn{\xi} and \eqn{(\mu, \sigma)}. The scale for
#'   \eqn{\xi} appears on the left and the scale for \eqn{(\mu, \sigma)} on the
#'   right.
#' @param vlines A numeric vector. If `vlines` is supplied then black dashed
#'   vertical lines are added to the plot at the values in `vlines` on the
#'   horizontal axis. This might be used to indicate the values of certain
#'   observations in a dataset.
#' @param ... For `plot.gev_influence`: to pass graphical parameters to the
#'   graphical functions [`matplot`][graphics::plot] and
#'   [`legend`][graphics::legend]. The parameters `col, lty` and `lwd` can be
#'    used to control line colour, type and width, with the parameters in the
#'    usual order, that is, \eqn{(\mu, \sigma, \xi)}.
#' @export
gev_influence <- function(z, mu = 0, sigma = 1, xi = 0) {
  # Check mu, sigma and xi
  if (any(length(mu) > 1, length(sigma) > 1, length(xi) > 1)) {
    stop("mu, sigma and xi must be scalar")
  }
  # Convert the normal quantile to the GEV(mu, sigma, xi) scale
  y <- nieve::qGEV(stats::pnorm(z), loc = mu, scale = sigma, shape = xi,
                   lower.tail = TRUE)
  # Calculate the GEV score at y
  loglik <- nieve::dGEV(y, loc = mu, scale = sigma, shape = xi, log = TRUE,
                        deriv = TRUE)
  score <- attr(loglik, "gradient")
  # Calculate the GEV expected information at (mu, sigma, xi)
  info <- gamlssx::gevExpInfo(scale = sigma, shape = xi)
  # Calculate the influence function value
  influence <- score %*% solve(info)
  influence <- cbind(z = z, y = y, influence)
  class(influence) <- c("gev_influence", class(influence))
  return(influence)
}

#' @rdname gev_influence
#' @export
plot.gev_influence <- function(x, xvar = c("z", "y"), sep_xi = TRUE, vlines,
                               ...) {
  # If cex.lab has been passed then extract it so that the label on the right
  # axis is of the correct size
  dots <- list(...)
  if (!is.null(dots$cex.lab)) {
    cex_lab <- dots$cex.lab
  } else {
    cex_lab <- 1
  }
  # Select the variable for the horizontal axis
  xvar <- match.arg(xvar)
  if (xvar == "z") {
    xvar <- x[, "z"]
    my_xlab <- "normal quantile"
  } else {
    xvar <- x[, "y"]
    my_xlab <- "GEV quantile"
  }
  y <- x[, -(1:2)]
  # Reset graphical parameters on exit
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  # Make room for the axis on the right
  if (!is.null(dots$main) && dots$main == "") {
    graphics::par(mar = c(4.5, 4.5, 1, 4.5) + 0.1)
  } else {
    graphics::par(mar = c(4.5, 4.5, 2.5, 4.5) + 0.1)
  }
  # Axis labels
  mu_sigma_label <- expression(paste("influence on  " * mu, "  and  ", sigma))
  xi_label <- expression("influence on  " * xi)
  # Create the main plotting function
  matplot_fn <- function(x, y, ..., type = "l", main = "GEV influence curves",
                         xlab = my_xlab, ylab = my_ylab,
                         lty = 1, lwd = 2, col = 1:3, side = 2, horiz) {
    if (side == 2) {
      plot_axes <- TRUE
    } else {
      plot_axes <- FALSE
    }
    graphics::matplot(x = x, y = y, axes = plot_axes, ..., type = type,
                      main = main, xlab = xlab, ylab = ylab,
                      lty = lty, lwd = lwd, col = col)
    if (side == 4) {
      graphics::axis(side = side, ...)
      graphics::mtext(mu_sigma_label, side = side, line = 3, cex = cex_lab)
    }
  }
  # Extract any colour vector in ...
  dots <- list(...)
  if (is.null(dots$col)) {
    col_vec <- 1:3
  } else {
    col_vec <- rep_len(dots$col, 3)
    dots$col <- NULL
  }
  # Similarly, extract any lty and lwed vectors in ...
  if (is.null(dots$lty)) {
    lty_vec <- rep(1, 3)
  } else {
    lty_vec <- rep_len(dots$lty, 3)
    dots$lty <- NULL
  }
  if (is.null(dots$lwd)) {
    lwd_vec <- rep(1, 3)
  } else {
    lwd_vec <- rep_len(dots$lwd, 3)
    dots$lwd <- NULL
  }
  # Plot the influence curves for xi (and mu and sigma if sep_xi is FALSE)
  if (sep_xi) {
    my_ylab <- xi_label
  } else {
    my_ylab <- "influence"
  }
  if (sep_xi) {
    matplot_args <- c(list(x = xvar, y = y[, 3]), dots,
                      list(col = col_vec[3], lty = lty_vec[3],
                           lwd = lwd_vec[3]))
  } else {
    matplot_args <- c(list(x = xvar, y = y[, 1:3]), dots,
                      list(col = col_vec, lty = lty_vec, lwd = lwd_vec))
  }
  do.call(matplot_fn, matplot_args)
  # If vlines is provided then add vertical lines to the plot
  if (!missing(vlines)) {
    graphics::abline(v = vlines, lty = 2)
  }
  if (sep_xi) {
    # Add the influence curves for mu and sigma with an axis on the right
    graphics::par(new = TRUE)
    matplot_args <- c(list(x = xvar, y = y[, 1:2]), dots,
                      list(col = col_vec[1:2], lty = lty_vec[1:2],
                           lwd = lwd_vec[1:2], side = 4))
    do.call(matplot_fn, matplot_args)
  }
  # Add the legend
  legend_text <- c(expression(mu), expression(sigma), expression(xi))
  legend_fn <- function(..., legend = legend_text, lty = 1, lwd = 2,
                        col = 1:3, main, xlab, ylab, xlim, ylim,
                        cex, cex.axis, cex.lab, cex.main, cex.sub,
                        bty = "n") {
    graphics::legend("top", ..., legend = legend, lty = lty_vec, lwd = lwd,
                     col = col, bty = bty)
  }
  legend_args <- c(dots, list(col = col_vec, lty = lty_vec, lwd = lwd_vec))
  do.call(legend_fn, legend_args)
  return(invisible())
}
