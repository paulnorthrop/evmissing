#' GEV influence curves for return levels
#'
#' Calculates influence function curves for maximum likelihood estimators of
#' 3 return levels based on Generalised Extreme Value (GEV) parameters.
#'
#' @param z A numeric vector. Values of normal quantiles \eqn{z} at which to
#'  calculate the GEV influence function. See **Details**.
#' @param mu,sigma,xi Numeric scalars supplying the values of the GEV
#'   parameters \eqn{\mu}, \eqn{\sigma} and \eqn{\xi}.
#' @param m A numeric vector of length 3 containing 3 unique return periods
#'   in years. All entries in `m` must be greater than 1.
#' @param npy A numeric scalar.  The number \eqn{n_{py}} of block maxima per
#'   year. If the blocks are of length 1 year then `npy = 1`.
#' @details See [`gev_influence`] for information about influence functions in
#'   general and influence curves for the parameters of a GEV distribution in
#'   particular. The GEV influence curves are reparameterised from
#'   \eqn{(\mu, \sigma, \xi)} to the required return levels.
#' @return `gev_influence_rl`: an object with class
#'   `c("gev_influence_rl", "matrix", "array")`, a `length(z)` by `5` numeric
#'   matrix. The first two columns contain the input values in `z` and the
#'   corresponding values of `y`. Columns 3-5 contain the values of the GEV
#'   influence function for the return levels in `m` respectively at the values
#'   of `z`.
#'
#'   `plot.gev_influence_rl`: nothing, only the plot is produced.
#' @references Hampel, F. R., Ronchetti, E. M., Rousseeuw, P. J., and
#'   Stahel, W. A. (2005). Robust Statistics. Wiley-Interscience, New York.
#'   \doi{10.1002/9781118186435}
#' @references Davison, A. C. and Smith, R. L. (1990). Models for exceedances
#'   over high thresholds. Journal of the Royal Statistical Society: Series B
#'   (Methodological), 52(3):393–425. \doi{10.1111/j.2517-6161.1990.tb01796.x}
#' @seealso [`gev_influence`], [`gev_return`]
#' @examples
#' # Influence curves based on the adjusted fit to the Plymouth ozone data
#' z <- seq(from = -3, to = 3, by = 0.01)
#' fit <- gev_mle(PlymouthOzoneMaxima)
#' pars <- coef(fit)
#' m <- c(25, 50, 100)
#' infp <- gev_influence_rl(z = z, mu = pars[1], sigma = pars[2], xi = pars[3],
#'                          m = m)
#' plot(infp)
#' @name gev_influence_rl
NULL

#' @rdname gev_influence_rl
#' @param x An object inheriting from class `"gev_influence_rl"`, a result of a
#'   call to [`gev_influence_rl`].
#' @param xvar A logical scalar.
#'   If `xvar = "z"` then the influence curves are plotted against the standard
#'   normal quantiles in `x[, "z"]`.
#'   If `xvar = "y"` then the influence curves are plotted against the
#'   corresponding GEV quantiles in `x[, "y"]`.
#' @param vlines A numeric vector. If `vlines` is supplied then black dashed
#'   vertical lines are added to the plot at the values in `vlines` on the
#'   horizontal axis. This might be used to indicate the values of certain
#'   observations in a dataset.
#' @param ... For `plot.gev_influence_rl`: to pass graphical parameters to the
#'   graphical functions [`matplot`][graphics::plot] and
#'   [`legend`][graphics::legend]. The parameters `col, lty` and `lwd` can be
#'   used to control line colour, type and width, with the return levels in
#'   the order that they were supplied in `m`.
#' @export
gev_influence_rl <- function(z, mu = 0, sigma = 1, xi = 0, m, npy = 1) {
  # Check mu, sigma and xi
  if (any(length(mu) > 1, length(sigma) > 1, length(xi) > 1)) {
    stop("mu, sigma and xi must be scalar")
  }
  # Check that m contains at least 3 unique values
  if (length(unique(m)) != 3 || any(m <= 1)) {
    stop("m must be a numeric vector or 3 unique values, all > 1.")
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
  # Set the annual probability of exceedance based on m and npy
  p <- 1 - (1 - 1 / m) ^ (1 / npy)
  yp <- -log(1 - p)
  # Calculate the GEV score at y
  temp <- nieve::qGEV(p, loc = mu, scale = sigma, shape = xi,
                      lower.tail = FALSE, deriv = TRUE)
  delta <- t(attr(temp, "gradient"))
  # Transform the VC matrix to the return level scale
  rl_var <- t(delta) %*% solve(info) %*% delta
  # Invert to obtain the information matrix for the return levels
  # rl_info <- solve(rl_var)
  # In fact, we need the inverse of this below, rl_var
  #
  # Infer the score for the return levels
  rl_score <- t(solve(delta) %*% t(score))
  # Calculate the influence function value
  rl_influence <- rl_score %*% rl_var
  rl_influence <- cbind(z = z, y = y, rl_influence)
  rls <- paste0(m, "-year return level")
  colnames(rl_influence) <- c("z", "y", rls)
  # Save the return levels as an attribute
  attr(rl_influence, "m") <- m
  class(rl_influence) <- c("gev_influence_rl", class(rl_influence))
  return(rl_influence)
}

#' @rdname gev_influence_rl
#' @export
plot.gev_influence_rl <- function(x, xvar = c("z", "y"), vlines, ...) {
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
  left_label <- "influence on return level"
  xi_label <- expression("influence on  " * xi)
  # Create the main plotting function
  matplot_fn <- function(x, y, ..., type = "l",
                         main = "GEV return level influence curves",
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
      graphics::mtext(left_label, side = side, line = 3, cex = cex_lab)
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
  # Plot the influence curves
  my_ylab <- "influence"
  matplot_args <- c(list(x = xvar, y = y[, 1:3]), dots,
                    list(col = col_vec, lty = lty_vec, lwd = lwd_vec))
  do.call(matplot_fn, matplot_args)
  # If vlines is provided then add vertical lines to the plot
  if (!missing(vlines)) {
    graphics::abline(v = vlines, lty = 2)
  }
  # Add the legend
  legend_text <- as.character(attr(x, "m"))
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
