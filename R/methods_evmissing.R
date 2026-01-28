#' Methods for objects of class `"evmissing"`
#'
#' Methods for objects of class `"evmissing"` returned from [`gev_mle`].
#' @param object An object inheriting from class `"evmissing"`, a result of a
#'   call to [`gev_mle`].
#' @param x An object returned by `summary.evmissing`.
#' @param ... Further arguments. Only used in the following cases.
#'
#'   * `plot.evmissing`: to pass graphical parameters to the graphical functions
#'     [`plot`][graphics::plot], [`matplot`][graphics::matplot],
#'     [`abline`][graphics::abline], [`lines`][graphics::lines],
#'     [`matlines`][graphics::matlines] and [`points`][graphics::points].
#'     In particular, `col`, `lty` and `lwd` may be used to control the colour,
#'     type and width of lines and `pch` the type of plotting symbol. All data
#'     points are coloured black in all plots, which cannot be changed.
#'   * `print.summary.evmissing`: to pass arguments to [`print`].
#' @details The plots produced by `plot.evmissing` are of a similar form to the
#'   visual diagnostics is the `ismev` package and described in Coles (2001),
#'   that is, a probability plot (`which = "pp"` or `which = 1`), a quantile
#'   plot (`which ="qq"` of `which = 2`), a return level plot
#'   (`which = "return"` or `which = 3`) and a histogram of block maxima with a
#'   fitted GEV density superimposed (`which = "density"` or `which = 4`).
#'   Pointwise confidence bands of level `level` are added to the probability
#'   plot and quantile plot.
#'
#'   The default setting for confidence intervals for a return level plot
#'   produced by `plot.evmissing` is `profile = TRUE`, which uses [`gev_return`]
#'   and [`confint.return_level`]. The plot takes longer to produce, but it
#'   avoids the unrealistic feature of the lower confidence limits decreasing
#'   as we extrapolate to long return periods.
#'
#'   If `adjust = TRUE` then the empirical values based on the observed block
#'   maxima are adjusted for the number of non-missing raw observations in
#'   each block based on the fitted GEV parameter values for reduced block
#'   sizes. Passing `adjust = FALSE` is not sensible, but, if there are missing
#'   data, then it can serve to show that making the adjustment is necessary to
#'   give the correct impression of how well the model has fitted the data.
#'
#' @references Coles, S. G. (2001) *An Introduction to Statistical
#'   Modeling of Extreme Values*, Springer-Verlag, London.
#'   \doi{10.1007/978-1-4471-3675-0_3}
#' @references Heffernan, J. E. and Stephenson, A. G. (2018). *ismev: An
#'   Introduction to Statistical Modeling of Extreme Values*. R package version
#'   1.42. \doi{10.32614/CRAN.package.ismev}
#' @return `coef.evmissing`: a numeric vector of length 3 with names
#'   `c("mu", "sigma", "xi")`.  The MLEs of the parameters \eqn{\mu},
#'   \eqn{\sigma} and \eqn{\xi}.
#'
#'   `vcov.evmissing`: a \eqn{3 \times 3}{3 x 3} matrix with row and column
#'   names `c("mu", "sigma", "xi")`. The estimated variance-covariance matrix
#'   for the model parameters \eqn{\mu}, \eqn{\sigma} and \eqn{\xi}.
#'
#'   `nobs.evmissing`: a numeric scalar. The number of maxima used in the model
#'   fit.
#'
#'   `logLik.evmissing`: an object of class \code{"logLik"}: a numeric scalar
#'   with value equal to the maximised log-likelihood. The returned object
#'   also has attributes `nobs`, the number of maxima used in the model fit
#'   and `"df"` (degrees of freedom), which is equal to the number of total
#'   number of parameters estimated (3).
#'
#'   `summary.evmissing`: an object with class `"summary.evmissing"` containing
#'   the original function call and a matrix of estimates and estimated
#'   standard errors with row names `c("mu", "sigma", "xi")`.  The object is
#'   printed by [`print.summary.evmissing`].
#'
#'   `print.summary.evmissing`: the argument `x` is returned, invisibly.
#'
#'   `confint.evmissing`: an object of class `c("confint_gev", "evmissing")`.
#'   A numeric matrix with 2 columns giving the lower and upper confidence
#'   limits for each parameter. These columns are labelled as `(1-level)/2` and
#'   `1-(1-level)/2`, expressed as a percentage, by default `2.5%` and `97.5%`.
#'   The row names are the names of the parameters supplied in `parm`.
#'   The ordering `"mu"`, `"sigma"`, `"xi"` is retained regardless of the
#'   ordering of the parameters in `parm`.
#'   If `profile = TRUE` then the returned object has extra attributes `crit`,
#'   `level` and `for_plot`. The latter is a named list of length 3 with
#'   components `mu`, `sigma` and `xi`. Each components is a 2-column numeric
#'   matrix. The first column (named `mu_values` etc) contains values of the
#'   parameter and the second column the corresponding values of the profile
#'   log-likelihood. The profile log-likelihood is equal to the attribute
#'   `crit` at the limits of the confidence interval. The attribute `level` is
#'   the input argument `level`.
#'
#'   `plot.evmissing`: if a return level plot has been requested, a 3-column
#'   matrix containing the values plotted in the return level plot. Column 2
#'   contains the estimated return levels and columns 1 and 3 the lower and
#'   upper confidence limits.
#' @examples
#' ## Plymouth ozone data
#'
#' # Make adjustment for the numbers of non-missing values per block
#' fit <- gev_mle(PlymouthOzoneMaxima)
#' coef(fit)
#' vcov(fit)
#' nobs(fit)
#' logLik(fit)
#' summary(fit)
#'
#' ## Model diagnostic plots
#'
#' # When profile = FALSE the return confidence limits are unrealistic
#' # for long return periods
#' plot(fit, profile = FALSE)
#'
#' # Create the return level plot only
#' # When profile = TRUE (the default) the confidence limits are fine
#' # but the plot takes longer
#' # For speed, we reduce the number, num, of points used to plot the curves
#' plot(fit, which = 3, num = 8)
#'
#' # If we do not reflect the adjustment in the plot then it gives a false
#' # impression of how well the model has fitted the data
#' plot(fit, adjust = FALSE, profile = FALSE)
#'
#' ## Confidence intervals
#'
#' # Confidence limits that are symmetric about the respective MLEs
#' confint(fit)
#'
#' # Calling confint to produce a smooth profile log-likelihood plot
#' x <- confint(fit, profile = TRUE)
#' x
#' plot(x, parm = "xi")
#'
#' # Doing this more quickly when we only want the approximate confidence limits
#' x <- confint(fit, profile = TRUE, mult = 32, faster = TRUE)
#' x
#' plot(x, parm = "xi", type = "b")
#' @seealso [`gev_mle`] and [`confint_gev_methods`].
#' @name evmissing_methods
NULL
## NULL

# ================================ coef.evmissing =========================== #

#' Extract model coefficients method for objects of class `"evmissing"`
#'
#' @rdname evmissing_methods
#' @export
coef.evmissing <- function(object, ...) {
  return(object$par)
}

# ================================ vcov.evmissing =========================== #

#' Calculate the variance-covariance matrix for an object of class `"evmissing"`
#'
#' @rdname evmissing_methods
#' @export
vcov.evmissing <- function(object, ...) {
  vc <- object$vcov
  vc_names <- names(coef(object))
  dim(vc) <- c(length(vc_names), length(vc_names))
  dimnames(vc) <- list(vc_names, vc_names)
  return(vc)
}

# ================================ nobs.evmissing =========================== #

#' Extract the number of observations for objects of class `"evmissing"`
#'
#' @rdname evmissing_methods
#' @export
nobs.evmissing <- function(object, ...) {
  return(length(object$maxima))
}

# ================================ logLik.evmissing ========================= #

#' Extract log-likelihood for objects of class `"evmissing"`
#'
#' @rdname evmissing_methods
#' @export
logLik.evmissing <- function(object, ...) {
  val <- object$loglik
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}

# =============================== summary.evmissing ========================= #

#' Summarising GEV fits
#'
#' @param digits An integer. Passed to [`signif`] to round the
#'   values in the summary.
#' @rdname evmissing_methods
#' @export
summary.evmissing <- function(object,
                              digits = max(3, getOption("digits") - 3L), ...) {
  res <- list()
  res$call <- attr(object, "call")
  mles <- signif(coef(object), digits = digits)
  ses <- signif(sqrt(diag(vcov(object))), digits = digits)
  res$matrix <- cbind(`Estimate` = mles, `Std. Error` = ses)
  rownames(res$matrix) <- names(mles)
  class(res) <- "summary.evmissing"
  return(res)
}

# ============================ print.summary.evmissing ====================== #

#' Print method for objects of class `"summary.evmissing"`
#'
#' @rdname evmissing_methods
#' @export
print.summary.evmissing <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  print(x$matrix, ...)
  return(invisible(x))
}

# ================================ confint.evmissing ======================== #

#' Confidence intervals for \code{"evmissing"} objects
#'
#' @param parm A character vector specifying the parameters for which
#'   confidence intervals are to be calculated. The default, `which = "all"`,
#'   produces confidence intervals for all the parameters, that is,
#'   \eqn{\mu}, \eqn{\sigma} and \eqn{\xi}. Otherwise, `parm` must be a subset
#'   of `c("mu", "sigma", "xi")`.
#' @param level The confidence level required.  A numeric scalar in (0, 1).
#' @param profile A logical scalar. If `TRUE` then confidence intervals
#'   based on a profile log-likelihood are returned.  If `FALSE` then intervals
#'   based on approximate large sample normal theory, which are symmetric about
#'   the MLE, are returned.
#' @param mult A positive numeric scalar. Controls the increment by which the
#'   parameter of interest is increased/decreased when profiling above/below
#'   its MLE. The increment is `mult * se / 100` where `se` is the estimated
#'   standard error of the estimator of the parameter. Decreasing `mult`
#'   profiles at more points but will be slower. The default, `mult = 2` should
#'   be sufficiently small to produce a smooth looking plot of the profile
#'   log-likelihood using [`plot.confint_gev`]. To estimate the confidence
#'   limits more quickly, the value of `mult` can be increased and/or the
#'   argument `faster` set to `TRUE`.
#' @param faster A logical scalar. If `faster = TRUE` then the profiling of the
#'   log-likelihood in search of a lower (upper) confidence limit is started at
#'   the corresponding symmetric lower (upper) confidence limit.
#' @param epsilon Only relevant if `profile = TRUE`. A numeric vector of values
#'   that determine the accuracy of the confidence limits. `epsilon` is
#'   recycled to the length of the parameter vector `parm`.
#'
#'   * If `epsilon[i] > 0` then this value is passed as the argument `epsilon`
#'     to the [`itp::itp`] function, which estimates the parameter values for
#'     which the profile log-likelihood for parameter `i` drops to the value
#'     that defines the confidence limits, once profiling has been successful
#'     in finding an interval within which this value lies.
#'
#'    * If `epsilon[i] < 0` monotonic cubic spline interpolation is used, which
#'      will tend to be faster.
#'
#'   * If `epsilon[i] = 0` then linear interpolation is used, which will be
#'     faster still.
#' @details For `confint.evmissing`, the default, `epsilon = -1`, should work well
#'   enough in most circumstances, but to achieve a specific accuracy set
#'   `epsilon` to be a small positive value, for example, `epsilon = 1e-4`.
#' @rdname evmissing_methods
#' @export
confint.evmissing <- function(object, parm = "all", level = 0.95,
                              profile = FALSE, mult = 2, faster = FALSE,
                              epsilon = 1e-4, ...) {
  # GEV parameter names
  parm_values <- c("mu", "sigma", "xi")
  # Convert numeric values for parm to parameter names
  if (is.numeric(parm)) {
    if (!all(is.element(parm, 1:3))) {
      stop("''parm'' must be a subset of {", paste(1:3, collapse = ","), "}")
    } else {
      parm <- parm_values[parm]
    }
  }
  # Remove any duplicates in parm
  parm <- unique(parm)
  # Check inputs
  check_values <- c("mu", "sigma", "xi", "all")
  p_message <- "c(''mu'', ''sigma'', ''xi'')"
  if (!all(is.element(parm, check_values))) {
    stop(paste("''parm'' must be ''all'', or a subset of", p_message))
  }
  # We use order() below to retain the mu, sigma, xi ordering
  if (length(parm) == 1 && parm == "all") {
    parm <- parm_values
  } else {
    parm <- parm[order(parm)]
  }
  # Logical vector indicating which parameters to include
  which_parm <- is.element(parm_values, parm)
  # Check the input confidence level
  if (level <= 0 | level >= 1) {
    stop("''level'' must be in (0, 1)")
  }

  # Calculate symmetric intervals

  z_val <- stats::qnorm(1 - (1 - level) / 2)
  mles <- coef(object)[which_parm]
  ses <- sqrt(diag(vcov(object)))[which_parm]
  sym_lower <- mles - z_val * ses
  sym_upper <- mles + z_val * ses
  ci_mat <- cbind(sym_lower, sym_upper)
  rownames(ci_mat) <- parm
  # Constrain any limits for sigma to [0, Infty]
  if (is.element("sigma", rownames(ci_mat))) {
    sigma_row <- which(rownames(ci_mat) == "sigma")
    ci_mat[sigma_row, 1] <- max(0, ci_mat[sigma_row, 1])
  }
  ci_sym_mat <- ci_mat

  # If profile log-likelihood-based intervals are required then calculate them
  # Was the fitted object produced by gev_weighted()?
  if (inherits(object, "weighted_mle")) {
    weighted_fit <- TRUE
  } else {
    weighted_fit <- FALSE
  }
  if (profile) {
    # The number of parameters
    n_parm <- length(parm)
    # Force epsilon to have length equal to the number of parameters
    epsilon <- rep_len(epsilon, n_parm)
    # Extract the parameter numbers to include
    parm_numbers <- (1:length(parm_values))[which_parm]
    # Recreate the list maxima_notNA
    maxima_notNA <- list(maxima = object$maxima, notNA = object$notNA,
                         n = object$n)
    # An empty list in which to store the profile log-likelihood values
    # Likewise for the values of the GEV parameters at the confidence limits
    for_plot <- list()
    lower_pars <- list()
    upper_pars <- list()
    # Loop over all parameters
    # In a small number of cases, particularly those where block maxima have
    # been discarded, there can be problems evaluating the log-likelihood at
    # the initial values used in the profiling, that is, the most recent
    # parameters are not admissible at the new value for the parameter being
    # profiled. This tends to happen when the MLE for xi is much less than 0
    # and the log-likelihood is relatively flat owing to a small sample size.
    # To avoid these problems, we reduce (halve) mult iteratively down to the
    # minimum value min_mult. If no fit is successful then we return NAs.
    min_mult <- 1
    if (mult < min_mult) {
      min_mult <- mult
    }
    save_mult <- mult
    for (i in 1:n_parm) {
      success <- FALSE
      if (faster) {
        while (mult >= min_mult) {
          # Set inc based on the estimated standard errors
          inc <- mult * ses / 100
          if (weighted_fit) {
            conf_list <- faster_profile_ci(negated_loglik_fn =
                                             weighted_negated_gev_loglik,
                                           which = parm_numbers[i], level = level,
                                           mle = coef(object),
                                           ci_sym_mat = ci_sym_mat,
                                           inc = inc[i],
                                           epsilon = epsilon[i],
                                           maxima = object$maxima,
                                           weights = object$weights)
          } else {
            conf_list <- faster_profile_ci(negated_loglik_fn = negated_gev_loglik,
                                           which = parm_numbers[i], level = level,
                                           mle = coef(object),
                                           ci_sym_mat = ci_sym_mat,
                                           inc = inc[i],
                                           epsilon = epsilon[i],
                                           maxima_notNA = maxima_notNA,
                                           adjust = object$adjust)
          }
          if (!is.null(conf_list$optim_error)) {
            mult <- mult / 2
          } else {
            success <- TRUE
            mult <- min_mult - 1
          }
        }
      } else {
        while (mult >= min_mult) {
          # Set inc based on the estimated standard errors
          inc <- mult * ses / 100
          if (weighted_fit) {
            conf_list <- profile_ci(negated_loglik_fn =
                                      weighted_negated_gev_loglik,
                                    which = parm_numbers[i], level = level,
                                    mle = coef(object), inc = inc[i],
                                    epsilon = epsilon[i],
                                    maxima = object$maxima,
                                    weights = object$weights)
          } else {
            conf_list <- profile_ci(negated_loglik_fn = negated_gev_loglik,
                                    which = parm_numbers[i], level = level,
                                    mle = coef(object), inc = inc[i],
                                    epsilon = epsilon[i],
                                    maxima_notNA = maxima_notNA,
                                    adjust = object$adjust)
          }
          if (!is.null(conf_list$optim_error)) {
            mult <- mult / 2
          } else {
            success <- TRUE
            mult <- min_mult - 1
          }
        }
      }
      if (success) {
        ci_mat[i, ] <- conf_list$par_prof[c(1, 3)]
        colnames(conf_list$for_plot)[1] <- paste0(parm[i], "_values")
        for_plot[[i]] <- conf_list$for_plot
        lower_pars[[i]] <- conf_list$lower_pars
        upper_pars[[i]] <- conf_list$upper_pars
      } else {
        ci_mat[i, ] <- matrix(NA, 1, 2)
        for_plot[[i]] <- NA
        lower_pars[[i]] <- NA
        upper_pars[[i]] <- NA
      }
      # Reset mult
      mult <- save_mult
    }
    names(for_plot) <- parm
  }
  # Format the matrix to be returned
  low <- paste0(100 * (1 - level)/ 2, "%")
  up <- paste0(100 - 100 * (1 - level)/ 2, "%")
  colnames(ci_mat) <- c(low, up)
  # If profile = TRUE save the profile log-likelihood values for plotting
  if (profile) {
    attr(ci_mat, "for_plot") <- for_plot
    attr(ci_mat, "crit") <- conf_list$crit
    attr(ci_mat, "level") <- level
    if (any(epsilon > 0) || !faster) {
      names(lower_pars) <- c("when profiling for mu", "when profiling for sigma",
                             "when profiling for xi")[which_parm]
      names(upper_pars) <- c("when profiling for mu", "when profiling for sigma",
                             "when profiling for xi")[which_parm]
    }
    attr(ci_mat, "lower_pars") <- lower_pars
    attr(ci_mat, "upper_pars") <- upper_pars
  }
  class(ci_mat) <- c("confint_gev", "evmissing")
  return(ci_mat)
}

#' @param adjust If `adjust = TRUE` then the diagnostic plots produced by
#'   `plot.evmissing` are adjusted for the number of non-missing observations
#'   contributing to each block maximum. Otherwise, no adjustment is made.
#' @param which If supplied, this must either be a character scalar, one of
#'   `"pp"`, `"qq"`, `"return"` or `"density"` or a numeric scalar in `1:4`,
#'   with `1` corresponding to `"pp"` etc. If `which` is missing then all four
#'   plots are produced in a 2 by 2 display.
#' @param m A numeric vector of return periods to label on the
#'   horizontal axis of the **return level plot**. Along with the data, the
#'   smallest and largest return period values in `m` influence the range of
#'   return periods for which return level estimates are plotted. All values in
#'   `m` must be greater than 1.
#' @param num An integer scalar. The number of return level estimates to
#'   calculate to produce the return level curve and pointwise confidence
#'   limits in the **return level plot**. The default setting is approximately
#'   5 times `log(max(m), base = 10)`. If `profile = TRUE` then reducing `num`
#'   will speed up the calculation of the confidence limits, at the expense of
#'   a reduction in smoothness of the curves.
#' @param npy A numeric scalar.  The number \eqn{n_{py}} of block maxima per
#'   year. If the blocks are of length 1 year then `npy = 1`. This is only
#'   used in the **return level plot**.
#' @rdname evmissing_methods
#' @export
plot.evmissing <- function(x, adjust = TRUE,
                           which = c("pp", "qq", "return", "density"),
                           m = c(2, 10, 100, 1000), level = 0.95,
                           profile = TRUE, num, npy = 1, ...) {
  # Choose which plots to produce
  if (!missing(which)) {
    if (is.character(which)) {
      which <- match.arg(which)
    } else if (is.numeric(which) && length(which) == 1 &&
               is.element(which, 1:4)) {
      which <- switch(which,
                      "1" = "pp",
                      "2" = "qq",
                      "3" = "return",
                      "4" = "density")
    } else {
      w_char <- paste(c("\"pp\"", "\"qq\"", "\"return\"", "\"density\""),
                      collapse =", ")
      w_num <- paste(1:4, collapse = ", ")
      stop("'which' should be one of: ", w_char, " or ", w_num)
    }
  }
  # Reset graphical parameters on exit
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  # Set the display according to the number of plots (4 or 1)
  if (length(which) == 4) {
    graphics::par(mfrow = c(2, 2))
  } else {
    graphics::par(mfrow = c(1, 1))
  }
  # Force the plots to be square and adjust the margins
  graphics::par(pty = "s", mar = c(4.5, 1, 3, 0))

  # Create the plots

  # Probability plot
  if (is.element("pp", which)) {
    gev_pp(x = x, adjust = adjust, level = level, ...)
  }

  # Quantile plot
  if (is.element("qq", which)) {
    gev_qq(x = x, adjust = adjust, level = level, ...)
  }

  # Return level plot
  if (is.element("return", which)) {
    if (any(m <= 1)) {
      stop("All return period values in 'm' must be > 1")
    }
    # If num has not been supplied then set it
    if (missing(num)) {
      max_r <- max(max(m), nobs(x) + 1)
      num <- ceiling(5 * log(max_r, base = 10))
    } else {
      num <- ceiling(num)
    }
    rl_cis <- gev_rl(x = x, adjust = adjust, m = m, level = level,
                     profile = profile, num = num, npy = npy, ...)
    val <- rl_cis
  } else {
    val <- NULL
  }

  # Density plot
  if (is.element("density", which)) {
    gev_his(x = x, adjust = adjust, ...)
  }

  return(invisible(val))
}
