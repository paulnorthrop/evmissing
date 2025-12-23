#' Methods for objects of class `"return_level"`
#'
#' Methods for objects of class `"return_level"` returned from [`gev_return`].
#' @param object,x An object inheriting from class `"return_level"`, a result
#'   of a call to [`gev_return`]. `object` is a named numeric vector of MLEs of
#'   return levels.
#'
#'   For `print.summary.return_level`, this is an object returned by the
#'   function `summary.return_level`.
#' @param ... Further arguments. Only used for `print.summary.return_level` to
#'   pass arguments to [`print`].
#' @return `print.return_level` and `coef.return_level`: a numeric vector
#'   containing the MLEs of return return levels.
#'
#'   `vcov.return_level`: a `length(object)` by `length(object)` matrix with
#'   row and column names indicating the return periods of the return levels.
#'   The estimated variance-covariance matrix for the return levels in
#'   `object`. The diagonal elements give the estimated variances associated
#'   with the individual return level estimates.
#'
#'   `summary.return_level`: an object containing the original function call
#'     and a matrix of estimates of return levels and associated estimated
#'     standard errors with row names indicating the respective return periods.
#'     The object is printed by [`print.summary.return_level`].
#'
#'   `print.summary.return_level`: the argument `x` is returned, invisibly.
#'
#'   `confint.return_level`: an object of class
#'   `c("confint_return_level", "evmissing")`. A numeric matrix with 2 columns
#'   giving the lower and upper confidence limits for each return level. These
#'   columns are labelled as `(1-level)/2` and `1-(1-level)/2`, expressed as a
#'   percentage, by default `2.5%` and `97.5%`. The row names indicate the
#'   return levels.
#'   If `profile = TRUE` then the returned object has extra attributes `crit`,
#'   `level` and `for_plot`. The latter is a named list of length  `parm` with
#'   components named after the return periods. Each components is a 2-column
#'   numeric matrix. The first column contains values of the return level and
#'   the second column the corresponding values of the profile log-likelihood.
#'   The profile log-likelihood is equal to the attribute `crit` at the limits
#'   of the confidence interval. The attribute `level` is the input argument
#'   `level`.
#' @examples
#' ## Plymouth ozone data
#'
#' # See ?gev_return for confidence intervals for return levels
#' fit <- gev_mle(PlymouthOzoneMaxima)
#' rl <- gev_return(fit, m = c(100, 200))
#' rl
#' vcov(rl)
#' summary(rl)
#' @seealso [`gev_mle`] and [`gev_return`], for examples of the use of
#'   [`confint.return_level`].
#' @name return_level_methods
NULL
## NULL

# ============================ coef.return_level ============================ #

#' Extract model coefficients method for objects of class `"return_level"`
#'
#' @rdname return_level_methods
#' @export
coef.return_level <- function(object, ...) {
  return(object[1:length(object)])
}

# ============================ print.return_level ============================ #

#' Print method for objects of class `"return_level"`
#'
#' @rdname return_level_methods
#' @export
print.return_level <- function(x, ...) {
  print(x[1:length(x)], ...)
  return(invisible(x))
}

# ============================= vcov.return_level =========================== #

#' Calculate the variance-covariance matrix associated with an estimated return
#' level for an object of class `"return_level"`
#'
#' @rdname return_level_methods
#' @export
vcov.return_level <- function(object, ...) {
  # Extract m, npy and the fitted GEV model object
  m <- attr(object, "m")
  npy <- attr(object, "npy")
  gev_mle_object <- attr(object, "gev_mle")
  # Extract the MLEs and their estimated VC matrix
  mle <- coef(gev_mle_object)
  mu <- mle[1]
  sigma <- mle[2]
  xi <- mle[3]
  gev_vc <- vcov(gev_mle_object)
  # Calculate the estimated VC matrix for all return levels in object
  p <- 1 - (1 - 1 / m) ^ (1 / npy)
  yp <- -log(1 - p)
  temp <- nieve::qGEV(p, loc = mu, scale = sigma, shape = xi,
                      lower.tail = FALSE, deriv = TRUE)
  delta <- t(attr(temp, "gradient"))
  rl_var <- t(delta) %*% gev_vc %*% delta
  vc_names <- names(coef(object))
  dimnames(rl_var) <- list(vc_names, vc_names)
  return(rl_var)
}

# =========================== summary.return_level ========================== #

#' Summarising GEV fits
#'
#' @param digits An integer. Passed to [`signif`] to round the
#'   values in the summary.
#' @rdname return_level_methods
#' @export
summary.return_level <- function(object, digits = max(3, getOption("digits") - 3L),
                           ...) {
  res <- list()
  res$call <- attr(object, "call")
  mles <- signif(coef(object), digits = digits)
  ses <- signif(sqrt(diag(vcov(object))), digits = digits)
  res$matrix <- cbind(`Estimate` = mles, `Std. Error` = ses)
  rownames(res$matrix) <- names(mles)
  class(res) <- "summary.return_level"
  return(res)
}

# ======================== print.summary.return_level ======================= #

#' Print method for objects of class `"summary.return_level"`
#'
#' @rdname return_level_methods
#' @export
print.summary.return_level <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  print(x$matrix, ...)
  return(invisible(x))
}

# ============================ confint.return_level ========================= #

#' Confidence intervals for \code{"return_level"} objects
#'
#' @param parm A numeric vector. For which components, that is, which return
#'   levels, in `object` we require a confidence interval.
#' @param level The confidence level required.  A numeric scalar in (0, 1).
#' @param profile A logical scalar. If `TRUE` then confidence intervals
#'   based on a profile log-likelihood are returned.  If `FALSE` then intervals
#'   based on approximate large sample normal theory, which are symmetric about
#'   the MLE, are returned.
#' @param mult A positive numeric scalar. Controls the increment by which the
#'   parameter of interest is increased/decreased when profiling above/below
#'   its MLE. The increment is `mult * se / 100` where `se` is the estimated
#'   standard error of the estimator of the return level. Decreasing `mult`
#'   profiles at more points but will be slower.
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
#'    * If `epsilon[i] < 0` quadratic interpolation is used, which will tend to
#'      be faster.
#'
#'   * If `epsilon[i] = 0` then linear interpolation is used, which will be
#'     faster still.
#' @details For `confint.return_level`, the default, `epsilon = -1`, should
#'   work well enough in most circumstances, but to achieve a specific accuracy
#'   set `epsilon` to be a small positive value, for example, `epsilon = 1e-4`.
#' @rdname return_level_methods
#' @export
confint.return_level <- function(object, parm = 1:length(object), level = 0.95,
                                 profile = FALSE, mult = 2, faster = FALSE,
                                 epsilon = 1e-4, ...) {
  # Check parm
  p_message <- 1:length(object)
  parm_values <- 1:length(object)
  if (!all(is.element(parm, parm_values))) {
    stop("''parm'' must be a subset of {", paste(parm_values, collapse = ","),
         "}")
  }
  parm <- parm[order(parm)]
  # Logical vector indicating which parameters to include
  which_parm <- is.element(parm_values, parm)
  # Extract the parameter numbers to include
  parm_numbers <- (1:length(parm_values))[which_parm]
  # Check the input confidence level
  if (level <= 0 | level >= 1) {
    stop("''level'' must be in (0, 1)")
  }

  # Calculate symmetric intervals

  z_val <- stats::qnorm(1 - (1 - level) / 2)
  mles <- coef(object)
  ses <- sqrt(diag(vcov(object)))
  sym_lower <- mles - z_val * ses
  sym_upper <- mles + z_val * ses
  ci_mat <- cbind(sym_lower, sym_upper)
  rownames(ci_mat) <- names(coef(object))
  ci_sym_mat <- ci_mat

  # If profile log-likelihood-based intervals are required then calculate them
  # Was the fitted object produced by gev_weighted()?
  if (inherits(attr(object, "gev_mle"), "weighted_mle")) {
    weighted_fit <- TRUE
  } else {
    weighted_fit <- FALSE
  }

  if (profile) {
    # Extract the fitted GEV object returned by gev_mle()
    gev_object <- attr(object, "gev_mle")
    # The number of parameters
    n_parm <- length(parm)
    # Force epsilon to have length equal to the total number of return levels
    epsilon <- rep_len(epsilon, length(object))
    # Extract from gev_object the required quantities
    if (weighted_fit) {
      maxima <- gev_object$maxima
      weights <- gev_object$weights
    } else {
      maxima_notNA <- list(maxima = gev_object$maxima, notNA = gev_object$notNA,
                           n = gev_object$n)
      adjust <- gev_object$adjust
    }
    # An empty list in which to store the profile log-likelihood values
    for_plot <- list()
    # Set inc based on the estimated standard errors for the return levels
    ses <- sqrt(diag(vcov(object)))
    inc <- mult * ses / 100
    # Extract the MLEs for the GEV parameters
    mle <- coef(gev_object)
    # Extract the MLEs of the return levels
    return_level_mle <- coef(object)
    # Extract the values of m and npy
    m <- attr(object, "m")
    npy <- attr(object, "npy")
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
    for (i in parm_numbers) {
      success <- FALSE
      mle_to_pass <- c(return_level_mle[i], mle[2:3])
      if (faster) {
        # Force epsilon to be positive
        while (mult >= min_mult) {
          ci_init <- return_level_profile_init(gev_object = gev_object,
                                               level = level, mult = mult,
                                               epsilon = 1e-4, m = m[i],
                                               npy = npy)
          # Reset inc
          inc <- mult * ses / 100
          if (weighted_fit) {
            conf_list <- faster_profile_ci(negated_loglik_fn =
                                             weighted_negated_gev_loglik_ret_levs,
                                           which = 1, level = level,
                                           mle = mle_to_pass,
                                           ci_sym_mat =
                                             ci_sym_mat[i, , drop = FALSE],
                                           inc = inc[i],
                                           epsilon = epsilon[i],
                                           maxima = maxima,
                                           weights = weights, m = m[i],
                                           npy = npy, ci_init = ci_init)
          } else {
            conf_list <- faster_profile_ci(negated_loglik_fn =
                                             negated_gev_loglik_ret_levs,
                                           which = 1, level = level,
                                           mle = mle_to_pass,
                                           ci_sym_mat =
                                             ci_sym_mat[i, , drop = FALSE],
                                           inc = inc[i],
                                           epsilon = epsilon[i],
                                           maxima_notNA = maxima_notNA,
                                           adjust = adjust, m = m[i],
                                           npy = npy, ci_init = ci_init)
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
          # Reset inc
          inc <- mult * ses / 100
          if (weighted_fit) {
            conf_list <- profile_ci(negated_loglik_fn =
                                      weighted_negated_gev_loglik_ret_levs,
                                    which = 1, level = level,
                                    mle = mle_to_pass, inc = inc[i],
                                    epsilon = epsilon[i],
                                    maxima = maxima,
                                    weights = weights, m = m[i],
                                    npy = npy)
          } else {
            conf_list <- profile_ci(negated_loglik_fn =
                                      negated_gev_loglik_ret_levs,
                                    which = 1, level = level,
                                    mle = mle_to_pass, inc = inc[i],
                                    epsilon = epsilon[i],
                                    maxima_notNA = maxima_notNA,
                                    adjust = adjust, m = m[i],
                                    npy = npy)
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
      } else {
        ci_mat[i, ] <- matrix(NA, 1, 2)
        for_plot[[i]] <- NA
      }
      # Reset mult
      mult <- save_mult
    }
    for_plot <- for_plot[parm_numbers]
    names(for_plot) <- names(coef(object))[parm]
  }
  # Format the matrix to be returned
  low <- paste0(100 * (1 - level)/ 2, "%")
  up <- paste0(100 - 100 * (1 - level)/ 2, "%")
  colnames(ci_mat) <- c(low, up)
  ci_mat <- ci_mat[parm_numbers, , drop = FALSE]
  rownames(ci_mat) <- names(coef(object))[parm]
  # If profile = TRUE save the profile log-likelihood values for plotting
  if (profile) {
    attr(ci_mat, "for_plot") <- for_plot
    attr(ci_mat, "crit") <- conf_list$crit
    attr(ci_mat, "level") <- level
  }
  class(ci_mat) <- c("confint_return_level", "evmissing")
  return(ci_mat)
}
