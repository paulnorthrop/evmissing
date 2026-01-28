#' Internal evmissing functions
#'
#' Internal evmissing functions
#' @details
#' These functions are not intended to be called by the user.
#' @name evmissing-internal
#' @keywords internal
NULL

# ============================ GEV log-likelihood =========================== #

#' @keywords internal
#' @rdname evmissing-internal
negated_gev_loglik <- function(parameters, maxima_notNA, adjust,
                               big_val = Inf) {
  # Extract the block maxima, the numbers of non-missing values in the blocks
  # and the the largest possible number of non-missing values in each block
  maxima <- maxima_notNA$maxima
  n_i <- maxima_notNA$notNA
  n <- maxima_notNA$n
  # Extract the GEV parameter values for a complete block
  mu <- parameters[1]
  sigma <- parameters[2]
  xi <- parameters[3]
  # Infer the GEV parameter values for all blocks
  # If adjust = TRUE then we adjust the location and scale parameters for
  # block i based on the proportion p_i of non-missing values in the block.
  # Otherwise, leave the location and scale parameters as they are.
  if (adjust) {
    p_i <- n_i / n
    mu <- mu + sigma * box_cox_vec(x = p_i, lambda = xi)
    sigma <- sigma * p_i ^ xi
  }
  # Check that the parameters are not out-of-bounds
  if (any(sigma <= 0) || any(1 + xi * (maxima - mu) / sigma <= 0)) {
    return(big_val)
  }
  # Use the nieve package to calculate the contributions to the log-likelihood
  loglik <- nieve::dGEV(x = maxima, loc = mu, scale = sigma, shape = xi,
                        log = TRUE)
  # Return the sum of the contributions
  return(-sum(loglik))
}

# ======================== Weighted GEV log-likelihood ====================== #

#' @keywords internal
#' @rdname evmissing-internal
weighted_negated_gev_loglik <- function(parameters, maxima, weights,
                                        big_val = Inf) {
  # Extract the GEV parameter values for a complete block
  mu <- parameters[1]
  sigma <- parameters[2]
  xi <- parameters[3]
  # Check that the parameters are not out-of-bounds
  if (any(sigma <= 0) || any(1 + xi * (maxima - mu) / sigma <= 0)) {
    return(big_val)
  }
  # Use the nieve package to calculate the contributions to the log-likelihood
  loglik <- nieve::dGEV(x = maxima, loc = mu, scale = sigma, shape = xi,
                        log = TRUE)
  # Multiply by the weights
  loglik <- weights * loglik
  # Return the sum of the contributions
  return(-sum(loglik))
}

# ============== Discard block maxima based on many missings ================ #

#' @keywords internal
#' @rdname evmissing-internal
discard_maxima <- function(maxima_notNA, discard) {
  # maxima_notNA is a list with the following components
  #  * maxima:  the block maxima
  #  * notNA : the number of non-missing raw values in the block
  #  * n : the block size, the largest possible number of non-missing values
  # We discard maxima from all blocks that have more than discard% of the raw
  # values missing
  retain <- 100 * (1 - maxima_notNA$notNA / maxima_notNA$n) <= discard
  return(lapply(maxima_notNA, FUN = function(x) x[retain]))
}

# ==================== GEV log-likelihood for return levels ================= #

#' @keywords internal
#' @rdname evmissing-internal
negated_gev_loglik_ret_levs <- function(parameters, maxima_notNA, adjust, m,
                                        npy, big_val = Inf) {
  # Extract the block maxima and the numbers of non-missing values in the blocks
  maxima <- maxima_notNA$maxima
  n_i <- maxima_notNA$notNA
  n <- maxima_notNA$n
  # Extract the parameter values: (return level zp, sigma, xi)
  zp <- parameters[1]
  sigma <- parameters[2]
  xi <- parameters[3]
  # Infer the value of the GEV location parameter mu
  # Set the annual probability of exceedance based on m and npy
  p <- 1 - (1 - 1 / m) ^ (1 / npy)
  mu <- zp + sigma * box_cox_vec(x = -log(1 - p), lambda = -xi)
  # Infer the GEV parameter values for all blocks
  # If adjust = TRUE then we adjust the location and scale parameters for
  # block i based on the proportion p_i of non-missing values in the block.
  # Otherwise, leave the location and scale parameters as they are.
  if (adjust) {
    p_i <- n_i / n
    mu <- mu + sigma * box_cox_vec(x = p_i, lambda = xi)
    sigma <- sigma * p_i ^ xi
  }
  # Check that the parameters are not out-of-bounds
  if (any(sigma <= 0) || any(1 + xi * (maxima - mu) / sigma <= 0)) {
    return(big_val)
  }
  # Use the nieve package to calculate the contributions to the log-likelihood
  loglik <- nieve::dGEV(x = maxima, loc = mu, scale = sigma, shape = xi,
                        log = TRUE)
  # Return the sum of the contributions
  return(-sum(loglik))
}

# =============== Weighted GEV log-likelihood for return levels ============= #

#' @keywords internal
#' @rdname evmissing-internal
weighted_negated_gev_loglik_ret_levs <- function(parameters, maxima, weights,
                                                 m, npy, big_val = Inf) {
  # Extract the parameter values: (return level zp, sigma, xi)
  zp <- parameters[1]
  sigma <- parameters[2]
  xi <- parameters[3]
  # Infer the value of the GEV location parameter mu
  # Set the annual probability of exceedance based on m and npy
  p <- 1 - (1 - 1 / m) ^ (1 / npy)
  mu <- zp + sigma * box_cox_vec(x = -log(1 - p), lambda = -xi)
  # Check that the parameters are not out-of-bounds
  if (any(sigma <= 0) || any(1 + xi * (maxima - mu) / sigma <= 0)) {
    return(big_val)
  }
  # Use the nieve package to calculate the contributions to the log-likelihood
  loglik <- nieve::dGEV(x = maxima, loc = mu, scale = sigma, shape = xi,
                        log = TRUE)
  # Multiply by the weights
  loglik <- weights * loglik
  # Return the sum of the contributions
  return(-sum(loglik))
}

# ======================== Faster GEV profiling function ==================== #

#' @keywords internal
#' @rdname evmissing-internal
faster_profile_ci <- function(negated_loglik_fn, which = 1, level, mle,
                              ci_sym_mat, inc, epsilon, ci_init, ...) {

  # Determine whether we are profiling with respect to a return level
  return_level <- grepl("level", names(mle)[1])

  # The -log-likelihood to profile over parameters in par other than par[which]
  profiling_fn <- function(par, par_which, ...) {
    # Place par_which in position which and the other parameters around it
    parameters <- rep_len(NA, length(par) + 1)
    parameters[which] <- par_which
    parameters[-which] <- par
    # Return the negated log-likelihood
    return(negated_loglik_fn(parameters, ...))
  }

  # The maximised log-likelihood and the MLE for the parameter of interest
  max_loglik <- -negated_loglik_fn(mle, ...)
  mle_which <- mle[which]
  # The horizontal line that determines the confidence limits
  conf_line <- max_loglik - 0.5 * stats::qchisq(level, 1)
  # Vectors to store values of the parameters (x1 and x2) and the values of
  # the profile log-likelihood (v1 and v2)
  v1 <- v2 <- x1 <- x2 <- NULL
  x2[1] <- x1[1] <- mle[which]
  v2[1] <- v1[1] <- max_loglik

  # Extract the sample maxima
  # Need to determine which function was used for the original fit:
  # gev_mle(), for which the data are in maxima_notNA, or
  # gev_Weighted(), for which the data are in maxima
  dots <- list(...)
  if (!is.null(dots$maxima_notNA)) {
    data <- dots$maxima_notNA$maxima
  } else {
    data <- dots$maxima
  }

  ### Upper tail ...

  # We start from the upper limit of the symmetric confidence interval, using
  # gev_profile_init() to set initial estimates of the GEV parameters other
  # than the parameter which

  # Call gev_profile_init() or use stored initial estimates for return levels
  if (return_level) {
    upper_init <- ci_init$upper_init
    par_which <- upper_init[1]
    init <- upper_init[-1]
  } else {
    par_which <- ci_sym_mat[which, 2]
    init <- call_gev_profile_init(data = data, which = which,
                                  par_value = par_which)
  }
  # Calculate the profile log-likelihood at the initial values
  # If this is greater than conf_line then we search upwards
  # If this is smaller than conf_line then we search downwards

  ii <- 2
  opt <- try(stats::optim(init, profiling_fn, par_which = par_which, ...),
             silent = TRUE)
  if (inherits(opt, "try-error")) {
    return(list(optim_error = attr(opt, "condition")))
  }
  my_val <- -opt$value
  x2[ii] <- par_which
  v2[ii] <- my_val

  # If my_val < conf_line then the profile log-likelihood has dropped below
  # the required level in one (big) step. We use linear interpolation between
  # this point and the MLE to move back (hopefully, just) above the level.
  # This should work because the profile log-likelihood should be convex.
  if (my_val < conf_line) {
    searched_upwards <- FALSE
    while_condition <- function(my_val) {
      return(my_val < conf_line)
    }
    temp <- stats::splinefun(c(v2[1], v2[2]),
                             c(x2[1], x2[2]),
                             method = "hyman")
    delta <- -(x2[2] - temp(conf_line))
  } else {
    searched_upwards <- TRUE
    while_condition <- function(my_val) {
      return(my_val > conf_line)
    }
    delta <- inc
  }
  sol <- opt$par
  while (while_condition(my_val)){
    par_which <- par_which + delta
    opt <- try(stats::optim(sol, profiling_fn, par_which = par_which, ...),
               silent = TRUE)
    if (inherits(opt, "try-error")) {
      return(list(optim_error = attr(opt, "condition")))
    }
    sol <- opt$par
    ii <- ii + 1
    x2[ii] <- par_which
    v2[ii] <- -opt$value
    my_val <- v2[ii]
  }
  # If we searched downwards then reorder the results
  if (!searched_upwards) {
    x2[-1] <- rev(x2[-1])
    v2[-1] <- rev(v2[-1])
  }
  # Save sol for possible use later by itp()
  sol_upper <- sol

  ### Lower tail ...

  # We start from the lower limit of the symmetric confidence interval, using
  # gev_profile_init() to set initial estimates of the GEV parameters other
  # than the parameter which

  # Call gev_profile_init() or use stored initial estimates for return levels
  if (return_level) {
    lower_init <- ci_init$lower_init
    par_which <- lower_init[1]
    init <- lower_init[-1]
  } else {
    par_which <- ci_sym_mat[which, 1]
    init <- call_gev_profile_init(data = data, which = which,
                                  par_value = par_which)
  }
  # Call gev_profile_init()
  # Calculate the profile log-likelihood at the initial values
  # If this is greater than conf_line then we search downwards
  # If this is smaller than conf_line then we search upwards

  ii <- 2
  opt <- try(stats::optim(init, profiling_fn, par_which = par_which, ...),
             silent = TRUE)
  if (inherits(opt, "try-error")) {
    return(list(optim_error = attr(opt, "condition")))
  }
  my_val <- -opt$value
  x1[ii] <- par_which
  v1[ii] <- my_val

  # If my_val < conf_line then the profile log-likelihood has dropped below
  # the required level in one (big) step. We use linear interpolation between
  # this point and the MLE to move back (hopefully, just) above the level.
  # This should work because the profile log-likelihood should be convex.
  if (my_val < conf_line) {
    searched_downwards <- FALSE
    while_condition <- function(my_val) {
      return(my_val < conf_line)
    }
    temp <- stats::splinefun(c(v1[1], v1[2]),
                             c(x1[1], x1[2]),
                             method = "hyman")
    delta <- -(temp(conf_line) - x1[2])
  } else {
    searched_downwards <- TRUE
    while_condition <- function(my_val) {
      return(my_val > conf_line)
    }
    delta <- inc
  }
  sol <- opt$par
  while (while_condition(my_val)){
    par_which <- par_which - delta
    opt <- try(stats::optim(sol, profiling_fn, par_which = par_which, ...),
               silent = TRUE)
    if (inherits(opt, "try-error")) {
      return(list(optim_error = attr(opt, "condition")))
    }
    sol <- opt$par
    ii <- ii + 1
    x1[ii] <- par_which
    v1[ii] <- -opt$value
    my_val <- v1[ii]
  }
  # If we searched upwards then reorder the results
  if (!searched_downwards) {
    x1[-1] <- rev(x1[-1])
    v1[-1] <- rev(v1[-1])
  }
  # Save sol for possible use later by itp()
  sol_lower <- sol

  # Find the limits of the confidence interval
  prof_lik <- c(rev(v1), v2)
  par_values <- c(rev(x1), x2)

  # Find where the curve crosses conf_line
  temp <- diff(prof_lik - conf_line > 0)
  # Find the upper limit of the confidence interval
  loc_upper <- which(temp == -1)
  x1up <- par_values[loc_upper]
  x2up <- par_values[loc_upper + 1]
  y1up <- prof_lik[loc_upper]
  y2up <- prof_lik[loc_upper + 1]
  # Find the lower limit of the confidence interval
  loc_lower <- which(temp == 1)
  x1low <- par_values[loc_lower]
  x2low <- par_values[loc_lower + 1]
  y1low <- prof_lik[loc_lower]
  y2low <- prof_lik[loc_lower + 1]

  # If epsilon = 0 use linear interpolation
  # If epsilon < 0 use monotonic cubic spline interpolation
  # If epsilon > 0 use monotonic cubic spline interpolation then itp::itp()

  lower_pars <- NULL
  upper_pars <- NULL
  up_lim <- x1up + (conf_line - y1up) * (x2up - x1up) / (y2up - y1up)
  low_lim <- x1low + (conf_line - y1low) * (x2low - x1low) / (y2low - y1low)
  if (epsilon != 0) {
    # Calculate the values of the profile log-likelihood at these limits and
    # use the 3 points (the bracketing points and this new point) to estimate
    # the confidence limits by monotonic cubic spline interpolation

    # Upper
    opt <- try(stats::optim(sol_upper, profiling_fn,
                            par_which = up_lim, ...), silent = TRUE)
    # If optim errors then reset the initial values and try again
    if (inherits(opt, "try-error")) {
      init <- call_gev_profile_init(data = data, which = which,
                                    par_value = up_lim)
      opt <- try(stats::optim(init, profiling_fn, par_which = up_lim, ...),
                 silent = TRUE)
    }
    up_new <- -opt$value
    temp <- stats::splinefun(c(y1up, up_new, y2up),
                             c(x1up, up_lim, x2up),
                             method = "hyman")
    save_up_lim <- up_lim
    up_lim <- temp(conf_line)

    # Lower
    opt <- try(stats::optim(sol_lower, profiling_fn,
                            par_which = low_lim, ...), silent = TRUE)
    # If optim errors then reset the initial values and try again
    if (inherits(opt, "try-error")) {
      init <- call_gev_profile_init(data = data, which = which,
                                    par_value = low_lim)
      opt <- try(stats::optim(init, profiling_fn, par_which = low_lim, ...),
                 silent = TRUE)
    }
    low_new <- -opt$value
    temp <- stats::splinefun(c(y1low, low_new, y2low),
                             c(x1low, low_lim, x2low),
                             method = "hyman")
    save_low_lim <- low_lim
    low_lim <- temp(conf_line)

    # Add the approximate solutions to par_values and prof_lik
    # Note that only the extreme ends of these vectors have prof_lik values
    # below conf_line
    # Only do this if epsilon < 0 to avoid messing up the ordering of points
    # when we do something similar below for the epsilon > 0 case
    if (epsilon < 0) {
      n <- length(par_values)
      par_values <- c(par_values[1:loc_lower], save_low_lim,
                      par_values[(loc_lower + 1):loc_upper], save_up_lim,
                      par_values[(loc_upper + 1):n])
      prof_lik <- c(prof_lik[1:loc_lower], low_new,
                    prof_lik[(loc_lower + 1):loc_upper],
                    up_new, prof_lik[(loc_upper + 1):n])
    }

    # If epsilon > 0 then use itp::itp(), creating a new bracket from 2 of the
    # 3 points available and set an initial estimate
    if (epsilon > 0) {
      itp_function <- function(x, par, ...) {
        opt <- try(stats::optim(par = par, fn = profiling_fn, par_which = x,
                                ...), silent = TRUE)
        if (inherits(opt, "try-error")) {
          init <- call_gev_profile_init(data = data, which = which,
                                        par_value = x)
          opt <- try(stats::optim(par = init, fn = profiling_fn, par_which = x,
                                  ...), silent = TRUE)
        }
        val <- -opt$value - conf_line
        attr(val, "gev_pars") <- opt$par
        return(val)
      }
      # Find the upper limit of the confidence interval
      if (up_new > conf_line) {
        interval <- c(x1up, save_up_lim)
      } else {
        interval <- c(save_up_lim, x2up)
      }
      upper <- itp::itp(f = itp_function, interval = interval,
                        par = sol_upper,
                        f.a = y1up - conf_line, f.b = y2up - conf_line,
                        epsilon = epsilon, ...)
      up_lim <- upper$root
      # Find the lower limit of the confidence interval
      if (low_new > conf_line) {
        interval <- c(x1low, save_low_lim)
      } else {
        interval <- c(save_low_lim, x2low)
      }
      lower <- itp::itp(f = itp_function, interval = interval,
                        par = sol_lower,
                        f.a = y1low - conf_line, f.b = y2low - conf_line,
                        epsilon = epsilon, ...)
      low_lim <- lower$root
      # Add the approximate solutions to par_values and prof_lik
      # Note that only the extreme ends of these vectors have prof_lik values
      # below conf_line
      n <- length(par_values)
      par_values <- c(par_values[1:loc_lower], low_lim,
                      par_values[(loc_lower + 1):loc_upper], up_lim,
                      par_values[(loc_upper + 1):n])
      prof_lik <- c(prof_lik[1:loc_lower], lower$f.root + conf_line,
                    prof_lik[(loc_lower + 1):loc_upper],
                    upper$f.root + conf_line, prof_lik[(loc_upper + 1):n])
      # Save the parameter values that apply to the solutions from itp::itp()
      lower_pars <- numeric(3)
      lower_pars[which] <- lower$root
      lower_pars[-which] <- attr(lower$f.root, "gev_pars")
      names(lower_pars) <- c("mu", "sigma", "xi")
      upper_pars <- numeric(3)
      upper_pars[which] <- upper$root
      upper_pars[-which] <- attr(upper$f.root, "gev_pars")
      names(upper_pars) <- c("mu", "sigma", "xi")
    }
  }

  par_prof <- c(lower = low_lim, mle_which, upper = up_lim)
  return(list(par_prof = par_prof, crit = conf_line,
              for_plot = cbind(par_values = par_values,
                               prof_loglik = prof_lik),
              lower_pars = lower_pars, upper_pars = upper_pars))
}

# ======================== General profiling function ======================= #

#' @keywords internal
#' @rdname evmissing-internal
profile_ci <- function(negated_loglik_fn, which = 1, level, mle, inc, epsilon,
                       ...) {

  # The -log-likelihood to profile over parameters in par other than par[which]
  profiling_fn <- function(par, par_which, ...) {
    # Place par_which in position which and the other parameters around it
    parameters <- rep_len(NA, length(par) + 1)
    parameters[which] <- par_which
    parameters[-which] <- par
    # Return the negated log-likelihood
    return(negated_loglik_fn(parameters, ...))
  }

  # The maximised log-likelihood and the MLE for the parameter of interest
  max_loglik <- -negated_loglik_fn(mle, ...)
  mle_which <- mle[which]
  # The horizontal line that determines the confidence limits
  conf_line <- max_loglik - 0.5 * stats::qchisq(level, 1)
  # Vectors to store values of the parameters (x1 and x2) and the values of
  # the profile log-likelihood (v1 and v2)
  v1 <- v2 <- x1 <- x2 <- NULL
  x2[1] <- x1[1] <- mle[which]
  v2[1] <- v1[1] <- max_loglik

  # Extract the sample maxima
  # Need to determine which function was used for the original fit:
  # gev_mle(), for which the data are in maxima_notNA, or
  # gev_Weighted(), for which the data are in maxima
  dots <- list(...)
  if (!is.null(dots$maxima_notNA)) {
    data <- dots$maxima_notNA$maxima
  } else {
    data <- dots$maxima
  }

  # Starting from the MLE, we search upwards and downwards until we pass the
  # cutoff for the 100level% confidence interval

  ### Upper tail ...
  par_which <- mle_which
  my_val <- max_loglik
  ii <- 1
  sol <- mle[-which]
  while (my_val > conf_line){
    par_which <- par_which + inc
    opt <- try(stats::optim(sol, profiling_fn, par_which = par_which, ...),
               silent = TRUE)
    if (inherits(opt, "try-error")) {
      return(list(optim_error = attr(opt, "condition")))
    }
    sol <- opt$par
    ii <- ii + 1
    x2[ii] <- par_which
    v2[ii] <- -opt$value
    my_val <- v2[ii]
  }
  # Save sol for possible use later by itp()
  sol_upper <- sol
  # Also save the values of (mu, sigma, xi) after the profile log-likelihood
  # has dropped below conf_line. These may be useful in providing initial
  # values for profiling with respect to a return level
  upper_pars <- numeric(3)
  upper_pars[which] <- par_which
  upper_pars[-which] <- sol_upper
  names(upper_pars) <- c("mu", "sigma", "xi")

  ### Lower tail ...
  par_which <- mle_which
  my_val <- max_loglik
  ii <- 1
  sol <- mle[-which]
  while (my_val > conf_line){
    par_which <- par_which - inc
    opt <- try(stats::optim(sol, profiling_fn, par_which = par_which, ...),
               silent = TRUE)
    if (inherits(opt, "try-error")) {
      return(list(optim_error = attr(opt, "condition")))
    }
    sol <- opt$par
    ii <- ii + 1
    x1[ii] <- par_which
    v1[ii] <- -opt$value
    my_val <- v1[ii]
  }
  # Save sol for possible use later by itp()
  sol_lower <- sol
  # Also save the values of (mu, sigma, xi) after the profile log-likelihood
  # has dropped below conf_line. These may be useful in providing initial
  # values for profiling with respect to a return level
  lower_pars <- numeric(3)
  lower_pars[which] <- par_which
  lower_pars[-which] <- sol_lower
  names(lower_pars) <- c("mu", "sigma", "xi")

  # Find the limits of the confidence interval
  prof_lik <- c(rev(v1), v2)
  par_values <- c(rev(x1), x2)

  # Find where the curve crosses conf_line
  temp <- diff(prof_lik - conf_line > 0)
  # Find the upper limit of the confidence interval
  loc_upper <- which(temp == -1)
  x1up <- par_values[loc_upper]
  x2up <- par_values[loc_upper + 1]
  y1up <- prof_lik[loc_upper]
  y2up <- prof_lik[loc_upper + 1]
  # Find the lower limit of the confidence interval
  loc_lower <- which(temp == 1)
  x1low <- par_values[loc_lower]
  x2low <- par_values[loc_lower + 1]
  y1low <- prof_lik[loc_lower]
  y2low <- prof_lik[loc_lower + 1]

  # If epsilon = 0 use linear interpolation
  # If epsilon < 0 use monotonic cubic spline interpolation
  # If epsilon > 0 use monotonic cubic spline interpolation then itp::itp()

  up_lim <- x1up + (conf_line - y1up) * (x2up - x1up) / (y2up - y1up)
  low_lim <- x1low + (conf_line - y1low) * (x2low - x1low) / (y2low - y1low)
  if (epsilon != 0) {
    # Calculate the values of the profile log-likelihood at these limits and
    # use the 3 points (the bracketing points and this new point) to estimate
    # the confidence limits by monotonic cubic spline interpolation

    # Upper
    opt <- try(stats::optim(sol_upper, profiling_fn,
                            par_which = up_lim, ...), silent = TRUE)
    # If optim errors then reset the initial values and try again
    if (inherits(opt, "try-error")) {
      init <- call_gev_profile_init(data = data, which = which,
                                    par_value = up_lim)
      opt <- try(stats::optim(init, profiling_fn, par_which = up_lim, ...),
                 silent = TRUE)
    }
    up_new <- -opt$value
    temp <- stats::splinefun(c(y1up, up_new, y2up),
                             c(x1up, up_lim, x2up),
                             method = "hyman")
    save_up_lim <- up_lim
    up_lim <- temp(conf_line)

    # Lower
    opt <- try(stats::optim(sol_lower, profiling_fn,
                            par_which = low_lim, ...), silent = TRUE)
    # If optim errors then reset the initial values and try again
    if (inherits(opt, "try-error")) {
      init <- call_gev_profile_init(data = data, which = which,
                                    par_value = low_lim)
      opt <- try(stats::optim(init, profiling_fn, par_which = low_lim, ...),
                 silent = TRUE)
    }
    low_new <- -opt$value
    temp <- stats::splinefun(c(y1low, low_new, y2low),
                             c(x1low, low_lim, x2low),
                             method = "hyman")
    save_low_lim <- low_lim
    low_lim <- temp(conf_line)

    # If epsilon > 0 then use itp::itp(), creating a new bracket from 2 of the
    # 3 points available and set an initial estimate
    if (epsilon > 0) {
      itp_function <- function(x, par, ...) {
        opt <- try(stats::optim(par = par, fn = profiling_fn, par_which = x,
                                ...), silent = TRUE)
        if (inherits(opt, "try-error")) {
          init <- call_gev_profile_init(data = data, which = which,
                                        par_value = x)
          opt <- try(stats::optim(par = init, fn = profiling_fn, par_which = x,
                                  ...), silent = TRUE)
        }
        val <- -opt$value - conf_line
        attr(val, "gev_pars") <- opt$par
        return(val)
      }
      # Find the upper limit of the confidence interval
      if (up_new > conf_line) {
        interval <- c(x1up, save_up_lim)
      } else {
        interval <- c(save_up_lim, x2up)
      }
      upper <- itp::itp(f = itp_function, interval = interval,
                        par = sol_upper,
                        f.a = y1up - conf_line, f.b = y2up - conf_line,
                        epsilon = epsilon, ...)
      up_lim <- upper$root
      # Find the lower limit of the confidence interval
      if (low_new > conf_line) {
        interval <- c(x1low, save_low_lim)
      } else {
        interval <- c(save_low_lim, x2low)
      }
      lower <- itp::itp(f = itp_function, interval = interval,
                        par = sol_lower,
                        f.a = y1low - conf_line, f.b = y2low - conf_line,
                        epsilon = epsilon, ...)
      low_lim <- lower$root
      # Add the approximate solutions to par_values and prof_lik
      # Note that only the extreme ends of these vectors have prof_lik values
      # below conf_line
      n <- length(par_values)
      par_values <- c(par_values[1:loc_lower], low_lim,
                      par_values[(loc_lower + 1):loc_upper], up_lim,
                      par_values[(loc_upper + 1):n])
      prof_lik <- c(prof_lik[1:loc_lower], lower$f.root + conf_line,
                    prof_lik[(loc_lower + 1):loc_upper],
                    upper$f.root + conf_line, prof_lik[(loc_upper + 1):n])
      # Save the parameter values that apply to the solutions from itp::itp()
      lower_pars <- numeric(3)
      lower_pars[which] <- lower$root
      lower_pars[-which] <- attr(lower$f.root, "gev_pars")
      names(lower_pars) <- c("mu", "sigma", "xi")
      upper_pars <- numeric(3)
      upper_pars[which] <- upper$root
      upper_pars[-which] <- attr(upper$f.root, "gev_pars")
      names(upper_pars) <- c("mu", "sigma", "xi")
    }
  }

  par_prof <- c(lower = low_lim, mle_which, upper = up_lim)
  return(list(par_prof = par_prof, crit = conf_line,
              for_plot = cbind(par_values = par_values,
                               prof_loglik = prof_lik),
              lower_pars = lower_pars, upper_pars = upper_pars))
}

# ========================== Box-Cox transformation ========================= #

#' @keywords internal
#' @rdname evmissing-internal
box_cox <- function (x, lambda = 1, gm = 1, lambda_tol = 1e-6) {
  #
  # Computes the Box-Cox transformation of a vector.
  #
  # Args:
  #   x          : A numeric vector. (Positive) values to be Box-Cox
  #                transformed.
  #   lambda     : A numeric scalar.  Transformation parameter.
  #   gm         : A numeric scalar.  Optional scaling parameter.
  #   lambda_tol : A numeric scalar.  For abs(lambda) < lambda_tol use
  #                a Taylor series expansion.
  #
  # Returns:
  #   A numeric vector.  The transformed value
  #     (x^lambda - 1) / (lambda * gm ^ (lambda - 1))
  #
  if (abs(lambda) > lambda_tol) {
    retval <- (x ^ lambda - 1) / lambda / gm ^ (lambda - 1)
  } else {
    i <- 0:3
    retval <- sum(log(x) ^ (i + 1) * lambda ^ i / factorial(i + 1))
    retval <- retval / gm ^ (lambda - 1)
  }
  retval
}

#' @keywords internal
#' @rdname evmissing-internal
box_cox_vec <- function(x, lambda = 1, lambda_tol = 1e-6) {
  #
  # Computes the Box-Cox transformation of a vector.  If lambda is very close
  # to zero then a first order Taylor series approximation is used.
  #
  # Args:
  #   x          : A numeric vector. (Non-negative) values to be Box-Cox
  #                transformed.
  #   lambda     : A numeric scalar.  Transformation parameter.
  #   lambda_tol : A numeric scalar.  For abs(lambda) < lambda_tol use
  #                a Taylor series expansion.
  # Returns:
  #   A numeric vector.  The transformed value
  #     (x^lambda - 1) / lambda
  #
  if (any(x < 0, na.rm = TRUE)) {
    stop("Invalid x: x must be non-negative")
  }
  max_len <- max(length(x), length(lambda))
  x <- rep_len(x, max_len)
  lambda <- rep_len(lambda, max_len)
  retval <- ifelse(abs(lambda) > lambda_tol, (x ^ lambda - 1) / lambda,
              ifelse(lambda == 0, log(x),
                ifelse(is.infinite(x),
                  ifelse(lambda < 0, -1 / lambda, Inf),
                    ifelse(x == 0, ifelse(lambda > 0, -1 / lambda, -Inf),
                      log(x) * (1 + log(x) * lambda / 2)))))
  return(retval)
}

# ============================== box_cox_deriv ============================== #

#' @keywords internal
#' @rdname evmissing-internal
box_cox_deriv <- function(x, lambda = 1, lambda_tol = 1 / 50,
                          poly_order = 3) {
  #
  # Computes the derivative with respect to lambda the Box-Cox
  # transformation.
  #
  # Args:
  #   x          : A numeric vector. (Positive) values to be Box-Cox
  #                transformed.
  #   lambda     : A numeric scalar.  Transformation parameter.
  #   lambda_tol : A numeric scalar.  For abs(lambda) < lambda_tol use
  #                a Taylor series expansion.
  #   poly_order : order of Taylor series polynomial in lambda used as
  #                an approximation if abs(lambda) < lambda_tol
  #
  # Returns:
  #   A numeric vector.  The derivative with respect to lambda of
  #     (x^lambda - 1) / lambda
  #
  lnx <- log(x)
  if (abs(lambda) > lambda_tol) {
    retval <- (lambda * x ^ lambda * lnx - x ^ lambda + 1) / lambda ^ 2
  } else {
    i <- 0:poly_order
    retval <- sum(lnx ^ (i + 2) * lambda ^ i / ((i + 2) * factorial(i)))
  }
  return(retval)
}

# =========================== GEV initial estimates ========================= #

#' @keywords internal
#' @rdname evmissing-internal
gev_init <- function(maxima_notNA, init_method = "quartiles") {
  # Check that init_method is appropriate
  init_method <- match.arg(init_method, c("quartiles", "moments"))
  # Extract the block maxima
  maxima <- maxima_notNA$maxima
  # If init_method = "quartiles" then base the initial estimates of location
  # and scale on the sample quantiles of block maxima, ignoring the underlying
  # numbers of non-missing raw data and a value of zero for the shape parameter.
  # If init_method = "moments" do similarly based on the sample mean and variance
  # of these maxima and an initial value of 0.1 for the shape parameter.
  if (init_method == "quartiles") {
    sample_quantiles <- stats::quantile(maxima, probs = 1:3/4)
    iqr <- sample_quantiles[3] - sample_quantiles[1]
    med <- sample_quantiles[2]
    sigma_init <- iqr / log(log(4) / log(4/3))
    mu_init <- med - sigma_init * log(1 / log(2))
    init <- c(mu_init, sigma_init, 0)
  } else {
    sigma_init <- sqrt(6 * stats::var(maxima)) / pi
    mu_init <- mean(maxima) - 0.57722 * sigma_init
    init <- c(mu_init, sigma_init, 0.1)
  }
  names(init) <- c("mu", "sigma", "xi")
  return(init)
}

# ================== Initial estimates for GEV profiling ==================== #

#' @keywords internal
#' @rdname evmissing-internal
gev_profile_init <- function(data, mu, sigma, xi) {
  # This function provides initial estimates of two of the GEV parameters
  # (mu, sigma, xi) given a user-supplied value for the other one.
  # The supplied value must be a numeric vector of length one.
  # This is NOT checked.

  # This may be helpful when profiling a GEV log-likelihood to find the
  # limits of a confidence interval for one of the parameters.
  # Suppose that a Wald-type symmetric confidence interval has already been
  # calculated for the parameter of interest. The limits of this interval could
  # be used to help start the search for the limits of the interval based on
  # the corresponding profile log-likelihood.

  # The estimates calculated below are based on the LRSE(EV) approach described
  # on page 99 of the first edition of the book Reiss and Thomas (1997)
  # Statistical Analysis of Extreme Values, Birkhauser, Basel, Switzerland.
  # The following uses the notation in this book.

  # We set q0 and q2 to correspond to the minimum and maximum of the observed
  # values in the input vector data. This enables us to ensure that the
  # resulting estimates respect the constraints on the GEV parameter space
  # given the observed values in data.

  # Which of the GEV parameters has been given? There should be exactly one
  mu_given <- !missing(mu)
  sigma_given <- !missing(sigma)
  xi_given <- !missing(xi)
  number_given <- sum(mu_given, sigma_given, xi_given)
  if (number_given != 1) {
    stop("Exactly one of mu, sigma, xi must be provided.")
  }

  # Calculate quantities that will be used throughout
  n <- length(data)
  # Values of q0, q1 and q2
  q0 <- 1 / (n + 1)
  q2 <- n / (n + 1)
  a <- sqrt(log(q2) / log(q0))
  q1 <- q0 ^ a
  # The corresponding order statistics
  n0 <- 1
  n1 <- round((n + 1) * q1)
  n2 <- n
  ns <- c(n0, n1, n2)
  qs <- c(q0, q1, q2)
  #  Extract these order statistics from the data
  x <- sort(data)[ns]

  # If xi is given then we set mu and sigma so that the fitted quantiles
  # for q0 and q2 equal the corresponding order statistics in x

  # If either mu or sigma is given then first we estimate xi. Then we use
  # either the sample minimum or maximum, depending on the sign of the estimate
  # of xi, to estimate the remaining unknown parameter.
  # If this sign is negative then we use the sample maximum.
  # If this sign is positive then we use the sample minimum.
  # This ensures that none of the data are outside the GEV support.

  if (xi_given) {
    fq0 <- box_cox(-log(q0), lambda = -xi)
    fq2 <- box_cox(-log(q2), lambda = -xi)
    sigma <- (x[3] - x[1]) / (fq0 - fq2)
    mu <- x[1] + sigma * fq0
  } else {
    r_hat <- (x[3] - x[2]) / (x[2] - x[1])
    xi <- -log(r_hat) / log(a)
    fq0 <- box_cox(-log(q0), lambda = -xi)
    fq2 <- box_cox(-log(q2), lambda = -xi)
    if (sigma_given) {
      if (xi < 0) {
        mu <- x[3] + sigma * fq2
      } else {
        mu <- x[1] + sigma * fq0
      }
    } else if (mu_given) {
      if (xi < 0) {
        sigma <- (mu - x[3]) / fq2
      } else {
        sigma <- (mu - x[1]) / fq0
      }
    }
  }
  val <- c(mu, sigma, xi)
  names(val) <- c("mu", "sigma", "xi")
  return(val)
}


# ======================== Calling gev_profile_init() ======================= #

#' @keywords internal
#' @rdname evmissing-internal
call_gev_profile_init <- function(data, which, par_value) {
  # which indicates the GEV parameter: mu, sigma or xi
  if (which == 1) {
    init <- gev_profile_init(data = data, mu = par_value)[-1]
  } else if (which == 2) {
    init <- gev_profile_init(data = data, sigma = par_value)[-2]
  } else {
    init <- gev_profile_init(data = data, xi = par_value)[-3]
  }
  return(init)
}


# =============== Initial estimates for return level profiling ============== #

#' @keywords internal
#' @rdname evmissing-internal
return_level_profile_init <- function(gev_object, level, mult, epsilon,
                                      m, npy) {
  # This function provides initial estimates of the GEV parameters sigma and xi
  # given a user-supplied value for the m-year return level.
  # The supplied value must be a numeric vector of length one.
  # This is NOT checked.

  # Call confint.evmissing() to obtain GEV parameter values for the
  # profile log-likelihood-based confidence limits.
  # These can help to produce good initial estimates for sigma and xi
  # Perhaps increase mult?
  # If faster = TRUE then we need epsilon > 0. Otherwise, estimates of
  # parameters from stats::optim() are not stored for use as initial estimates
  gev_prof <- confint(gev_object, level = level, profile = TRUE,
                      mult = mult, faster = TRUE, epsilon = epsilon)
  # Extract the GEV parameter values at the lower and upper limits
  lower_pars <- attr(gev_prof, "lower_pars")
  upper_pars <- attr(gev_prof, "upper_pars")
  # Calculate the m-year return levels implied by these GEV parameter values
  return_level_fn <- function(pars) {
    mu <- pars["mu"]
    sigma <- pars["sigma"]
    xi <- pars["xi"]
    p <- 1 - (1 - 1 / m) ^ (1 / npy)
    val <- nieve::qGEV(p, loc = mu, scale = sigma, shape = xi,
                       lower.tail = FALSE)
    return(val)
  }
  lower_return_levels <- lapply(lower_pars, FUN = return_level_fn)
  upper_return_levels <- lapply(upper_pars, FUN = return_level_fn)
  # Pick the most extreme (furthest from the MLE) cases for lower and upper
  where_min <- which.min(lower_return_levels)
  where_max <- which.max(upper_return_levels)
  lower_init <- unlist(lower_pars[where_min])
  upper_init <- unlist(upper_pars[where_max])
  # Replace the values of mu with the corresponding values of return level
  lower_init[1] <- lower_return_levels[where_min]
  upper_init[1] <- upper_return_levels[where_max]
  names(lower_init) <- c(paste0(m, "-year ", "return level"), "sigma", "xi")
  names(upper_init) <- c(paste0(m, "-year ", "return level"), "sigma", "xi")
  return(list(lower_init = unlist(lower_init),
              upper_init = unlist(upper_init)))
}

#' @keywords internal
#' @rdname evmissing-internal
gev_pp <- function (x, adjust, level, ...){
  # Extract the MLEs, block maxima and sample size
  mle <- coef(x)
  maxima <- x$maxima
  b <- nobs(x)
  # Standard uniform quantiles
  model <- (1:b) / (b + 1)

  # Adjust the GEV parameters for each block maximum based on the number of
  # non-missing raw observations
  n_i <- x$notNA
  n <- x$n
  # Extract the GEV parameter values for a complete block
  mu <- mle["mu"]
  sigma <- mle["sigma"]
  xi <- mle["xi"]
  # Infer the GEV parameter values for all blocks
  if (adjust) {
    p_i <- n_i / n
    mu <- mu + sigma * box_cox_vec(x = p_i, lambda = xi)
    sigma <- sigma * p_i ^ xi
  }

  # Calculate the GEV cdf for each block maximum
  # We need to retain the order of the block maxima for this calculation
  empirical <- nieve::pGEV(q = maxima,
                           loc = mu, scale = sigma, shape = xi,
                           lower.tail = TRUE)
  # Calculate uniform 100*level% pointwise confidence bands
  bands <- uniform_confidence_bands(b = b, k = 1:b, level = level)
  # Create the plot
  pp_and_qq_plot_fn(x = model, y = sort(empirical), bands = bands,
                    which = "pp", ..., xlab = "uniform quantiles",
                    main = "probability plot")
  return(invisible())
}

#' @keywords internal
#' @rdname evmissing-internal
gev_qq <- function (x, adjust, level, ...){
  # Extract the MLEs, block maxima and sample size
  mle <- coef(x)
  maxima <- x$maxima
  b <- nobs(x)

  # Extract the GEV parameter values for a complete block
  mu_full <- mle["mu"]
  sigma_full <- mle["sigma"]
  xi <- mle["xi"]

  # GEV quantiles for a 'full' block
  model <- nieve::qGEV(p = (1:b) / (b + 1),
                       loc = mu_full, scale = sigma_full, shape = xi,
                       lower.tail = TRUE)

  # If adjust = TRUE transform maxima to have a 'full-block' GEV distribution.
  # Adjust the GEV parameters for each block maximum based on the number of
  # non-missing raw observations.
  # Otherwise, use the maxima untransformed.
  n_i <- x$notNA
  n <- x$n
  # Infer the GEV parameter values for all blocks
  if (adjust) {
    p_i <- n_i / n
    mu <- mu_full + sigma_full * box_cox_vec(x = p_i, lambda = xi)
    sigma <- sigma_full * p_i ^ xi
    # Calculate the GEV cdf for each block maximum
    # We need to retain the order of the block maxima for this calculation
    empirical <- nieve::pGEV(q = maxima,
                             loc = mu, scale = sigma, shape = xi,
                             lower.tail = TRUE)
    empirical <- nieve::qGEV(p = empirical,
                             loc = mu_full, scale = sigma_full, shape = xi,
                             lower.tail = TRUE)
  } else {
    empirical <- maxima
  }

  # Calculate fitted GEV 100*level% pointwise confidence bands
  bands <- gev_confidence_bands(b = b, k = 1:b, level = level, mle = mle)
  # Set ylim to include the bands
  ylim <- c(min(bands), max(bands))
  # Create the plot
  pp_and_qq_plot_fn(x = model, y = sort(empirical), bands = bands,
                    which = "qq", ylim = ylim, ..., xlab = "GEV quantiles",
                    main = "quantile plot")
  return(invisible())
}

#' @keywords internal
#' @rdname evmissing-internal
gev_rl <- function (x, adjust, m, level, profile, num, npy, ...){
  # Extract the MLEs, block maxima and sample size and covariance matrix
  mle <- coef(x)
  maxima <- x$maxima
  b <- nobs(x)
  covariance <- vcov(x)

  # Extract the GEV parameter values for a complete block
  mu_full <- mle["mu"]
  sigma_full <- mle["sigma"]
  xi <- mle["xi"]

  # If adjust = TRUE transform maxima to have a 'full-block' GEV distribution.
  # Adjust the GEV parameters for each block maximum based on the number of
  # non-missing raw observations.
  # Otherwise, use the maxima untransformed.
  n_i <- x$notNA
  n <- x$n
  # Infer the GEV parameter values for all blocks
  if (adjust) {
    p_i <- n_i / n
    mu <- mu_full + sigma_full * box_cox_vec(x = p_i, lambda = xi)
    sigma <- sigma_full * p_i ^ xi
    # Calculate the GEV cdf for each block maximum
    # We need to retain the order of the block maxima for this calculation
    empirical <- nieve::pGEV(q = maxima,
                             loc = mu, scale = sigma, shape = xi,
                             lower.tail = TRUE)
    empirical <- nieve::qGEV(p = empirical,
                             loc = mu_full, scale = sigma_full, shape = xi,
                             lower.tail = TRUE)
  } else {
    empirical <- maxima
  }

  # Base the return periods to plot and label on m
  min_r <- min(min(m), (b + 1) / b)
  max_r <- max(max(m), b + 1)

  # Space the values on the horizontal axis regularly on the
  # log_10(-1 / log(f)) scale, that is, the
  # log_10(-1 / log(1 - 1 / return_period)) scale
  min_f_scale <- log(-1 / log(1 - 1 / min_r), base = 10)
  max_f_scale <- log(-1 / log(1 - 1 / max_r), base = 10)
  f_values <- 10 ^ seq(from = min_f_scale, to = max_f_scale, length.out = num)
  # Infer the values of f required below
  f <- exp(-1 / f_values)

  # Return level estimates
  ret_levs <- gev_return(x, m = 1 / (1 - f), npy = npy)
  # Return level confidence intervals
  # If profile = TRUE then the arguments mult and faster are set to improve the
  # speed of calculation
  cis <- confint(ret_levs, level = level, profile = profile,
                 mult = 32, faster = TRUE)
  lower_limits <- cis[, 1]
  upper_limits <- cis[, 2]
  # Create the plot
  xlim <- range(-1 / log(f), m)
  ylim <- c(min(empirical, ret_levs, lower_limits),
            max(empirical, ret_levs, upper_limits))
  horizontal_axis_values <- - 1 / log(1 - 1 / m)
  y_mat <- cbind(lower_limits, ret_levs, upper_limits)
  return_level_plot_fn(x = -1 / log(f), y = y_mat, empirical = empirical,
                       at = horizontal_axis_values, labels = m,
                       xlim = xlim, ylim = ylim, ...)
  return(invisible(y_mat))
}

#' @keywords internal
#' @rdname evmissing-internal
gev_his <- function (x, adjust, ...){
  # Extract the MLEs, block maxima and sample size
  mle <- coef(x)
  maxima <- x$maxima
  b <- nobs(x)

  # Extract the GEV parameter values for a complete block
  mu_full <- mle["mu"]
  sigma_full <- mle["sigma"]
  xi <- mle["xi"]

  # GEV quantiles for a 'full' block
  model <- nieve::qGEV(p = (1:b) / (b + 1),
                       loc = mu_full, scale = sigma_full, shape = xi,
                       lower.tail = TRUE)

  # If adjust = TRUE transform maxima to have a 'full-block' GEV distribution.
  # Adjust the GEV parameters for each block maximum based on the number of
  # non-missing raw observations.
  # Otherwise, use the maxima untransformed.
  n_i <- x$notNA
  n <- x$n
  # Infer the GEV parameter values for all blocks
  if (adjust) {
    p_i <- n_i / n
    mu <- mu_full + sigma_full * box_cox_vec(x = p_i, lambda = xi)
    sigma <- sigma_full * p_i ^ xi
    # Calculate the GEV cdf for each block maximum
    # We need to retain the order of the block maxima for this calculation
    empirical <- nieve::pGEV(q = maxima,
                             loc = mu, scale = sigma, shape = xi,
                             lower.tail = TRUE)
    empirical <- nieve::qGEV(p = empirical,
                             loc = mu_full, scale = sigma_full, shape = xi,
                             lower.tail = TRUE)
  } else {
    empirical <- maxima
  }

  # The relevant GEV parameters are those for the full block: mle above and
  # stored in (mu_full, sigma_full, xi)
  density_plot_fn(x = empirical, mle = mle, mu_full, sigma_full,
                  ...)
  return(invisible())
}

#' @keywords internal
#' @rdname evmissing-internal
pp_and_qq_plot_fn <- function(x, y, bands, which, ...,
                              xlab, ylab = "empirical", main) {
  # Plot the points
  dots <- list(...)
  dots <- c(dots, list(x = x, y = y, xlab = xlab, ylab = ylab, main = main))
  # Force all the points to be black and unaffected by lwd
  dots$col <- NULL
  dots$lwd <- NULL
  # Add the lines of equality and pointwise confidence bands
  line_dots <- list(...)
  y_values <- cbind(bands[, 1], x, bands[,2])
  line_dots <- c(line_dots, list(x = x, y = y_values))
  # Set ylim
  if (which == "pp") {
    dots$ylim = c(0, 1)
  } else {
    dots$ylim <- range(y_values)
  }
  do.call(graphics::plot, dots)
  # Set default line colours of blue for CI limits and black for return levels
  if (is.null(line_dots$col)) {
    line_dots$col <- c(4, 1, 4)
  }
  if (is.null(line_dots$lwd)) {
    line_dots$lwd <- 2
  }
  if (is.null(line_dots$lty)) {
    line_dots$lty <- 1
  }
  do.call(graphics::matlines, line_dots)
  # Replot the points to put them in the foreground
  do.call(graphics::points, dots)
  return(invisible())
}

#' @keywords internal
#' @rdname evmissing-internal
return_level_plot_fn <- function(x, y, empirical, at, labels, ..., type = "l",
                                 xlab = "return period", ylab = "return level",
                                 main = "return level plot", xaxt = "n") {
  # Create the plot area
  dots <- list(...)
  dots <- c(dots, list(x = x, y = y, type = "n", log = "x", xlab = xlab,
                       ylab = ylab, main = main, xaxt = xaxt))
  # Force all the points to be black and unaffected by lwd
  dots$col <- NULL
  dots$lwd <- NULL
  do.call(graphics::matplot, dots)
  graphics::axis(1, at = at, labels = labels)
  # Add the return level estimates and CI limits
  line_dots <- list(...)
  line_dots <- c(line_dots, list(x = x, y = y))
  # Set default line colours of blue for CI limits and black for return levels
  if (is.null(line_dots$col)) {
    line_dots$col <- c(4, 1, 4)
  }
  if (is.null(line_dots$lwd)) {
    line_dots$lwd <- 2
  }
  if (is.null(line_dots$lty)) {
    line_dots$lty <- 1
  }
  # Plot the return level estimates and CIs
  do.call(graphics::matlines, line_dots)
  # Add the observations
  n <- length(empirical)
  dots <- list(...)
  dots$x <- -1/log((1:n)/(n + 1))
  dots$y <- sort(empirical)
  dots$col <- "black"
  do.call(graphics::points, dots)
  return(invisible())
}

#' @keywords internal
#' @rdname evmissing-internal
density_plot_fn <- function(x, mle, mu_full, sigma_full, ...,
                            xlab = "block maxima", ylab = "density",
                            main = "density plot") {
  # Call hist() but don't produce the plot yet
  h <- graphics::hist(x, plot = FALSE)
  xi <- mle[3]
  # Set xlim to cover at least the middle 98% of the fitted GEV density
  gev_qs <- nieve::qGEV(p = c(0.01, 0.99), loc = mu_full, scale = sigma_full,
                        shape = xi)
  if (xi < 0) {
    x_values <- seq(min(h$breaks, gev_qs[1]),
                    min(max(h$breaks, gev_qs[2]),
                        (mle[1] - mle[2] / mle[3] - 0.001)),
                    length = 200)
  } else {
    x_values <- seq(max(min(h$breaks, gev_qs[1]),
                        (mle[1] - mle[2] / mle[3] + 0.001)),
                    max(h$breaks, gev_qs[2]),
                    length = 200)
  }
  gev_dens <- nieve::dGEV(x = x_values, loc = mu_full, scale = sigma_full,
                          shape = xi)
  ylim <- c(0, max(gev_dens, h$density))
  dots <- list(...)
  dots <- c(dots, list(x = x, xlab = xlab, ylab = ylab, main = main,
                       ylim = ylim, xlim = range(x_values), prob = TRUE))
  dots$lwd <- 1
  dots$col <- "lightgrey"
  dots$lty <- 1
  # Plot the histogram
  do.call(graphics::hist, dots)
  # Add the fitted GEV density
  line_dots <- list(...)
  line_dots <- c(line_dots, list(x = x_values, y = gev_dens))
  # Set the default to be a black solid line of width 2
  line_dots$col <- "black"
  line_dots$lty <- 1
  if (is.null(line_dots$lwd)) {
    line_dots$lwd <- 2
  }
  do.call(graphics::lines, line_dots)
  # Add the observations
  dots <- list(...)
  dots$x <- x
  dots$y  <- rep(0, length(x))
  dots$col <- "black"
  do.call(graphics::points, dots)
  return(invisible())
}

#' @keywords internal
#' @rdname evmissing-internal
uniform_confidence_bands <- function(b, k, level) {
  p <- (1 - level) / 2
  bands <- matrix(NA, ncol = 2, nrow = length(k))
  bands[, 1] <- stats::qbeta(p = p, shape1 = k, shape2 = b - k + 1)
  bands[, 2] <- stats::qbeta(p = 1 - p, shape1 = k, shape2 = b - k + 1)
  return(bands)
}

#' @keywords internal
#' @rdname evmissing-internal
gev_confidence_bands <- function(b, k, level, mle) {
  p <- (1 - level) / 2
  gev_bands <- matrix(NA, ncol = 2, nrow = length(k))
  gev_bands[, 1] <- stats::qbeta(p = p, shape1 = k, shape2 = b - k + 1)
  gev_bands[, 2] <- stats::qbeta(p = 1 - p, shape1 = k, shape2 = b - k + 1)
  gev_bands[, 1] <- nieve::qGEV(p = gev_bands[, 1], loc = mle["mu"],
                                scale = mle["sigma"], shape = mle["xi"],
                                lower.tail = TRUE)
  gev_bands[, 2] <- nieve::qGEV(p = gev_bands[, 2], loc = mle["mu"],
                                scale = mle["sigma"], shape = mle["xi"],
                                lower.tail = TRUE)
  return(gev_bands)
}

#' @keywords internal
#' @rdname evmissing-internal
merge_two_lists <- function(list1, list2) {
  merged_list <- list1
  merged_list[names(list2)[!names(list2) %in% names(list1)]] <-
    list2[!names(list2) %in% names(list1)]
  return(merged_list)
}
