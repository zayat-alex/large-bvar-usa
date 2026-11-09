# tilting.R: soft conditions on forecasts by entropic tilting.
#
# Why it exists: a hard condition fixes a cell of the forecast exactly
# (R/forecast_conditional.R). Often a user's belief is softer:
# "4-quarter PCE inflation three years ahead is around 2%, with the usual
# uncertainty". Crump et al. (2021), Section 4.2, impose such beliefs by
# entropic tilting (Robertson, Tallman and Whiteman 2005; Kruger, Clark and
# Ravazzolo 2017): keep the forecast paths as they are, and change only how
# much each path counts. Among all reweightings that meet the targets, take
# the one closest to equal weights in Kullback-Leibler divergence. The
# model's own spread around the target survives, reweighted, instead of
# being collapsed to a point.
#
# Notation (paper eqs. 14-16): J forecast paths y^(j), each with weight 1/J;
# m targets written as moment conditions E[g(y)] = gbar, with g a function
# of a path (g_j = g(y^(j)), an m-vector). The tilted weights solve
#   min_w sum_j w_j log(J w_j)  s.t.  sum_j w_j g_j = gbar,  sum_j w_j = 1   (14)
# and are
#   w_j = exp(gamma' g_j) / sum_k exp(gamma' g_k)                             (15)
#   gamma = argmin_gamma sum_j exp(gamma' (g_j - gbar))                        (16)
# The objective of (16) is convex in gamma, so Newton's method finds its
# minimum when the targets can be met (gbar strictly inside the convex hull
# of the g_j); otherwise there is no minimum and the tilting stops.
#
# Targets: a target is a statistic of a quantity s(y) computed on each path:
#   median of s at x:  g = 1[s <= x], gbar = 0.5 (the paper's anchor);
#   mean of s at x:    g = s,          gbar = x.
# The quantity is a variable in one forecast quarter, as a level (model
# units) or as an exact percent change (as path_from_growth() in
# R/scenarios.R): yoy = 100 (x_t / x_{t-4} - 1), qoq_ann =
# 100 [(x_t / x_{t-1})^4 - 1], or the growth of an annual average,
# annual_avg = 100 (sum_{i=0..3} x_{t-i} / sum_{i=4..7} x_{t-i} - 1), with
# x = exp(y / 100). In a Q4 quarter annual_avg is the growth of the
# calendar-year average over the previous year's, the form in which
# outlooks such as the Survey of Professional Forecasters quote annual
# real GDP growth. A base quarter after the last observed quarter T is
# taken from the same path; up to T, from the observed history.
#
# The formulas are set out in docs/methodology.md (section 9.4); code
# comments cite its equation labels (M30)-(M33).
#
# Inputs:  forecast paths (an array J x H x n in model units), the history,
#          a table of targets.
# Outputs: weights, effective sample size, weighted quantiles.


# ---- Quantities and moment conditions ----

#' Value of one target's quantity on every path.
#'
#' @param paths Array J x H x n of forecast paths in model units.
#' @param variable Short name of the variable.
#' @param h Forecast quarter (1 = T + 1).
#' @param measure "level", "yoy", "qoq_ann" or "annual_avg" (growth only
#'   for log series; checked when the targets are read).
#' @param history Data matrix in model units (quarters x n); its last rows
#'   provide base quarters at or before T.
#' @return Numeric vector of length J.
tilt_quantity <- function(paths, variable, h, measure, history) {
  y <- paths[, h, variable]
  if (measure == "level") return(y)
  if (measure == "annual_avg") return(annual_average_growth(paths, variable, h, history))
  lag <- switch(measure, yoy = 4L, qoq_ann = 1L, stop("unknown measure: ", measure))
  # Base quarter T + h - lag: from the path if it is a forecast quarter,
  # otherwise from the history (row T + h - lag, counting from its end).
  base <- if (h - lag >= 1) paths[, h - lag, variable] else history[nrow(history) + h - lag, variable]
  d <- (y - base) / 100                          # ln(x_t / x_{t-lag})
  if (measure == "yoy") 100 * (exp(d) - 1) else 100 * (exp(4 * d) - 1)
}

#' Growth of a four-quarter average over the previous four quarters, on
#' every path.
#'
#' annual_avg in quarter t = 100 (sum_{i=0..3} x_{t-i} / sum_{i=4..7} x_{t-i} - 1),
#' x = exp(y / 100). In a Q4 quarter this is calendar-year average growth.
#' The eight quarters t-7..t come from the path when they are forecast
#' quarters and from the observed history otherwise, as for yoy. The levels
#' are divided by exp(y_t / 100) before summing (the factor cancels in the
#' ratio), so large log levels cannot overflow.
#'
#' @param paths,variable,h,history As in tilt_quantity().
#' @return Numeric vector of length J.
annual_average_growth <- function(paths, variable, h, history) {
  J <- dim(paths)[1]
  if (nrow(history) + h - 7 < 1) {
    stop("annual_avg at forecast quarter ", h, " needs ", 8 - h, " observed quarters; the history has ",
         nrow(history))
  }
  # Step 1: the eight log levels y_{t-7}, ..., y_t of each path (J x 8).
  y8 <- matrix(NA_real_, J, 8)
  for (i in 0:7) {
    k <- h - i                                   # forecast quarter of t - i (k <= 0: observed)
    y8[, 8 - i] <- if (k >= 1) paths[, k, variable] else history[nrow(history) + k, variable]
  }
  # Step 2: levels relative to x_t, then the ratio of the two annual sums.
  x8 <- exp((y8 - y8[, 8]) / 100)
  100 * (rowSums(x8[, 5:8, drop = FALSE]) / rowSums(x8[, 1:4, drop = FALSE]) - 1)
}

#' Moment conditions g_j and gbar for a set of targets.
#'
#' Implements the targets of docs/methodology.md (M30).
#'
#' @param paths Array J x H x n of forecast paths.
#' @param targets data.frame with columns variable, h, measure, statistic
#'   ("median" or "mean") and target (x), one row per target
#'   (soft_targets() in R/scenarios.R).
#' @param history Data matrix (see tilt_quantity()).
#' @return List: g (J x m matrix), gbar (length m), quantity (J x m matrix
#'   of s(y), kept for reporting).
tilt_moments <- function(paths, targets, history) {
  m <- nrow(targets)
  J <- dim(paths)[1]
  s <- matrix(NA_real_, J, m)
  g <- matrix(NA_real_, J, m)
  gbar <- numeric(m)
  for (k in seq_len(m)) {
    s[, k] <- tilt_quantity(paths, targets$variable[k], targets$h[k], targets$measure[k], history)
    if (targets$statistic[k] == "median") {
      g[, k] <- as.numeric(s[, k] <= targets$target[k])
      gbar[k] <- 0.5
    } else if (targets$statistic[k] == "mean") {
      g[, k] <- s[, k]
      gbar[k] <- targets$target[k]
    } else {
      stop("unknown statistic: ", targets$statistic[k])
    }
  }
  list(g = g, gbar = gbar, quantity = s)
}


# ---- Weights ----

#' Entropic tilting weights: solve (16) by Newton's method.
#'
#' Works with f(gamma) = log sum_j exp(gamma' (g_j - gbar)), the log of the
#' objective of (16): the same minimizer, better scaled. With w(gamma) the
#' weights (15),
#'   gradient  = sum_j w_j (g_j - gbar)       (the miss of the targets),
#'   Hessian   = weighted covariance of the g_j under w,
#' so the gradient is zero exactly when the weighted moments meet the
#' targets. Each step solves Hessian * step = -gradient and halves the step
#' until f decreases. Every computation uses the log-sum-exp form, so large
#' gamma does not overflow.
#'
#' Implements docs/methodology.md (M30), (M31).
#'
#' @param g Matrix J x m of moment functions.
#' @param gbar Targets (length m).
#' @param tol Largest acceptable |weighted mean of g - gbar|.
#' @param maxit Largest number of Newton steps.
#' @return List: weights (length J, positive, sum 1), gamma, iterations,
#'   miss (weighted mean of g minus gbar), kl (sum w log(J w)). Stops if a
#'   target cannot be met.
tilt_weights <- function(g, gbar, tol = 1e-10, maxit = 100) {
  J <- nrow(g)
  m <- ncol(g)
  if (m == 0) return(list(weights = rep(1 / J, J), gamma = numeric(0), iterations = 0L,
                          miss = numeric(0), kl = 0))
  # Step 1: a necessary condition for each target on its own: gbar strictly
  # between the smallest and largest g_j (for a median: some paths on each
  # side of the value). Joint infeasibility shows up as non-convergence.
  for (k in seq_len(m)) {
    if (!(min(g[, k]) < gbar[k] && gbar[k] < max(g[, k]))) {
      stop("target ", k, " cannot be met by reweighting: the draws lie entirely on one side of it")
    }
  }
  dev <- sweep(g, 2, gbar)                       # g_j - gbar, J x m
  weights_at <- function(gamma) {
    a <- as.vector(dev %*% gamma)
    a_max <- max(a)
    e <- exp(a - a_max)
    list(w = e / sum(e), f = a_max + log(sum(e)))   # f = log sum exp(a)
  }

  # Step 2: Newton iterations from gamma = 0 (equal weights).
  gamma <- rep(0, m)
  cur <- weights_at(gamma)
  it <- 0L
  repeat {
    grad <- as.vector(crossprod(dev, cur$w))
    if (max(abs(grad)) < tol) break
    if (it >= maxit) stop("tilting did not converge in ", maxit, " steps (largest miss ",
                          signif(max(abs(grad)), 3), "): the targets may be jointly infeasible")
    hess <- crossprod(dev, dev * cur$w) - tcrossprod(grad)
    step <- tryCatch(solve(hess, -grad), error = function(e) {
      stop("tilting failed: the targets are redundant or conflict (singular Hessian)", call. = FALSE)
    })
    # Step 3: backtracking line search, so that every step lowers f.
    t <- 1
    repeat {
      cand <- weights_at(gamma + t * step)
      if (cand$f <= cur$f + 1e-4 * t * sum(grad * step) || t < 1e-10) break
      t <- t / 2
    }
    gamma <- gamma + t * step
    cur <- cand
    it <- it + 1L
  }
  w <- cur$w
  list(weights = w, gamma = gamma, iterations = it,
       miss = as.vector(crossprod(dev, w)), kl = sum(w[w > 0] * log(J * w[w > 0])))
}

#' Effective sample size of a set of weights (Kish): 1 / sum w_j^2.
#'
#' Equal weights give J; one weight of 1 gives 1. It says how many equally
#' weighted paths the tilted sample is worth.
#'
#' Implements docs/methodology.md (M33).
#'
#' @param w Weights (summing to 1).
#' @return Scalar.
tilt_ess <- function(w) 1 / sum(w^2)

#' Tilt a set of forecast paths to meet soft targets.
#'
#' Implements docs/methodology.md (M30)-(M31), (M33).
#'
#' @param paths Array J x H x n.
#' @param targets Target table (see tilt_moments()).
#' @param history Data matrix (see tilt_quantity()).
#' @param tol,maxit Solver settings (tilt_weights()).
#' @return List: weights, gamma, iterations, ess, kl, check (the targets
#'   with the achieved weighted statistic and the untilted one).
tilt_paths <- function(paths, targets, history, tol = 1e-10, maxit = 100) {
  mom <- tilt_moments(paths, targets, history)
  sol <- tilt_weights(mom$g, mom$gbar, tol, maxit)
  achieved <- untilted <- numeric(nrow(targets))
  for (k in seq_len(nrow(targets))) {
    if (targets$statistic[k] == "median") {
      # Share of weight at or below the target value (0.5 when met).
      achieved[k] <- sum(sol$weights * mom$g[, k])
      untilted[k] <- mean(mom$g[, k])
    } else {
      achieved[k] <- sum(sol$weights * mom$quantity[, k])
      untilted[k] <- mean(mom$quantity[, k])
    }
  }
  check <- data.frame(targets, goal = mom$gbar, achieved = achieved, untilted = untilted,
                      untilted_median = apply(mom$quantity, 2, stats::median))
  list(weights = sol$weights, gamma = sol$gamma, iterations = sol$iterations,
       ess = tilt_ess(sol$weights), kl = sol$kl, check = check)
}


# ---- Weighted quantiles ----

#' Weighted quantiles (weighted Hyndman-Fan type 7).
#'
#' R's default quantile (type 7) interpolates between two order statistics.
#' Akinshin (2023, "Weighted quantile estimators", Section 6) writes it as
#' sum_i W_i x_(i) with W_i = F(t_i) - F(t_{i-1}), t_i = i / n, and F the
#' CDF of the uniform distribution on [(h - 1)/n, h/n], h = (n - 1) p + 1.
#' The weighted version replaces t_i by the cumulative normalized weights of
#' the sorted values and n by Kish's effective sample size
#' n* = (sum w)^2 / sum w^2. With equal weights it equals quantile(x, p)
#' exactly; values with zero weight drop out.
#'
#' Implements docs/methodology.md (M32).
#'
#' @param x Numeric vector.
#' @param w Non-negative weights (NULL: equal), not all zero.
#' @param probs Probabilities.
#' @return Numeric vector, one value per probability.
weighted_quantiles <- function(x, w = NULL, probs) {
  if (is.null(w)) return(stats::quantile(x, probs, names = FALSE, type = 7))
  if (length(w) != length(x) || any(w < 0) || sum(w) <= 0) stop("weights must be non-negative, not all zero")
  o <- order(x)
  x <- x[o]
  w <- w[o] / sum(w[o])
  n_eff <- 1 / sum(w^2)
  t <- c(0, cumsum(w))
  t[length(t)] <- 1                               # guard against rounding in the sum
  vapply(probs, function(p) {
    h <- (n_eff - 1) * p + 1
    F <- pmin(pmax(t * n_eff - h + 1, 0), 1)      # CDF of U((h - 1)/n*, h/n*)
    sum(diff(F) * x)
  }, numeric(1))
}
