# recovery_checks.R: simulation recovery. Simulate data from a known VAR,
# estimate the model on them, and check how often the posterior bands
# contain the truth.
#
# Why it exists: the estimator, the forecast simulator and the impulse
# responses are each checked against exact results elsewhere
# (reports/forecast_checks.qmd and the tests). This check runs the whole chain together, from
# data to hyperparameters, coefficient draws, forecast bands and impulse
# responses, in a setting where the right answer is known. It catches errors
# that pass the separate checks (draws paired with the wrong covariance, a
# scaling slip, bands that are too narrow), and it shows how strongly the
# prior shrinks the estimates.
#
# Design (one replication):
#   1. the true model is a fixed VAR(p), coefficients B and covariance Sigma;
#   2. simulate T0 + H quarters from it, starting from p given quarters: the
#      p starting quarters and the first T0 - p simulated quarters are the
#      sample (T0 rows, as in the data), the last H are the true future;
#   3. estimate the BVAR on the sample (R/estimate_bvar.R);
#   4. coverage: is the truth inside the 68% and 90% posterior bands (16-84
#      and 5-95 percent quantiles)? For forecasts 1..H quarters ahead (the
#      predictive density, with shocks, against the simulated future), for
#      the recursive impulse response (against the true model's response),
#      and for the coefficients and Sigma (against B and Sigma).
#
# How to read coverage: the bands are Bayesian credible bands, not
# confidence intervals, so nominal coverage over repeated samples from one
# fixed model is not guaranteed. For the predictive density it should be
# close; for coefficients the prior shrinks cross-variable lags towards zero,
# so some undercoverage is expected there.
#
# Inputs:  the true (B, Sigma), starting quarters, an estimation setup's
#          series and prior settings, MCMC settings, seeds.
# Outputs: per replication, coverage flags and summaries; across
#          replications, coverage rates, shrinkage bias and figures. Run by
#          scripts/check_recovery.R, shown in reports/recovery_checks.qmd.


# ---- Simulation ----

#' Simulate a path from a known VAR.
#'
#' Uses forecast_unconditional() (R/forecast.R; docs/methodology.md (M18))
#' with the single "draw" (B, Sigma): the same recursion and the same
#' shocks e_t = L u_t, so the simulated data come from exactly the model the
#' forecasts assume.
#'
#' @param B Coefficients K x n (constant in row 1).
#' @param Sigma Covariance n x n.
#' @param init Starting quarters, a matrix with at least p rows and the
#'   variable names as columns.
#' @param lags Number of lags p.
#' @param n_quarters Number of quarters to simulate after `init`.
#' @param seed Seed of the shocks.
#' @return Matrix n_quarters x n.
simulate_var <- function(B, Sigma, init, lags, n_quarters, seed) {
  sim <- forecast_unconditional(array(B, c(1, dim(B))), array(Sigma, c(1, dim(Sigma))), init, lags,
                                n_quarters, seed = seed)
  sim[1, , ]
}


# ---- Coverage ----

#' Is the truth inside the 68% and 90% bands of the draws?
#'
#' @param draws Array with draws on the first dimension.
#' @param truth Array of the remaining dimensions (or a vector of their
#'   length).
#' @return List of two logical arrays shaped like `truth`: in68 (between
#'   the 16% and 84% quantiles) and in90 (between 5% and 95%).
band_coverage <- function(draws, truth) {
  d <- dim(draws)[-1]
  q <- apply(draws, seq_along(dim(draws))[-1], stats::quantile, c(0.05, 0.16, 0.84, 0.95), names = FALSE)
  q <- matrix(q, 4)                     # rows: 5, 16, 84, 95 percent; columns: cells
  lo90 <- array(q[1, ], d); lo68 <- array(q[2, ], d)
  hi68 <- array(q[3, ], d); hi90 <- array(q[4, ], d)
  truth <- array(truth, d)
  list(in68 = truth >= lo68 & truth <= hi68, in90 = truth >= lo90 & truth <= hi90)
}

#' Group of each coefficient of B (K x n), for summaries.
#'
#' Row 1 is the constant; row 1 + (s - 1) n + i holds lag s of variable i
#' (R/forecast.R). Groups: constant; own first lag (lag 1 of the equation's
#' own variable, prior mean 1 or 0); own longer lags (lags 2..p, prior mean
#' 0); cross lags (other variables, any lag, prior mean 0).
#'
#' @param n,lags Number of variables and lags.
#' @return Character matrix K x n.
coefficient_groups <- function(n, lags) {
  g <- matrix("cross lags", 1 + n * lags, n)
  g[1, ] <- "constant"
  for (s in seq_len(lags)) {
    for (j in seq_len(n)) g[1 + (s - 1) * n + j, j] <- if (s == 1) "own first lag" else "own longer lags"
  }
  g
}


# ---- One replication ----

#' Simulate, estimate and record coverage, once.
#'
#' @param rep Replication number (its seeds are seed + rep and derived).
#' @param truth List: B, Sigma, init (p starting quarters), irf (true
#'   recursive response, (H_irf + 1) x n).
#' @param series Series configuration of the model's variables.
#' @param prior_cfg Prior settings (config `prior`).
#' @param lags Number of lags p.
#' @param n_sample Number of quarters in the sample, including init.
#' @param horizon Forecast horizon H.
#' @param mcmc List n_draw, n_burn, n_thin.
#' @param recursive List slow, policy, fast (the recursive identification).
#' @param irf_horizon Horizon of the impulse responses.
#' @param seed Base seed.
#' @return List: rep, ok, error (if the fit failed), forecast coverage
#'   (H x n, in68, in90), irf coverage ((H_irf + 1) x n), coefficient
#'   coverage (K x n), Sigma coverage (n x n), posterior means of B and
#'   Sigma, acceptance, largest root of the posterior mean, the sample's
#'   last quarters, the true future and forecast quantiles (for figures),
#'   seconds.
recovery_replication <- function(rep, truth, series, prior_cfg, lags, n_sample, horizon, mcmc,
                                 recursive, irf_horizon, seed) {
  t0 <- proc.time()[["elapsed"]]
  vars <- series$short_name
  s_rep <- seed + rep * 10   # simulation seed + 0, MCMC + 1, forecasts + 2

  # Step 1: simulate the sample and the true future.
  init <- truth$init
  sim <- simulate_var(truth$B, truth$Sigma, init, lags, n_sample - nrow(init) + horizon, seed = s_rep)
  y <- rbind(init, sim[seq_len(n_sample - nrow(init)), , drop = FALSE])
  rownames(y) <- NULL
  future <- sim[n_sample - nrow(init) + seq_len(horizon), , drop = FALSE]

  # Step 2: estimate. A failed fit is recorded, not dropped.
  setup <- estimation_setup(y, series, prior_cfg, lags)
  fit <- tryCatch(estimate_bvar(setup, mcmc$n_draw, mcmc$n_burn, mcmc$n_thin, seed = s_rep + 1),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(rep = rep, ok = FALSE, error = conditionMessage(fit),
                seconds = proc.time()[["elapsed"]] - t0))
  }

  # Step 3: forecasts (predictive density, with shocks) against the future.
  fc <- forecast_unconditional(fit$beta, fit$sigma, y, lags, horizon, seed = s_rep + 2)
  cov_fc <- band_coverage(fc, future)

  # Step 4: recursive impulse response against the true response.
  ir <- irf_draws(fit$beta, fit$sigma, lags, irf_horizon, vars, character(0), recursive)[[1]]
  cov_irf <- band_coverage(ir, truth$irf)

  # Step 5: coefficients and Sigma against the truth.
  cov_B <- band_coverage(fit$beta, truth$B)
  cov_S <- band_coverage(fit$sigma, truth$Sigma)
  B_mean <- apply(fit$beta, c(2, 3), mean)

  list(rep = rep, ok = TRUE,
       forecast = cov_fc, irf = cov_irf, B = cov_B, Sigma = cov_S,
       B_mean = B_mean, Sigma_mean = apply(fit$sigma, c(2, 3), mean),
       acceptance = fit$sampler$acceptance, lambda_median = stats::median(fit$hyper[, "lambda"]),
       max_root_mean = companion_max_modulus(B_mean, lags),
       last_sample = y[nrow(y) - 7:0, , drop = FALSE], future = future,
       fc_quantiles = apply(fc, c(2, 3), stats::quantile, c(0.05, 0.16, 0.5, 0.84, 0.95), names = FALSE),
       seconds = proc.time()[["elapsed"]] - t0)
}


# ---- Summaries across replications ----

#' Binomial standard error of a coverage rate.
#'
#' @param p Nominal rate.
#' @param n Number of replications.
coverage_se <- function(p, n) sqrt(p * (1 - p) / n)

#' Coverage rates of one object, per cell, across replications.
#'
#' @param reps List of successful recovery_replication() results.
#' @param what "forecast", "irf", "B" or "Sigma".
#' @return List: c68, c90 (arrays shaped like the object: share of
#'   replications whose band contains the truth).
coverage_rates <- function(reps, what) {
  avg <- function(band) Reduce(`+`, lapply(reps, function(r) r[[what]][[band]] * 1)) / length(reps)
  list(c68 = avg("in68"), c90 = avg("in90"))
}

#' Coverage of the coefficients and Sigma by group, and shrinkage bias.
#'
#' Bias is the average over replications of (posterior mean - truth); for
#' the lag groups it is also shown relative to the average absolute true
#' coefficient of the group, to show how far the prior pulls them to zero.
#'
#' @param reps List of successful recovery_replication() results.
#' @param truth List with B and Sigma.
#' @param lags Number of lags p.
#' @return data.frame: group, cells, coverage 68 and 90, mean bias, RMSE of
#'   the posterior mean, mean |truth|.
parameter_coverage_table <- function(reps, truth, lags) {
  n <- ncol(truth$B)
  groups <- coefficient_groups(n, lags)
  cB <- coverage_rates(reps, "B")
  cS <- coverage_rates(reps, "Sigma")
  errB <- lapply(reps, function(r) r$B_mean - truth$B)
  errS <- lapply(reps, function(r) r$Sigma_mean - truth$Sigma)
  row <- function(name, sel, c68, c90, err, true_vals) {
    e <- sapply(err, function(x) x[sel])            # cells x reps
    data.frame(group = name, cells = sum(sel), cover68 = mean(c68[sel]), cover90 = mean(c90[sel]),
               mean_bias = mean(e), rmse = sqrt(mean(e^2)), mean_abs_truth = mean(abs(true_vals[sel])))
  }
  out <- lapply(c("constant", "own first lag", "own longer lags", "cross lags"), function(g) {
    row(g, groups == g, cB$c68, cB$c90, errB, truth$B)
  })
  diag_sel <- diag(n) == 1
  off_sel <- lower.tri(truth$Sigma)
  out <- c(out, list(row("Sigma, variances", diag_sel, cS$c68, cS$c90, errS, truth$Sigma),
                     row("Sigma, covariances", off_sel, cS$c68, cS$c90, errS, truth$Sigma)))
  do.call(rbind, out)
}


# ---- Figures ----

#' Coverage by horizon, one panel per variable, with nominal levels and a
#' band of two binomial standard errors around them.
#'
#' @param rates Output of coverage_rates() for "forecast" or "irf" (horizon
#'   x variable).
#' @param n_reps Number of replications.
#' @param x Horizons for the x axis.
#' @param xlab Axis label.
#' @param variables Variable names (panel titles).
plot_coverage <- function(rates, n_reps, x, xlab, variables) {
  # Bottom outer margin keeps room for legend_coverage() below the panels.
  op <- graphics::par(mfrow = c(2, ceiling(length(variables) / 2)), mar = c(3.5, 3, 2, 0.5),
                      mgp = c(2, 0.6, 0), oma = c(2, 0, 0, 0))
  on.exit(graphics::par(op))
  for (j in seq_along(variables)) {
    graphics::plot(NA, xlim = range(x), ylim = c(0, 1), xlab = xlab, ylab = "coverage", bty = "n",
                   main = variables[j], font.main = 1, col.axis = ink[["secondary"]])
    for (p in c(0.68, 0.90)) {
      se <- coverage_se(p, n_reps)
      graphics::rect(min(x), p - 2 * se, max(x), p + 2 * se, col = "#0000000f", border = NA)
      graphics::abline(h = p, col = ink[["grid"]], lty = 2)
    }
    graphics::lines(x, rates$c90[, j], col = check_colours[["blue"]], lwd = 2, type = "b", pch = 16, cex = 0.6)
    graphics::lines(x, rates$c68[, j], col = check_colours[["orange"]], lwd = 2, type = "b", pch = 16, cex = 0.6)
  }
}

#' Legend for plot_coverage(), drawn across the bottom of the whole figure.
legend_coverage <- function() {
  op <- graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  on.exit(graphics::par(op))
  graphics::plot.new()
  graphics::legend("bottom", horiz = TRUE, bty = "n", xpd = NA, cex = 0.9,
                   legend = c("90% band", "68% band", "nominal +/- 2 s.e."),
                   col = c(check_colours[["blue"]], check_colours[["orange"]], "#00000030"),
                   lwd = c(2, 2, 8))
}

#' Forecast fan of one replication against its true future, one panel per
#' variable; the forecast starts at the last sample quarter.
#'
#' @param r One recovery_replication() result.
#' @param variables Variable names.
plot_recovery_fan <- function(r, variables) {
  op <- graphics::par(mfrow = c(2, ceiling(length(variables) / 2)), mar = c(3.5, 3, 2, 0.5),
                      mgp = c(2, 0.6, 0))
  on.exit(graphics::par(op))
  k <- nrow(r$last_sample); H <- nrow(r$future)
  for (j in seq_along(variables)) {
    last <- r$last_sample[k, j]
    q <- cbind(last, r$fc_quantiles[, , j])                  # 5 x (H + 1), joined to the last quarter
    x <- 0:H; xh <- seq(-k + 1, 0)
    yl <- range(c(q, r$last_sample[, j], r$future[, j]))
    graphics::plot(xh, r$last_sample[, j], type = "l", xlim = c(-k + 1, H), ylim = yl, bty = "n",
                   lwd = 1.5, col = ink[["primary"]], xlab = "quarters from the end of the sample", ylab = "",
                   main = variables[j], font.main = 1, col.axis = ink[["secondary"]])
    graphics::polygon(c(x, rev(x)), c(q[1, ], rev(q[5, ])), col = "#2a78d633", border = NA)
    graphics::polygon(c(x, rev(x)), c(q[2, ], rev(q[4, ])), col = "#2a78d666", border = NA)
    graphics::lines(x, q[3, ], col = check_colours[["blue"]], lwd = 2)
    graphics::lines(x, c(last, r$future[, j]), col = ink[["primary"]], lwd = 1.5, lty = 2)
    graphics::abline(v = 0, col = ink[["grid"]])
  }
}
