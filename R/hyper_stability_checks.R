# hyper_stability_checks.R: do the hyperparameters change when the model is
# estimated on data that end earlier, and does it matter for the results?
#
# Why it exists: forecasts from past origins (comparison with the Survey of
# Professional Forecasters, historical scenario replays) need the model
# estimated on the data available at each origin. Re-estimating the
# hyperparameters (lambda, mu, psi) at every origin costs a mode search and
# a Metropolis run each; Crump et al. (2021, Section 2.4) do so. If the
# hyperparameters estimated on earlier samples are close to the current ones,
# and fixing them leaves forecasts and impulse responses unchanged, the
# current hyperparameters can be reused at past origins and only the
# coefficients (B, Sigma) redrawn, which is cheap. This file supplies the
# pieces of that check; scripts/check_hyper_stability.R runs it.
#
# Design (one origin):
#   1. cut the estimation sample at the origin quarter;
#   2. find the mode of the hyperparameter posterior on that sample, and the
#      Hessian there (docs/methodology.md section 7, steps 1-2; no MCMC);
#   3. compare the modes: theta = log(h), so a difference is a log ratio,
#      read in units of the posterior standard deviation of theta (the
#      current run's MCMC draws, and the normal approximation from the
#      Hessian at each origin);
#   4. on the same sample, draw coefficients given (i) the origin's mode and
#      (ii) the current mode, with the same random numbers, and compare the
#      forecasts and impulse responses they give. A second draw under (ii)
#      with other random numbers measures how much the comparison moves by
#      simulation noise alone.
#
# Inputs:  the estimation data (COVID method none), the current run's
#          estimation.rds and draws.rds, the run configuration.
# Outputs: data.frames and figures used by reports/hyper_stability.qmd.


# ---- Sample ----

#' Estimation data up to and including an origin quarter.
#'
#' @param est_data data.frame with `quarter` (Date) and the model variables.
#' @param origin Quarter label, e.g. "2016Q2": the last quarter of the sample.
#' @return est_data restricted to quarters up to the origin.
sample_to_origin <- function(est_data, origin) {
  q <- parse_quarter(origin)
  if (!q %in% est_data$quarter) stop("origin ", origin, " is not a quarter of the data")
  est_data[est_data$quarter <= q, , drop = FALSE]
}


# ---- Mode search with the starts in parallel ----

#' find_hyper_mode() with each starting point in its own R process.
#'
#' The starts are independent optimizations, so running them side by side
#' gives the same result as find_hyper_mode() (each start is deterministic)
#' in the time of the slowest start.
#'
#' @param setup Output of estimation_setup().
#' @param start_factors Starting points (config mcmc$mode_search$start_factors).
#' @param maxit Iteration limit of each optimization.
#' @param r_dir Folder of the project's R files (sourced in each process).
#' @return As find_hyper_mode(): table (one row per start), theta_mode, best;
#'   plus seconds (wall-clock time).
find_hyper_mode_parallel <- function(setup, start_factors, maxit, r_dir = "R") {
  t0 <- proc.time()[["elapsed"]]
  cl <- parallel::makePSOCKcluster(length(start_factors))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  # Each process uses this project's package library and functions.
  r_files <- normalizePath(list.files(r_dir, pattern = "\\.R$", full.names = TRUE))
  parallel::clusterCall(cl, function(lib, files) {
    .libPaths(lib)
    for (f in files) source(f)
    NULL
  }, .libPaths(), r_files)
  runs <- parallel::clusterApply(cl, start_factors, function(f, setup, maxit)
    find_hyper_mode(setup, f, maxit), setup = setup, maxit = maxit)

  table <- do.call(rbind, lapply(runs, `[[`, "table"))
  best  <- which.max(table$log_post_theta)
  list(table = table, theta_mode = runs[[best]]$theta_mode, best = best,
       seconds = proc.time()[["elapsed"]] - t0)
}


# ---- Comparing the modes ----

#' Hyperparameter modes at an origin against the current ones.
#'
#' All differences are on theta = log(h), the scale of the mode search and
#' the sampler: a difference of 0.1 means h is about 10 percent higher.
#' The shift is read in two yardsticks: the current posterior s.d. of theta
#' (from the 1,000 retained MCMC draws, the most accurate measure), and the
#' origin's approximate s.d. from the Hessian (the normal approximation;
#' a shorter sample gives a wider posterior). A shift well below 1 s.d. is
#' within the current posterior's uncertainty; several s.d. means the
#' posterior is centred elsewhere.
#'
#' @param mode_past,proposal_past find_hyper_mode() and hyper_proposal()
#'   results at the origin.
#' @param mode_now,proposal_now The same for the current run (estimation.rds).
#' @param hyper_now Current hyperparameter draws (draws.rds$hyper; original scale).
#' @param setup_past estimation_setup() at the origin (for the bounds).
#' @param variables Variable names, in model order (to label psi_j).
#' @return data.frame, one row per hyperparameter.
hyper_shift_table <- function(mode_past, proposal_past, mode_now, proposal_now, hyper_now,
                              setup_past, variables) {
  nm <- names(mode_now$theta_mode)
  theta_past <- mode_past$theta_mode[nm]
  theta_now  <- mode_now$theta_mode[nm]
  sd_now_mcmc <- apply(log(hyper_now[, nm, drop = FALSE]), 2, sd)
  shift <- theta_past - theta_now
  label <- nm
  is_psi <- grepl("^psi[0-9]+$", nm)
  label[is_psi] <- paste0("psi_", variables[as.integer(sub("psi", "", nm[is_psi]))])
  # Within 1 percent of a search bound, on the original scale.
  h_past <- exp(theta_past)
  near_bound <- h_past < setup_past$lower[nm] * 1.01 | h_past > setup_past$upper[nm] / 1.01
  data.frame(hyperparameter = label,
             mode_now = exp(theta_now), mode_past = h_past,
             log_ratio = shift,
             sd_now_mcmc = sd_now_mcmc, sd_now_hessian = proposal_now$approx_sd[nm],
             sd_past_hessian = proposal_past$approx_sd[nm],
             shift_in_sd_now = shift / sd_now_mcmc,
             shift_in_sd_past = shift / proposal_past$approx_sd[nm],
             near_bound_past = near_bound,
             row.names = NULL)
}


# ---- Coefficient draws at fixed hyperparameters ----

#' Draws of (B, Sigma) given one fixed vector of hyperparameters.
#'
#' With h fixed, the posterior of (B, Sigma) is normal-inverse-Wishart in
#' closed form (docs/methodology.md (M16), (M17)), so no sampler is needed.
#'
#' @param theta Named hyperparameters on the log scale (lambda, soc, psi1..psin).
#' @param setup estimation_setup() of the sample to draw on.
#' @param n_draws Number of draws.
#' @param seed Seed, set before the first draw.
#' @return List: beta (draws x K x n), sigma (draws x n x n).
draws_at_mode <- function(theta, setup, n_draws, seed) {
  h <- complete_hyper(exp(theta), setup)
  post <- niw_posterior(h, setup)
  set.seed(seed)
  draw_coefficients(post, n_draws)
}


# ---- Comparing forecasts and impulse responses ----

#' Compare two sets of draws by their median and 90 percent band.
#'
#' For each variable and horizon: the 5, 50 and 95 percent quantiles under
#' the reference (`ref`, the current hyperparameters) and the alternative
#' (`alt`, the origin's hyperparameters); the difference of medians as a
#' share of the reference's 90 percent band width; the ratio of band widths.
#' The same two measures for `ref_seed` (the reference with other random
#' numbers) against `ref` give the level reached by simulation noise alone.
#'
#' @param ref,alt,ref_seed Arrays draws x horizons x n (forecasts or IRFs).
#' @param variables Variables to compare.
#' @param horizons Horizon labels (one per position of the second dimension).
#' @param base Optional named vector subtracted per variable (the last
#'   observation, so forecasts read as change since the origin).
#' @return data.frame, one row per variable and horizon.
compare_bands <- function(ref, alt, ref_seed, variables, horizons, base = NULL) {
  rows <- list()
  for (v in variables) {
    j <- match(v, dimnames(ref)[[3]])
    b <- if (is.null(base)) 0 else base[[v]]
    for (k in seq_along(horizons)) {
      q <- lapply(list(ref = ref, alt = alt, seed = ref_seed),
                  function(a) quantile(a[, k, j] - b, c(0.05, 0.5, 0.95), names = FALSE))
      width <- q$ref[3] - q$ref[1]
      rows[[length(rows) + 1]] <- data.frame(
        variable = v, horizon = horizons[k],
        q05_ref = q$ref[1], q50_ref = q$ref[2], q95_ref = q$ref[3],
        q05_alt = q$alt[1], q50_alt = q$alt[2], q95_alt = q$alt[3],
        median_diff_share = if (width > 0) (q$alt[2] - q$ref[2]) / width else NA_real_,
        width_ratio = if (width > 0) (q$alt[3] - q$alt[1]) / width else NA_real_,
        median_diff_share_seed = if (width > 0) (q$seed[2] - q$ref[2]) / width else NA_real_,
        width_ratio_seed = if (width > 0) (q$seed[3] - q$seed[1]) / width else NA_real_)
    }
  }
  do.call(rbind, rows)
}

#' Summarize compare_bands() over windows of horizons, per variable.
#'
#' @param cmp Output of compare_bands() (numeric horizons).
#' @param windows Named list of horizon ranges, e.g. list(`1-12` = 1:12).
#' @return data.frame: variable, window, largest |median difference| (share
#'   of the band width), range of width ratios, and the same for the
#'   seed-noise comparison. Cells with a zero-width band (slow variables on
#'   impact in the recursive IRF, zero by construction) are skipped.
summarize_bands <- function(cmp, windows) {
  rows <- list()
  for (v in unique(cmp$variable)) for (w in names(windows)) {
    d <- cmp[cmp$variable == v & cmp$horizon %in% windows[[w]] & !is.na(cmp$median_diff_share), ]
    if (!nrow(d)) next
    rows[[length(rows) + 1]] <- data.frame(
      variable = v, horizons = w,
      max_abs_median_diff = max(abs(d$median_diff_share)),
      width_ratio_min = min(d$width_ratio), width_ratio_max = max(d$width_ratio),
      seed_max_abs_median_diff = max(abs(d$median_diff_share_seed)),
      seed_width_ratio_min = min(d$width_ratio_seed), seed_width_ratio_max = max(d$width_ratio_seed))
  }
  do.call(rbind, rows)
}


# ---- Comparing in growth rates and percentage points ----
#
# Bands of a BVAR are wide, so a change that is small against the band can
# still be large in economic terms (half a point of GDP growth). The
# measures below read the same draws in the units a forecast is quoted in,
# path by path (the median of a growth rate is not the growth rate of the
# median level), and scale each change by the variable's own history.
#
# Notation: model units are 100 x log for log series (config/series.csv,
# transform "log") and the raw value for rates and indexes ("raw").
#   annual measure:    log series: YoY growth, y_h - y_{h-4} (percent, log
#                      approximation); raw series: the level (pp or index points)
#   quarterly measure: log series only: QoQ annualized growth, 4 (y_h - y_{h-1})
# Yardstick from history: the same annual measure as a change, i.e. YoY
# growth for log series and the 4-quarter change for raw series; its s.d.
# and its 10-90 percentile range.

#' Forecast paths in the annual or quarterly measure.
#'
#' @param paths Array draws x H x n, horizons 1..H, model units.
#' @param history Matrix of observed data (quarters x n); its last 4 rows
#'   give the quarters before the forecast.
#' @param is_log Logical vector (one per variable): log series?
#' @param measure "annual" or "quarterly".
#' @return Array draws x H x n in the measure (NA for raw series in "quarterly").
growth_measure <- function(paths, history, is_log, measure = c("annual", "quarterly")) {
  measure <- match.arg(measure)
  H <- dim(paths)[2]
  lag <- if (measure == "annual") 4 else 1
  out <- paths
  for (j in seq_len(dim(paths)[3])) {
    if (!is_log[j]) {
      if (measure == "quarterly") out[, , j] <- NA_real_
      next
    }
    hist_j <- history[nrow(history) - (lag - 1):0, j]          # the lag quarters before h = 1
    p <- matrix(paths[, , j], nrow = dim(paths)[1])   # draws x H, also for one draw
    ext <- cbind(matrix(hist_j, nrow(p), lag, byrow = TRUE), p)
    d <- ext[, lag + seq_len(H), drop = FALSE] - ext[, seq_len(H), drop = FALSE]
    out[, , j] <- if (measure == "quarterly") 4 * d else d
  }
  out
}

#' Impulse responses in the annual measure.
#'
#' A response r_h (h = 0..H) is the difference from the no-shock path, so the
#' effect on YoY growth of a log series is r_h - r_{h-4}, with r = 0 before
#' impact. Raw series keep their level response.
#'
#' @param irf Array draws x (H + 1) x n.
#' @param is_log Logical vector, one per variable.
#' @return Array of the same shape.
irf_annual_measure <- function(irf, is_log) {
  out <- irf
  n_h <- dim(irf)[2]
  for (j in which(is_log)) {
    r <- matrix(irf[, , j], nrow = dim(irf)[1])   # draws x horizons, also for one draw
    prev <- cbind(matrix(0, nrow(r), 4), r)[, seq_len(n_h), drop = FALSE]
    out[, , j] <- r - prev
  }
  out
}

#' Typical size of each variable's annual change in history.
#'
#' @param y Data matrix (quarters x n, model units), rows named "1973Q1", ...
#' @param is_log Logical vector, one per variable.
#' @param end Last quarter used (label); the default leaves out 2020 onward,
#'   whose extreme quarters would dominate the s.d.
#' @return data.frame: variable, measure ("YoY growth, %" or "4-quarter
#'   change"), sd, p10, p90, range_10_90, and qoq_sd, qoq_range_10_90 (QoQ
#'   annualized growth, log series), over the quarters up to `end`.
historical_scale <- function(y, is_log, end = "2019Q4") {
  keep <- seq_len(match(end, rownames(y)))
  y <- y[keep, , drop = FALSE]
  rows <- lapply(seq_len(ncol(y)), function(j) {
    d4 <- diff(y[, j], lag = 4)
    q <- quantile(d4, c(0.1, 0.9), names = FALSE)
    q1 <- 4 * diff(y[, j])
    qq <- quantile(q1, c(0.1, 0.9), names = FALSE)
    data.frame(variable = colnames(y)[j],
               measure = if (is_log[j]) "YoY growth, %" else "4-quarter change",
               sd = sd(d4), p10 = q[1], p90 = q[2], range_10_90 = q[2] - q[1],
               qoq_sd = if (is_log[j]) sd(q1) else NA_real_,
               qoq_range_10_90 = if (is_log[j]) qq[2] - qq[1] else NA_real_)
  })
  do.call(rbind, rows)
}

#' Difference of medians between two sets of paths, and the noise level.
#'
#' @param ref,alt,ref_seed Arrays draws x horizons x n in the same measure.
#' @param variables Variables to compare.
#' @param horizons Horizon labels.
#' @return data.frame: variable, horizon, median_ref, median_alt,
#'   diff (alt - ref), diff_seed (ref_seed - ref), in the measure's units.
compare_medians <- function(ref, alt, ref_seed, variables, horizons) {
  rows <- list()
  for (v in variables) {
    j <- match(v, dimnames(ref)[[3]])
    for (k in seq_along(horizons)) {
      m <- vapply(list(ref, alt, ref_seed), function(a) median(a[, k, j]), numeric(1))
      rows[[length(rows) + 1]] <- data.frame(variable = v, horizon = horizons[k],
                                             median_ref = m[1], median_alt = m[2],
                                             diff = m[2] - m[1], diff_seed = m[3] - m[1])
    }
  }
  out <- do.call(rbind, rows)
  out[!is.na(out$diff), , drop = FALSE]
}

#' Largest change per variable and window, in units and against history.
#'
#' @param cmp Output of compare_medians().
#' @param scale Output of historical_scale().
#' @param windows Named list of horizon ranges.
#' @param measure "annual" (yardstick columns sd, range_10_90) or "quarterly"
#'   (qoq_sd, qoq_range_10_90).
#' @return data.frame: variable, horizons, max_abs_diff (units),
#'   max_abs_noise (units), diff_in_sd (max_abs_diff / historical s.d.),
#'   diff_in_range (max_abs_diff / historical 10-90 range), noise_in_sd.
summarize_medians <- function(cmp, scale, windows, measure = c("annual", "quarterly")) {
  measure <- match.arg(measure)
  sd_column    <- if (measure == "annual") "sd" else "qoq_sd"
  range_column <- if (measure == "annual") "range_10_90" else "qoq_range_10_90"
  rows <- list()
  for (v in unique(cmp$variable)) for (w in names(windows)) {
    d <- cmp[cmp$variable == v & cmp$horizon %in% windows[[w]], ]
    if (!nrow(d)) next
    s <- scale[scale$variable == v, ]
    mx <- max(abs(d$diff)); nz <- max(abs(d$diff_seed))
    rows[[length(rows) + 1]] <- data.frame(
      variable = v, horizons = w, max_abs_diff = mx, max_abs_noise = nz,
      diff_in_sd = mx / s[[sd_column]], diff_in_range = mx / s[[range_column]],
      noise_in_sd = nz / s[[sd_column]])
  }
  do.call(rbind, rows)
}


# ---- Plots (base graphics; colours of R/estimate_checks.R) ----

#' Shift of each hyperparameter's mode, in current posterior s.d.
#'
#' @param shift Output of hyper_shift_table().
#' @param origin Origin label (title).
plot_hyper_shift <- function(shift, origin) {
  n <- nrow(shift)
  x <- shift$shift_in_sd_now
  lim <- max(3, abs(x), na.rm = TRUE) * c(-1, 1)
  par(mar = c(4.5, 9, 2.5, 1), col.axis = ink[["secondary"]], col.lab = ink[["primary"]])
  plot(x, rev(seq_len(n)), xlim = lim, yaxt = "n", pch = 19, col = check_colours[["blue"]],
       xlab = "mode at origin minus current mode, in current posterior s.d. of log h",
       ylab = "", bty = "n", main = paste("Hyperparameter modes, origin", origin),
       font.main = 1, col.main = ink[["primary"]])
  abline(v = c(-2, -1, 1, 2), col = ink[["grid"]], lty = c(2, 1, 1, 2))
  abline(v = 0, col = ink[["secondary"]])
  axis(2, at = rev(seq_len(n)), labels = shift$hyperparameter, las = 1, cex.axis = 0.7)
}

#' Median and 68/90 percent bands under the two sets of hyperparameters.
#'
#' Current hyperparameters in blue (shaded bands), the origin's in orange
#' (lines), for one variable; x = horizon. Forecasts start at the last
#' observation (change 0), so the lines join the history at the origin.
#'
#' @param ref,alt Arrays draws x horizons x n.
#' @param variable Variable name.
#' @param x Horizon values (one per position of the second dimension).
#' @param base Value subtracted (last observation), or 0.
#' @param title Panel title.
#' @param start_at_zero TRUE for forecasts: prepend horizon 0 with value 0.
plot_band_pair <- function(ref, alt, variable, x, base, title, start_at_zero) {
  j <- match(variable, dimnames(ref)[[3]])
  qs <- lapply(list(ref, alt), function(a)
    t(apply(a[, , j, drop = FALSE] - base, 2, quantile, c(0.05, 0.16, 0.5, 0.84, 0.95))))
  if (start_at_zero) {
    x <- c(0, x)
    qs <- lapply(qs, function(q) rbind(rep(0, 5), q))
  }
  ylim <- range(unlist(qs))
  plot(NA, xlim = range(x), ylim = ylim, bty = "n", xlab = "quarters", ylab = "",
       main = title, font.main = 1, cex.main = 0.9, col.main = ink[["primary"]],
       col.axis = ink[["secondary"]])
  abline(h = 0, col = ink[["grid"]])
  q <- qs[[1]]
  polygon(c(x, rev(x)), c(q[, 1], rev(q[, 5])), col = "#2a78d633", border = NA)
  polygon(c(x, rev(x)), c(q[, 2], rev(q[, 4])), col = "#2a78d666", border = NA)
  lines(x, q[, 3], col = check_colours[["blue"]], lwd = 2)
  q <- qs[[2]]
  lines(x, q[, 3], col = check_colours[["orange"]], lwd = 2)
  for (k in c(1, 2, 4, 5)) lines(x, q[, k], col = check_colours[["orange"]], lwd = 1,
                                 lty = if (k %in% c(1, 5)) 3 else 2)
}
