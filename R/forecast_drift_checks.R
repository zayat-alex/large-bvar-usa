# forecast_drift_checks.R: why are the forecasts low? Diagnostics of the
# long-run drift and of the starting point, for each COVID treatment.
#
# Why it exists: the baseline forecasts settle at low inflation and keep
# lowering interest rates. A VAR in levels with a constant and roots at or
# near 1 (the Minnesota and sum-of-coefficients priors pull towards a unit
# root in each equation; Crump et al. 2021, Section 2.2) lets every variable
# drift at a rate set by the constant, which reflects the trends of the
# estimation sample: since 1981, falling inflation and interest rates. The
# paper's own unconditional inflation forecast drifts down (its Figure 16),
# which its 2 percent anchor corrects (Section 4.2). These diagnostics
# measure the drift, compare forecasts made from recent origins with what
# happened, and show how much of today's forecast comes from unusual recent
# data (consumer sentiment far below its history, an oil price spike).
#
# Every diagnostic uses each posterior draw's no-shock forecast (the mean
# path given the draw's coefficients; forecast_unconditional(shocks = FALSE))
# and reports the median over draws, so simulated shocks add no noise.
# Measures: YoY growth (log approximation, y_h - y_{h-4}) for log series,
# levels for rates (growth_measure(), R/hyper_stability_checks.R).
#
# Inputs:  a run's draws and forecast history.
# Outputs: data.frames for reports/covid_comparison.qmd.


#' Median over draws of the no-shock forecast, in the annual measure.
#'
#' @param draws The run's draws (beta, sigma, lags).
#' @param history Forecast history (quarters x n, model units).
#' @param horizon Quarters to forecast.
#' @param is_log Logical vector, one per variable.
#' @return Matrix horizon x n: median YoY growth (log series) or level (rates).
no_shock_median <- function(draws, history, horizon, is_log) {
  p <- forecast_unconditional(draws$beta, draws$sigma, history, draws$lags, horizon, seed = 1, shocks = FALSE)
  apply(growth_measure(p, history, is_log, "annual"), c(2, 3), median)
}

#' Drift at long horizons: change per year of the no-shock median.
#'
#' @param med Output of no_shock_median() (at least `to` quarters).
#' @param variables Variables reported.
#' @param from,to Horizons between which the change is measured.
#' @return data.frame: variable, value at a few horizons, change per year
#'   between `from` and `to`.
drift_table <- function(med, variables, from = 60, to = 80) {
  hs <- intersect(c(4, 12, 20, 40, 60, 80), seq_len(nrow(med)))
  out <- data.frame(variable = variables, t(med[hs, variables, drop = FALSE]), row.names = NULL)
  names(out)[-1] <- paste0("h", hs)
  out$change_per_year <- (med[to, variables] - med[from, variables]) / ((to - from) / 4)
  out
}

#' Average change per year over 1973-2019, first to last value, in the same measure.
#'
#' @param y Data (quarters x n), rows named by quarter.
#' @param variables Variables.
#' @param is_log Logical vector, one per variable (named by variable).
#' @param end Last quarter.
#' @return Named numeric vector.
history_change_per_year <- function(y, variables, is_log, end = "2019Q4") {
  y <- y[seq_len(match(end, rownames(y))), , drop = FALSE]
  vapply(variables, function(v) {
    x <- if (is_log[[v]]) diff(y[, v], lag = 4) else y[, v]
    (x[length(x)] - x[1]) / (length(x) / 4)
  }, numeric(1))
}

#' Forecasts from earlier origins with the same estimates, against what happened.
#'
#' For each origin, the history is cut at the origin and each draw's
#' no-shock forecast is made with the run's coefficients. These are not
#' real-time forecasts: the coefficients were estimated on data that may
#' include the outcomes (in-sample for `none` and `exclude_window` after
#' 2022), and the data are the latest vintage.
#'
#' @param draws,history,is_log As in no_shock_median().
#' @param origins Origin quarters (labels).
#' @param horizons Horizons compared.
#' @param variables Variables reported.
#' @return data.frame: origin, horizon, variable, forecast, realized, error
#'   (realized minus forecast: positive = forecast too low).
origin_errors <- function(draws, history, is_log, origins, horizons, variables) {
  rows <- list()
  for (o in origins) {
    k <- match(o, rownames(history))
    h <- history[seq_len(k), , drop = FALSE]
    med <- no_shock_median(draws, h, max(horizons), is_log)
    for (hh in horizons) {
      if (k + hh > nrow(history)) next
      real <- vapply(variables, function(v)
        if (is_log[[v]]) history[k + hh, v] - history[k + hh - 4, v] else history[k + hh, v], numeric(1))
      rows[[length(rows) + 1]] <- data.frame(origin = o, horizon = hh, variable = variables,
                                             forecast = med[hh, variables], realized = real,
                                             row.names = NULL)
    }
  }
  out <- do.call(rbind, rows)
  out$error <- out$realized - out$forecast
  out
}

#' Effect of replacing recent data on the no-shock forecast.
#'
#' @param draws,history,is_log As in no_shock_median().
#' @param variable Variable whose recent values are replaced.
#' @param value Replacement value (model units).
#' @param n_quarters Number of last quarters replaced.
#' @param horizons,variables Reported cells.
#' @return data.frame: horizon, variable, effect (counterfactual minus actual).
counterfactual_effect <- function(draws, history, is_log, variable, value, n_quarters, horizons, variables) {
  base <- no_shock_median(draws, history, max(horizons), is_log)
  h <- history
  h[nrow(h) - seq_len(n_quarters) + 1, variable] <- value
  cf <- no_shock_median(draws, h, max(horizons), is_log)
  do.call(rbind, lapply(horizons, function(hh)
    data.frame(horizon = hh, variable = variables, effect = cf[hh, variables] - base[hh, variables],
               row.names = NULL)))
}
