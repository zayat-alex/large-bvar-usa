# forecast_review.R: is the baseline forecast reasonable? (spec 9.4)
#
# Why it exists: the engine's forecasts are checked for correctness against
# exact results elsewhere (reports/forecast_checks.qmd). This review asks
# whether they are economically sensible: against the latest Survey of
# Professional Forecasters median, against the ranges the variables have
# taken in history, at long horizons (where a model in levels can drift),
# with and without the inflation anchor, under each COVID treatment, and
# how much of what is published is simulation noise.
#
# Notation: model units are 100 x log for log series (config/series.csv,
# transform "log") and the raw value for rates and indexes. Forecast paths
# are arrays draws x H x n, horizons 1..H after the last observed quarter T.
#   exact annualized QoQ growth: 100 ((x_h / x_{h-1})^4 - 1) = 100 (exp(4 (y_h - y_{h-1}) / 100) - 1)
#   YoY growth (log approximation): y_h - y_{h-4}
#
# Inputs:  engine outputs (outputs/<run_id>/engine/), the run's draws and
#          forecast history, the SPF median files (data/raw/spf/).
# Outputs: data.frames and figures for reports/forecast_review.qmd.


# ---- Growth measures ----

#' Exact annualized quarter-on-quarter growth of log series, path by path.
#'
#' The SPF reports annualized percent changes, so the comparison uses the
#' exact formula rather than the log approximation 4 (y_h - y_{h-1}).
#'
#' @param paths Array draws x H x n (model units).
#' @param history Matrix of observed data (quarters x n); its last row is quarter T.
#' @param is_log Logical vector, one per variable.
#' @return Array draws x H x n; NA for raw series.
qoq_annualized_exact <- function(paths, history, is_log) {
  out <- paths
  for (j in seq_len(dim(paths)[3])) {
    if (!is_log[j]) { out[, , j] <- NA_real_; next }
    p <- matrix(paths[, , j], nrow = dim(paths)[1])
    d <- p - cbind(history[nrow(history), j], p[, -ncol(p), drop = FALSE])
    out[, , j] <- 100 * (exp(4 * d / 100) - 1)
  }
  out
}


# ---- SPF ----

#' SPF median forecasts of one survey, as horizons of our forecast.
#'
#' Column 2 of each SPF sheet is the survey quarter, columns 3-6 the next
#' four quarters (SPF documentation): with our last observed quarter the
#' quarter before the survey, column k is our horizon k - 1.
#'
#' @param level_file,growth_file The SPF median files.
#' @param year,quarter The survey.
#' @param map data.frame: variable (ours), file ("level" or "growth"),
#'   sheet, prefix (column prefix, e.g. "drgdp", "UNEMP").
#' @return data.frame: variable, horizon (1-5), spf.
read_spf_survey <- function(level_file, growth_file, year, quarter, map) {
  rows <- lapply(seq_len(nrow(map)), function(k) {
    m <- map[k, ]
    f <- if (m$file == "level") level_file else growth_file
    x <- as.data.frame(readxl::read_excel(f, sheet = m$sheet, na = "#N/A"))
    r <- x[x$YEAR == year & x$QUARTER == quarter, , drop = FALSE]
    if (nrow(r) != 1) stop("survey ", year, "Q", quarter, " not found in sheet ", m$sheet)
    data.frame(variable = m$variable, horizon = 1:5,
               spf = as.numeric(unlist(r[paste0(m$prefix, 2:6)])))
  })
  do.call(rbind, rows)
}

#' Our forecast distribution against the SPF median, per variable and horizon.
#'
#' @param measure Array draws x H x n of our forecasts in the SPF's measure.
#' @param spf Output of read_spf_survey().
#' @return spf with our q05, q16, q50, q84, q95, the difference
#'   spf - q50, and the SPF median's percentile in our distribution.
compare_spf <- function(measure, spf) {
  rows <- lapply(seq_len(nrow(spf)), function(k) {
    v <- spf$variable[k]; h <- spf$horizon[k]
    x <- measure[, h, v]
    q <- quantile(x, c(0.05, 0.16, 0.5, 0.84, 0.95), names = FALSE)
    data.frame(spf[k, ], q05 = q[1], q16 = q[2], q50 = q[3], q84 = q[4], q95 = q[5],
               spf_minus_median = spf$spf[k] - q[3], spf_percentile = 100 * mean(x <= spf$spf[k]))
  })
  do.call(rbind, rows)
}


# ---- History-based ranges ----

#' Ranges of each variable in history: YoY growth for log series, level for rates.
#'
#' @param y Data matrix (quarters x n, model units), rows named by quarter.
#' @param is_log Logical vector, one per variable.
#' @param end Last quarter used (default leaves out the 2020-21 swings).
#' @param probs Percentiles of the range.
#' @return data.frame: variable, measure, mean, low, high.
historical_ranges <- function(y, is_log, end = "2019Q4", probs = c(0.05, 0.95)) {
  y <- y[seq_len(match(end, rownames(y))), , drop = FALSE]
  do.call(rbind, lapply(seq_len(ncol(y)), function(j) {
    x <- if (is_log[j]) diff(y[, j], lag = 4) else y[, j]
    q <- quantile(x, probs, names = FALSE)
    data.frame(variable = colnames(y)[j], measure = if (is_log[j]) "YoY growth, %" else "level",
               mean = mean(x), low = q[1], high = q[2])
  }))
}

#' Share of forecast draws outside the historical range, and the median.
#'
#' @param annual Array draws x H x n in the annual measure (growth_measure()).
#' @param ranges Output of historical_ranges().
#' @param horizons Horizons checked.
#' @return data.frame: variable, horizon, median, low, high, share_below,
#'   share_above, median_outside.
range_check <- function(annual, ranges, horizons) {
  rows <- list()
  for (v in ranges$variable) for (h in horizons) {
    x <- annual[, h, v]; r <- ranges[ranges$variable == v, ]
    m <- median(x)
    rows[[length(rows) + 1]] <- data.frame(variable = v, horizon = h, median = m, low = r$low, high = r$high,
                                           share_below = mean(x < r$low), share_above = mean(x > r$high),
                                           median_outside = m < r$low || m > r$high)
  }
  do.call(rbind, rows)
}

#' Share of draws below zero, by horizon.
#'
#' @param paths Array draws x H x n.
#' @param variables Variables (e.g. the yields).
#' @return data.frame: variable, horizon, share_negative.
negative_share <- function(paths, variables) {
  do.call(rbind, lapply(variables, function(v)
    data.frame(variable = v, horizon = seq_len(dim(paths)[2]),
               share_negative = colMeans(paths[, , v] < 0))))
}

#' Quantiles of a measure at chosen horizons.
#'
#' @param measure Array draws x H x n.
#' @param variables,horizons Selection.
#' @param weights Optional weights per path (tilted forecasts); quantiles weighted
#'   as by the engine (weighted_quantiles(), R/tilting.R).
#' @return data.frame: variable, horizon, q05, q16, q50, q84, q95.
measure_quantiles <- function(measure, variables, horizons, weights = NULL) {
  probs <- c(0.05, 0.16, 0.5, 0.84, 0.95)
  rows <- list()
  for (v in variables) for (h in horizons) {
    x <- measure[, h, v]
    q <- if (is.null(weights)) quantile(x, probs, names = FALSE) else weighted_quantiles(x, weights, probs)
    rows[[length(rows) + 1]] <- data.frame(variable = v, horizon = h, q05 = q[1], q16 = q[2], q50 = q[3],
                                           q84 = q[4], q95 = q[5])
  }
  do.call(rbind, rows)
}

# ---- Simulation noise ----

#' Difference of medians between two forecast runs that differ only in the seed.
#'
#' @param a,b Arrays draws x H x n in the same measure.
#' @param variables Variables.
#' @return data.frame: variable, horizon, median_a, diff (b - a).
seed_difference <- function(a, b, variables) {
  do.call(rbind, lapply(variables, function(v) {
    ma <- apply(a[, , v, drop = FALSE], 2, median); mb <- apply(b[, , v, drop = FALSE], 2, median)
    data.frame(variable = v, horizon = seq_along(ma), median_a = ma, diff = mb - ma)
  }))
}


# ---- Plots (base graphics; colours of R/estimate_checks.R) ----

#' Our forecast bands against the SPF median, one panel per variable.
#'
#' @param cmp Output of compare_spf() for one COVID method.
#' @param variables Variables shown.
#' @param units Panel label per variable.
#' @param title Figure title.
plot_spf_panels <- function(cmp, variables, units, title) {
  op <- par(mfrow = c(ceiling(length(variables) / 4), 4), mar = c(3, 3, 2.2, 0.8), oma = c(2.5, 0, 2.5, 0),
            col.axis = ink[["secondary"]], mgp = c(1.8, 0.5, 0))
  on.exit(par(op))
  for (v in variables) {
    d <- cmp[cmp$variable == v, ]
    ylim <- range(c(d$q05, d$q95, d$spf), na.rm = TRUE)
    plot(NA, xlim = c(1, 5), ylim = ylim, bty = "n", xlab = "quarters ahead", ylab = "",
         main = paste0(v, " (", units[[v]], ")"), font.main = 1, cex.main = 0.95, col.main = ink[["primary"]],
         cex.axis = 0.8)
    polygon(c(d$horizon, rev(d$horizon)), c(d$q05, rev(d$q95)), col = "#2a78d633", border = NA)
    polygon(c(d$horizon, rev(d$horizon)), c(d$q16, rev(d$q84)), col = "#2a78d666", border = NA)
    lines(d$horizon, d$q50, col = check_colours[["blue"]], lwd = 2)
    points(d$horizon, d$spf, pch = 19, col = check_colours[["orange"]])
    lines(d$horizon, d$spf, col = check_colours[["orange"]], lty = 2)
  }
  mtext(title, outer = TRUE, line = 0.8, cex = 1.05, col = ink[["primary"]])
  mtext("Blue: model median, 68% and 90% bands. Orange dots: SPF median (2026Q3 survey).",
        side = 1, outer = TRUE, line = 1, cex = 0.8, col = ink[["secondary"]])
}

#' Long-horizon medians and bands of one measure, every method, against history.
#'
#' @param q Named list (method) of measure_quantiles() outputs over horizons 1..H.
#' @param variables Variables shown.
#' @param ranges historical_ranges() (mean and 5-95 range drawn as lines).
#' @param title Figure title.
plot_long_horizon <- function(q, variables, ranges, title) {
  op <- par(mfrow = c(ceiling(length(variables) / 3), 3), mar = c(3, 3, 2.2, 0.8), oma = c(2.5, 0, 2.5, 0),
            col.axis = ink[["secondary"]], mgp = c(1.8, 0.5, 0))
  on.exit(par(op))
  cols <- c(check_colours[["blue"]], check_colours[["orange"]], check_colours[["aqua"]])
  for (v in variables) {
    r <- ranges[ranges$variable == v, ]
    d <- lapply(q, function(x) x[x$variable == v, ])
    ylim <- range(c(unlist(lapply(d, function(x) c(x$q05, x$q95))), r$low, r$high))
    plot(NA, xlim = range(d[[1]]$horizon), ylim = ylim, bty = "n", xlab = "quarters ahead", ylab = "",
         main = paste0(v, " (", r$measure, ")"), font.main = 1, cex.main = 0.95, col.main = ink[["primary"]],
         cex.axis = 0.8)
    abline(h = c(r$low, r$high), col = ink[["grid"]], lwd = 2)
    abline(h = r$mean, col = ink[["secondary"]], lty = 3)
    for (k in seq_along(d)) {
      lines(d[[k]]$horizon, d[[k]]$q50, col = cols[k], lwd = 2)
      lines(d[[k]]$horizon, d[[k]]$q05, col = cols[k], lty = 3)
      lines(d[[k]]$horizon, d[[k]]$q95, col = cols[k], lty = 3)
    }
  }
  mtext(title, outer = TRUE, line = 0.8, cex = 1.05, col = ink[["primary"]])
  mtext(paste0("Median (solid) and 90% band (dotted): ",
               paste(c("blue", "orange", "green")[seq_along(q)], names(q), collapse = ", "),
               ". Grey: 1973-2019 5th and 95th percentiles; dotted grey: 1973-2019 mean."),
        side = 1, outer = TRUE, line = 1, cex = 0.8, col = ink[["secondary"]])
}


# ---- Readable forecast table ----

#' Baseline forecast table: median and 68 percent band per variable and quarter.
#'
#' Log series as year-on-year percent changes, computed exactly per path,
#' 100 (exp((y_h - y_{h-4}) / 100) - 1), so the numbers compare with
#' published forecasts; rates and indexes as levels.
#'
#' @param paths Array draws x H x n (model units).
#' @param history Observed data (quarters x n); its last 4 rows precede h = 1.
#' @param is_log Logical vector, one per variable.
#' @param variables Variables in the table.
#' @param quarters Labels of the forecast quarters (length H).
#' @return data.frame: variable, measure, quarter, q16, q50, q84.
forecast_table <- function(paths, history, is_log, variables, quarters) {
  yoy <- growth_measure(paths, history, is_log, "annual")       # y_h - y_{h-4}, log series
  rows <- list()
  for (v in variables) {
    j <- match(v, dimnames(paths)[[3]])
    x <- if (is_log[j]) 100 * (exp(yoy[, , j] / 100) - 1) else paths[, , j]
    q <- apply(matrix(x, nrow = dim(paths)[1]), 2, quantile, c(0.16, 0.5, 0.84), names = FALSE)
    rows[[v]] <- data.frame(variable = v, measure = if (is_log[j]) "YoY % change" else "level",
                            quarter = quarters, q16 = q[1, ], q50 = q[2, ], q84 = q[3, ])
  }
  do.call(rbind, rows)
}
