# tilting_checks.R: checks of entropic tilting (R/tilting.R) on the
# model's own posterior draws, for reports/tilting_checks.qmd.
#
# Why it exists: the tests (tests/testthat/test-tilting.R) check the
# method on small cases with known answers. This file applies it to the
# full model: the inflation anchor of Crump et al. (2021), Section 4.2
# (Figure 16b: the median of 4-quarter headline PCE inflation 12 quarters
# ahead tilted to 2%), how many paths the tilted result really rests on
# (effective sample size, weight concentration), how the result depends on
# the number of paths per posterior draw, what the anchor implies for other
# variables, and a scenario that mixes hard and soft conditions.
#
# Inputs:  engine inputs (engine_inputs()), settings, scenario files.
# Outputs: summary tables and quantiles for the report, and figures.


#' 4-quarter percent change of a log series on every path and quarter.
#'
#' @param paths Array J x H x n (model units).
#' @param variable Short name of a log series.
#' @param history Data matrix (base quarters up to T).
#' @return Matrix J x H of exact yoy percent changes.
yoy_paths <- function(paths, variable, history) {
  sapply(seq_len(dim(paths)[2]), function(h) tilt_quantity(paths, variable, h, "yoy", history))
}

#' Observed 4-quarter percent change over the last k quarters of the history.
#'
#' @param history Data matrix (model units), row names = quarters.
#' @param variable Short name of a log series.
#' @param k Number of quarters.
#' @return Named numeric vector (quarters).
yoy_history <- function(history, variable, k) {
  y <- history[, variable]
  n <- length(y)
  idx <- (n - k + 1):n
  stats::setNames(100 * (exp((y[idx] - y[idx - 4]) / 100) - 1), rownames(history)[idx])
}

#' Quantiles per column, weighted or not.
#'
#' @param x Matrix J x H.
#' @param weights NULL or J weights.
#' @param probs Probabilities.
#' @return Matrix H x length(probs).
column_quantiles <- function(x, weights, probs) {
  t(apply(x, 2, weighted_quantiles, weights, probs))
}

#' How concentrated tilting weights are.
#'
#' Paths are stored draw by draw (all paths of draw 1 first), so the weight
#' of a posterior draw is the sum of its paths' weights.
#'
#' @param w Tilting weights (length J).
#' @param n_paths Paths per posterior draw.
#' @return One-row data.frame: J, ESS and its share, largest weight relative
#'   to 1/J, the share of paths and of posterior draws that carry half of
#'   the weight (0.5 when the weights are equal).
weight_concentration <- function(w, n_paths) {
  J <- length(w)
  half <- function(v) {
    v <- sort(v, decreasing = TRUE)
    which(cumsum(v) >= 0.5)[1] / length(v)
  }
  draw_w <- tapply(w, rep(seq_len(J / n_paths), each = n_paths), sum)
  data.frame(paths = J, ess = tilt_ess(w), ess_share = tilt_ess(w) / J, max_weight_x_J = max(w) * J,
             paths_for_half = half(w), draws_for_half = half(draw_w))
}

#' Run one tilted scenario and collect what the report needs.
#'
#' @param conditions Hard conditions (H x n, NA = free).
#' @param targets Soft targets (soft_targets()).
#' @param inputs Output of engine_inputs().
#' @param settings Output of engine_settings().
#' @param baseline_mean The baseline's per-draw means (draws x H x n).
#' @param fan_vars Variables whose 4-quarter growth (log series) or level
#'   (rates) fans are kept.
#' @param series Series configuration (transforms).
#' @param level_vars Log series to show as levels (model units) instead of
#'   growth rates (for example oil, whose scenario cells are prices).
#' @return List: check (targets met), concentration, seconds, fixed_max_error,
#'   fans (per variable: untilted and tilted quantiles H x 5, the measure),
#'   at (untilted and tilted medians of the fan variables in each target
#'   quarter).
tilting_case <- function(conditions, targets, inputs, settings, baseline_mean, fan_vars, series,
                         level_vars = character(0)) {
  s <- engine_scenario(conditions, inputs, settings, baseline_mean, targets)
  w <- s$tilt$weights
  probs <- c(0.05, 0.16, 0.25, 0.5, 0.75, 0.84, 0.95)
  fans <- lapply(stats::setNames(fan_vars, fan_vars), function(v) {
    as_growth <- series$transform[match(v, series$short_name)] == "log" && !v %in% level_vars
    x <- if (as_growth) yoy_paths(s$paths, v, inputs$history) else s$paths[, , v]
    list(measure = if (as_growth) "yoy" else "level",
         untilted = column_quantiles(x, NULL, probs), tilted = column_quantiles(x, w, probs))
  })
  h_t <- unique(targets$h)
  at <- do.call(rbind, lapply(fan_vars, function(v) {
    f <- fans[[v]]
    data.frame(variable = v, measure = f$measure, quarter = inputs$quarters[h_t], h = h_t,
               untilted_median = f$untilted[h_t, 4], tilted_median = f$tilted[h_t, 4],
               untilted_q16_q84 = f$untilted[h_t, 6] - f$untilted[h_t, 2],
               tilted_q16_q84 = f$tilted[h_t, 6] - f$tilted[h_t, 2])
  }))
  list(check = s$tilt$check, iterations = s$tilt$iterations, kl = s$tilt$kl,
       concentration = weight_concentration(w, settings$tilting$paths_per_draw),
       seconds = s$seconds, fixed_max_error = s$fixed_max_error, fans = fans, at = at,
       probs = probs)
}

#' Fan chart of untilted and tilted quantiles, joined to the history.
#'
#' Grey: untilted 90% band and median; blue: tilted 90% and 68% bands and
#' median. Both start at the last observation. Dots mark soft targets
#' (medians or means) and, optionally, hard cells.
#'
#' @param fan Element of tilting_case()$fans.
#' @param hist_values Observed values of the same measure (last quarters).
#' @param title Panel title.
#' @param target_h,target Quarters (h) and values of targets to mark.
#' @param hard_h,hard Quarters and values of hard cells to mark.
plot_tilt_fan <- function(fan, hist_values, title, target_h = NULL, target = NULL,
                          hard_h = NULL, hard = NULL) {
  qu <- fan$untilted; qt <- fan$tilted
  H <- nrow(qt); k <- length(hist_values)
  last <- hist_values[k]
  x <- c(0, seq_len(H))
  qu <- rbind(rep(last, ncol(qu)), qu)
  qt <- rbind(rep(last, ncol(qt)), qt)
  ylim <- range(c(qu, qt, hist_values, target, hard), na.rm = TRUE)
  plot(seq(-k + 1, 0), hist_values, type = "l", xlim = c(-k + 1, H), ylim = ylim, bty = "n",
       col = ink[["primary"]], lwd = 1.5, xlab = "quarters from the last observation", ylab = "",
       main = title, font.main = 1, cex.main = 0.9, col.main = ink[["primary"]],
       col.axis = ink[["secondary"]], col.lab = ink[["secondary"]])
  abline(v = 0, col = ink[["grid"]])
  polygon(c(x, rev(x)), c(qu[, 1], rev(qu[, 7])), col = "#52514e22", border = NA)
  polygon(c(x, rev(x)), c(qt[, 1], rev(qt[, 7])), col = "#2a78d633", border = NA)
  polygon(c(x, rev(x)), c(qt[, 2], rev(qt[, 6])), col = "#2a78d666", border = NA)
  lines(x, qu[, 4], col = ink[["secondary"]], lwd = 1.5, lty = 2)
  lines(x, qt[, 4], col = check_colours[["blue"]], lwd = 2)
  if (length(target_h)) points(target_h, target, pch = 21, bg = check_colours[["orange"]], col = ink[["primary"]], cex = 1.2)
  if (length(hard_h)) points(hard_h, hard, pch = 19, col = ink[["primary"]], cex = 0.9)
}

#' Legend for plot_tilt_fan(), across the bottom of the figure.
legend_tilt_fan <- function(hard = FALSE) {
  op <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  on.exit(par(op))
  plot.new()
  lab <- c("untilted 90%", "untilted median", "tilted 90% / 68%", "tilted median", "soft target",
           if (hard) "hard cells")
  legend("bottom", horiz = TRUE, bty = "n", cex = 0.85, text.col = ink[["secondary"]], legend = lab,
         fill = c("#52514e22", NA, "#2a78d655", NA, NA, if (hard) NA), border = NA,
         col = c(NA, ink[["secondary"]], NA, check_colours[["blue"]], ink[["primary"]], if (hard) ink[["primary"]]),
         pt.bg = c(NA, NA, NA, NA, check_colours[["orange"]], if (hard) NA),
         lty = c(NA, 2, NA, 1, NA, if (hard) NA), lwd = c(NA, 1.5, NA, 2, NA, if (hard) NA),
         pch = c(NA, NA, NA, NA, 21, if (hard) 19))
}
