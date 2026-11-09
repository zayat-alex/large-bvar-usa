# chart_data.R: the numbers behind one chart, taken from the bundle.
#
# A chart shows one variable in one transform: its recent history, then for
# each scenario the reader has switched on the forecast median and,
# optionally, the 68% and 90% credible bands (quantiles 16-84 and 5-95 of
# the simulated paths, computed path by path in the export). Every forecast
# starts at the last observed quarter, at the observed value, so the lines
# join the history without a gap.
#
# In the difference view a chart shows, for each scenario, its median minus
# the baseline's median, quarter by quarter: a difference between two
# conditional forecasts, not a causal effect. It is a difference of
# medians, not the median of a difference, so it has no band; the spread of
# each scenario is on its own fan chart.


#' Transform in which a variable is shown.
#'
#' The export marks each variable: "level" (rates, raw indices, spreads),
#' "natural" (log series read as levels, such as the oil price) or "growth"
#' (other log series, shown in the reader's growth setting).
#'
#' @param vars The bundle's variables table.
#' @param variable Variable name.
#' @param growth Reader's growth setting: "yoy", "qoq_ann" or "natural".
#' @return "level", "natural", "yoy" or "qoq_ann".
display_transform <- function(vars, variable, growth) {
  d <- vars$display[match(variable, vars$variable)]
  if (is.na(d)) stop("unknown variable: ", variable)
  if (d == "growth") growth else d
}

#' History and forecasts of one variable, ready to plot.
#'
#' @param bundle The bundle.
#' @param variable Variable name.
#' @param transform Transform (display_transform()).
#' @param objects Ids of the scenarios to include.
#' @param history_quarters Quarters of history to keep (Inf = all).
#' @return List: history (data.frame quarter, date, value), last (its last
#'   row), forecasts (named list by object of data.frames quarter, date,
#'   q05, q16, q50, q84, q95, the first row being the last observation).
chart_series <- function(bundle, variable, transform, objects, history_quarters = Inf) {
  h <- bundle$history[bundle$history$variable == variable & bundle$history$transform == transform,
                      c("quarter", "date", "value")]
  if (!nrow(h)) stop("no history for ", variable, " (", transform, ")")
  h <- h[order(h$date), ]
  if (is.finite(history_quarters)) h <- utils::tail(h, history_quarters)
  last <- h[nrow(h), ]
  qcols <- c("q05", "q16", "q50", "q84", "q95")
  f <- bundle$forecasts[bundle$forecasts$variable == variable & bundle$forecasts$transform == transform, ]
  forecasts <- lapply(stats::setNames(objects, objects), function(obj) {
    x <- f[f$object == obj, c("quarter", "date", qcols)]
    if (!nrow(x)) stop("no forecast of ", variable, " (", transform, ") for ", obj)
    x <- x[order(x$date), ]
    # Anchor: the last observation, with every quantile equal to it.
    anchor <- data.frame(quarter = last$quarter, date = last$date, q05 = last$value, q16 = last$value,
                         q50 = last$value, q84 = last$value, q95 = last$value)
    rbind(anchor, x)
  })
  rownames(h) <- NULL
  list(history = h, last = last, forecasts = forecasts)
}

#' First day of a quarter ("2026Q3" -> 2026-07-01).
quarter_date <- function(q) as.Date(sprintf("%s-%02d-01", substr(q, 1, 4), 3L * as.integer(substr(q, 6, 6)) - 2L))

#' What a scenario assumes, as points to mark on a chart.
#'
#' Hard scenarios fix some cells exactly (the conditions table): they are
#' assumptions, not model results, so the chart marks them. A tilted
#' scenario's soft targets (the targets table) are marked where the chart
#' shows the target's measure: a median target in yoy on a yoy chart, for
#' example. A fixed level cannot be marked on a growth chart (its growth
#' also depends on the quarter before), so fixed cells appear on level
#' charts only.
#'
#' @param bundle The bundle.
#' @param variable,transform The chart's variable and transform.
#' @param ids Ids of the scenarios drawn.
#' @return data.frame object, quarter, date, value, kind ("fixed" or
#'   "target"); zero rows if nothing to mark.
chart_markers <- function(bundle, variable, transform, ids) {
  out <- data.frame(object = character(0), quarter = character(0), date = as.Date(character(0)),
                    value = numeric(0), kind = character(0))
  cn <- bundle$conditions
  if (!is.null(cn) && nrow(cn)) {
    x <- cn[cn$variable == variable & cn$transform == transform & cn$object %in% ids, ]
    if (nrow(x)) out <- rbind(out, data.frame(object = x$object, quarter = x$quarter, date = quarter_date(x$quarter),
                                              value = x$value, kind = "fixed"))
  }
  tg <- bundle$targets
  if (!is.null(tg) && nrow(tg) && all(c("measure", "statistic", "target") %in% names(tg))) {
    x <- tg[tg$variable == variable & tg$measure == transform & tg$statistic == "median" & tg$object %in% ids, ]
    if (nrow(x)) out <- rbind(out, data.frame(object = x$object, quarter = x$quarter, date = quarter_date(x$quarter),
                                              value = x$target, kind = "target"))
  }
  out
}

#' The numbers of a tab's charts, as one table (for the Data table and the
#' CSV download).
#'
#' Forecast view: for each variable, the observed history in the chart's
#' window, then each shown scenario's median and its 68% and 90% bands, in
#' the transform the chart shows. Difference view: each shown scenario's
#' median minus the baseline's. The observed history of the variables in
#' APP_NO_HISTORY_DOWNLOAD is left out (the publishers' terms).
#'
#' @param bundle The bundle.
#' @param variables The tab's variables.
#' @param s Settings (show, baseline, mode, growth, history).
#' @param scenarios The listed scenarios (id, label, colour).
#' @return data.frame: Variable, Unit, Series, Quarter, Value, and in the
#'   forecast view the band columns.
tab_data <- function(bundle, variables, s, scenarios) {
  difference <- identical(s$mode, "difference") && !is.null(s$baseline) && s$baseline %in% scenarios$id
  shown <- scenarios[scenarios$id %in% s$show, , drop = FALSE]
  rows <- list()
  for (v in variables) {
    tr <- display_transform(bundle$variables, v, s$growth)
    ids <- if (difference) union(shown$id, s$baseline) else shown$id
    series <- chart_series(bundle, v, tr, ids, APP_HISTORY[[s$history]])
    unit <- unit_label(v, tr, difference)
    if (difference) {
      base_label <- scenarios$label[scenarios$id == s$baseline]
      d <- median_differences(series, s$baseline)
      for (k in seq_len(nrow(shown))) {
        x <- d[[shown$id[k]]]
        if (is.null(x)) next
        x <- x[-1, ]                                                  # drop the anchor (zero)
        rows[[length(rows) + 1]] <- data.frame(Variable = var_label(v), Unit = unit,
                                               Series = paste(shown$label[k], "minus", base_label),
                                               Quarter = x$quarter, Value = x$value)
      }
    } else {
      if (!v %in% APP_NO_HISTORY_DOWNLOAD) {
        h <- series$history
        rows[[length(rows) + 1]] <- data.frame(Variable = var_label(v), Unit = unit, Series = "Observed",
                                               Quarter = h$quarter, Value = h$value, `68% low` = NA, `68% high` = NA,
                                               `90% low` = NA, `90% high` = NA, check.names = FALSE)
      }
      for (k in seq_len(nrow(shown))) {
        f <- series$forecasts[[shown$id[k]]][-1, ]                    # drop the anchor (the last observation)
        rows[[length(rows) + 1]] <- data.frame(Variable = var_label(v), Unit = unit,
                                               Series = paste(shown$label[k], "(median)"),
                                               Quarter = f$quarter, Value = f$q50, `68% low` = f$q16,
                                               `68% high` = f$q84, `90% low` = f$q05, `90% high` = f$q95,
                                               check.names = FALSE)
      }
    }
  }
  if (!length(rows)) return(data.frame(Variable = character(0), Unit = character(0), Series = character(0),
                                       Quarter = character(0), Value = numeric(0)))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Differences of medians from the baseline.
#'
#' @param series Output of chart_series(), including the baseline.
#' @param baseline Id of the baseline.
#' @return Named list by object (baseline excluded) of data.frames quarter,
#'   date, value = median minus the baseline's median (0 at the anchor).
median_differences <- function(series, baseline) {
  base <- series$forecasts[[baseline]]
  if (is.null(base)) stop("the baseline ", baseline, " is not among the chart's forecasts")
  others <- setdiff(names(series$forecasts), baseline)
  lapply(stats::setNames(others, others), function(obj) {
    x <- series$forecasts[[obj]]
    if (!identical(x$quarter, base$quarter)) stop("forecast quarters differ between ", obj, " and ", baseline)
    data.frame(quarter = x$quarter, date = x$date, value = x$q50 - base$q50)
  })
}
