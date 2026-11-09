# site_data.R: the static website's data, taken from the export bundle.
#
# The static site (site/) shows the same numbers as the Shiny app (app/),
# drawn in the browser instead of on a server. It reads two JSON files
# written here from the app's bundle:
#   core.json  settings, variables, scenarios, observed history, forecast
#              quantiles, fixed cells and targets (loaded with the page);
#   irf.json   impulse-response quantiles (loaded when the IRF tab opens).
# Labels, colours, tabs and defaults come from the app's own settings
# (app/R/settings.R, app/R/irf_data.R), so the two sites cannot drift
# apart. Only what the app can display is written: for each variable the
# transforms its charts can show, quantiles 5-16-50-84-95, and responses up
# to the longest horizon the IRF tab offers.
#
# Layout of the numbers (columnar, to keep the files small):
#   history[[variable]][[transform]]           list(start = "1973Q1", values)
#   forecasts[[object]][[variable]][[transform]] list(q05, q16, q50, q84, q95),
#       one value per forecast quarter (meta$forecast_quarters);
#   responses[[shock]][[variable]][[transform]]  list(q05, ..., q95),
#       one value per horizon 0, 1, ..., horizon_max.


# Numbers are written with 17 significant digits, which is enough to read
# back exactly the double R holds, so the browser rounds for display the
# same number the app rounds. Fewer digits can move a shown digit: at 15,
# 3210.7574999999988 (shown as 3210.757) is written as 3210.7575 and shown
# as 3210.758; rounding to six decimals moves -0.5135001 to -0.5135, shown
# as -0.513 instead of -0.514.

#' Quantile columns the site uses (68% and 90% bands and the median).
SITE_QUANTILES <- c("q05", "q16", "q50", "q84", "q95")

#' Transforms in which the app can show a variable.
#'
#' @param display The variable's display mark in the bundle: "level",
#'   "natural" or "growth" (shown in the reader's growth setting).
#' @return Character vector of transforms.
site_transforms <- function(display) {
  switch(display, level = "level", natural = "natural", growth = c("yoy", "qoq_ann", "natural"),
         stop("unknown display mark: ", display))
}

#' Quarter labels from one to another ("2026Q3", ...), consecutive.
#'
#' @param from,to Quarters as "YYYYQn".
#' @return Character vector.
site_quarter_seq <- function(from, to) {
  idx <- function(q) as.integer(substr(q, 1, 4)) * 4L + as.integer(substr(q, 6, 6)) - 1L
  k <- seq(idx(from), idx(to))
  sprintf("%dQ%d", k %/% 4L, k %% 4L + 1L)
}

#' Run settings and data facts the site shows (footer, About page).
#'
#' @param bundle The bundle.
#' @return List.
site_meta <- function(bundle) {
  m <- bundle$meta
  list(run_id = m$run_id,
       vintage = sub("^FRED-QD_(.*)\\.csv$", "\\1", basename(m$data$qd_vintage_file)),
       api_realtime_date = m$data$api_realtime_date,
       first_quarter = m$data$first_quarter,
       last_observed = m$data$last_observed,
       lags = m$model$lags,
       covid_window = as.list(unlist(m$model$covid_window)),
       horizon_exported = m$forecast$horizon_exported,
       n_draws = m$forecast$n_draws,
       paths_per_draw = m$forecast$paths_per_draw,
       forecast_quarters = as.list(site_quarter_seq(m$forecast$quarters[1], m$forecast$quarters[2])))
}

#' The app's presentation settings (app/R/settings.R, app/R/irf_data.R).
#'
#' @param app Environment holding the app's functions and settings.
#' @return List.
site_settings <- function(app) {
  hist <- app$APP_HISTORY
  list(palette = as.list(app$APP_PALETTE), palette_extra = app$APP_PALETTE_EXTRA,
       overview_default = as.list(app$APP_OVERVIEW_DEFAULT), overview_max = app$APP_OVERVIEW_MAX,
       # Inf (all history) is written as null.
       history = lapply(names(hist), function(n) list(name = n, quarters = if (is.finite(hist[[n]])) hist[[n]] else NULL)),
       no_history_download = as.list(app$APP_NO_HISTORY_DOWNLOAD),
       topics = lapply(app$APP_TOPICS, function(t) list(id = t$id, title = t$title, variables = as.list(t$variables))),
       irf_default_shock = app$APP_IRF_DEFAULT_SHOCK,
       irf_default_responses = lapply(app$APP_IRF_DEFAULT_RESPONSES, as.list),
       irf_general_responses = as.list(app$APP_IRF_GENERAL_RESPONSES),
       irf_horizons = as.list(app$APP_IRF_HORIZONS))
}

#' Variables with their labels, units, sources and display marks.
#'
#' @param bundle The bundle.
#' @param app The app's environment.
#' @return data.frame, one row per variable in the app's order.
site_variables <- function(bundle, app) {
  v <- bundle$variables
  av <- app$APP_VARIABLES
  src <- app$APP_SOURCES
  vars <- intersect(av$variable, v$variable)
  i <- match(vars, v$variable)
  data.frame(variable = vars, label = av$label[match(vars, av$variable)], description = v$label[i],
             unit = av$unit[match(vars, av$variable)], kind = v$kind[i], model_transform = v$model_transform[i],
             display = v$display[i], source = src$source[match(vars, src$variable)],
             note = src$note[match(vars, src$variable)], stringsAsFactors = FALSE)
}

#' Observed history, per variable and transform.
#'
#' @param bundle The bundle.
#' @param variables Output of site_variables().
#' @return Nested list (see the file header). Stops if a history has a gap,
#'   since the site rebuilds the quarters from the first one.
site_history <- function(bundle, variables) {
  h <- bundle$history
  out <- list()
  for (k in seq_len(nrow(variables))) {
    v <- variables$variable[k]
    for (tr in site_transforms(variables$display[k])) {
      x <- h[h$variable == v & h$transform == tr, ]
      if (!nrow(x)) stop("no history of ", v, " (", tr, ")")
      x <- x[order(x$date), ]
      if (!identical(x$quarter, site_quarter_seq(x$quarter[1], x$quarter[nrow(x)]))) {
        stop("the history of ", v, " (", tr, ") skips quarters")
      }
      out[[v]][[tr]] <- list(start = x$quarter[1], values = x$value)
    }
  }
  out
}

#' Forecast quantiles, per object, variable and transform.
#'
#' @param bundle The bundle.
#' @param variables Output of site_variables().
#' @param quarters The forecast quarters (site_meta()$forecast_quarters).
#' @return Nested list (see the file header). Stops if a forecast does not
#'   cover exactly the forecast quarters.
site_forecasts <- function(bundle, variables, quarters) {
  f <- bundle$forecasts
  out <- list()
  for (obj in bundle$objects$id) {
    for (k in seq_len(nrow(variables))) {
      v <- variables$variable[k]
      for (tr in site_transforms(variables$display[k])) {
        x <- f[f$object == obj & f$variable == v & f$transform == tr, ]
        x <- x[order(x$date), ]
        if (!identical(x$quarter, quarters)) stop("the forecast of ", v, " (", tr, ") for ", obj, " has other quarters")
        out[[obj]][[v]][[tr]] <- lapply(stats::setNames(SITE_QUANTILES, SITE_QUANTILES), function(q) x[[q]])
      }
    }
  }
  out
}

#' Cells fixed by hard scenarios and median targets of tilted ones, in the
#' transforms the site can show (the markers of the app's chart_markers()).
#'
#' @param bundle The bundle.
#' @param variables Output of site_variables().
#' @return List: conditions (data.frame object, variable, transform,
#'   quarter, value) and targets (object, variable, transform, quarter,
#'   value); a target's transform is its measure.
site_markers <- function(bundle, variables) {
  shown <- do.call(rbind, lapply(seq_len(nrow(variables)), function(k) {
    data.frame(variable = variables$variable[k], transform = site_transforms(variables$display[k]))
  }))
  key <- paste(shown$variable, shown$transform)
  cn <- bundle$conditions
  cn <- cn[paste(cn$variable, cn$transform) %in% key, ]
  tg <- bundle$targets
  tg <- tg[tg$statistic == "median" & paste(tg$variable, tg$measure) %in% key, ]
  list(conditions = data.frame(object = cn$object, variable = cn$variable, transform = cn$transform,
                               quarter = cn$quarter, value = cn$value),
       targets = data.frame(object = tg$object, variable = tg$variable, transform = tg$measure,
                            quarter = tg$quarter, value = tg$target))
}

#' Impulse-response quantiles up to the longest horizon the site offers.
#'
#' @param bundle The bundle.
#' @param app The app's environment.
#' @return List: horizon_max, shocks (the app's irf_shocks() table) and
#'   responses (see the file header).
site_irf <- function(bundle, app) {
  hmax <- max(app$APP_IRF_HORIZONS)
  r <- bundle$irf[bundle$irf$horizon <= hmax, ]
  responses <- list()
  for (id in unique(r$identification)) {
    x <- r[r$identification == id, ]
    for (v in unique(x$variable)) {
      for (tr in unique(x$transform[x$variable == v])) {
        y <- x[x$variable == v & x$transform == tr, ]
        y <- y[order(y$horizon), ]
        if (!identical(as.integer(y$horizon), 0:as.integer(hmax))) stop("responses of ", v, " to ", id, " skip horizons")
        responses[[id]][[v]][[tr]] <- lapply(stats::setNames(SITE_QUANTILES, SITE_QUANTILES),
                                             function(q) y[[q]])
      }
    }
  }
  list(horizon_max = hmax, shocks = app$irf_shocks(bundle), responses = responses)
}

#' Everything in core.json.
#'
#' @param bundle The bundle.
#' @param app The app's environment.
#' @return List.
site_core <- function(bundle, app) {
  meta <- site_meta(bundle)
  variables <- site_variables(bundle, app)
  markers <- site_markers(bundle, variables)
  list(meta = meta, settings = site_settings(app), variables = variables,
       scenarios = app$scenario_catalog(bundle),
       history = site_history(bundle, variables),
       forecasts = site_forecasts(bundle, variables, unlist(meta$forecast_quarters)),
       conditions = markers$conditions, targets = markers$targets)
}

#' Write core.json and irf.json.
#'
#' @param bundle The bundle.
#' @param app The app's environment.
#' @param dir Output folder.
#' @return Paths of the two files, invisibly.
write_site_data <- function(bundle, app, dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  paths <- file.path(dir, c("core.json", "irf.json"))
  # digits = I(17): 17 significant digits (see the note at the top of the
  # file); data frames are written row by row, NA as null.
  write <- function(x, path) {
    jsonlite::write_json(x, path, auto_unbox = TRUE, digits = I(17), na = "null", null = "null", dataframe = "rows")
  }
  write(site_core(bundle, app), paths[1])
  write(site_irf(bundle, app), paths[2])
  invisible(paths)
}
