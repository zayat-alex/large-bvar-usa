# export.R: stage 5 of the pipeline, the export. Turns the engine's
# outputs for the published run (<run_dir>/engine/) into the small tables
# the website reads, saved as one file, <run_dir>/bundle/bundle.rds, with a
# readable copy of its metadata in <run_dir>/bundle/meta.json.
#
# Why it exists: the engine keeps every simulated path (10,000 paths x 20
# quarters x 31 variables per object, about 50 MB each). The website needs
# only their quantiles, in the units a reader looks at, for the quarters it
# shows. The app runs no model code; everything it plots is computed here.
#
# The rule: transform first, summarize last. A display transform (a growth
# rate, a spread) is applied to every path (or IRF draw), and the
# quantiles are taken afterwards. Quantiles of growth rates computed from
# quantiles of levels would be wrong: the 5th percentile path of a level is
# not the path of the 5th percentile growth rate.
#
# Display transforms (docs/methodology.md, section 15), with x = exp(y / 100)
# for a log series stored as y = 100 ln x (docs/data.md, section 2):
#   level    model units (100 ln x for log series, natural units for rates);
#   natural  x itself (billions of chained dollars, index points, $/barrel);
#   qoq_ann  100 [(x_t / x_{t-1})^4 - 1], percent change at an annual rate;
#   yoy      100 (x_t / x_{t-4} - 1), percent change from a year earlier.
# These are the exact published conventions (as path_from_growth() in
# R/scenarios.R and the tilting measures in R/tilting.R), not the log
# approximations 4 Δy and Δ4 y. Growth rates exist for log series only. In
# the first forecast quarters the base quarter t-1 or t-4 is observed and
# comes from the history. Spreads (differences of two yields, percentage
# points) are formed path by path and have the level transform only.
#
# Objects exported: the baseline (unconditional forecast) and the scenarios
# listed in config/run.yml (`scenarios`, all of them). A scenario is hard (fixed
# cells, a conditional forecast) or tilted (soft targets imposed by entropic
# tilting, Crump et al. 2021, Section 4.2), as in R/engine_run.R:
#   - forecast quantiles: over paths; weighted by the tilting weights for a
#     tilted object (weighted_quantiles(), R/tilting.R).
# Differences between objects are not exported: the app takes them as the
# difference between two objects' medians, for whichever pair the reader
# chooses, so any object can serve as the reference. (The engine's
# per-draw scenario effects, Crump et al. 2021, Section 2.3, stay in the
# engine files.)
# Impulse responses: quantiles over draws of the response in model units
# per unit shock, and, for log series, of the change in the response over
# four quarters (the response of the 4-quarter log change, in percentage
# points: an approximation, since a response is a deviation, not a level).
#
# Inputs:  <run_dir>/engine/ (baseline.rds, scenario_<id>.rds, irf.rds,
#          irf_extra.rds, meta.json), <run_dir>/covid.rds (the observed history),
#          config/run.yml (export, quantiles), config/series.csv,
#          scenarios/scenarios.yml (labels).
# Outputs: <run_dir>/bundle/bundle.rds and meta.json.


# ---- Settings ----

#' Export settings from the run configuration, checked.
#'
#' @param cfg Parsed config/run.yml.
#' @param series Series configuration (read_series_config()).
#' @return List: run (COVID method of the published run), horizon (forecast
#'   quarters exported), objects (every scenario id in cfg$scenarios), spreads (data.frame id, a,
#'   b, label: the spread is a minus b). Stops on an invalid setting.
export_settings <- function(cfg, series) {
  e <- cfg$export
  if (is.null(e)) stop("export must be set in config/run.yml")
  if (!is.character(e$run) || length(e$run) != 1 || !e$run %in% unlist(cfg$covid$methods)) {
    stop("export$run must be one of covid$methods")
  }
  h <- e$horizon
  if (!is.numeric(h) || length(h) != 1 || h < 1 || h != round(h) || h > cfg$forecast$horizon) {
    stop("export$horizon must be a whole number between 1 and forecast$horizon")
  }
  # Every scenario the engine builds is exported; the website lets the
  # reader choose which ones to list.
  if (!is.null(e$objects)) stop("export$objects is no longer used: every scenario in `scenarios` is exported")
  objects <- as.character(unlist(cfg$scenarios))
  if (anyDuplicated(objects) || "baseline" %in% objects) stop("scenario ids must be distinct and not \"baseline\"")
  spreads <- data.frame(id = character(0), a = character(0), b = character(0), label = character(0))
  for (s in e$spreads) {
    of <- unlist(s$of)
    if (!is.character(s$id) || length(of) != 2 || !is.character(s$label)) {
      stop("each export$spreads entry needs id, of (two variables) and label")
    }
    spreads <- rbind(spreads, data.frame(id = s$id, a = of[1], b = of[2], label = s$label))
  }
  if (anyDuplicated(spreads$id) || any(spreads$id %in% series$short_name)) {
    stop("spread ids must be distinct and differ from the variable names")
  }
  tr <- series$transform[match(c(spreads$a, spreads$b), series$short_name)]
  if (any(is.na(tr)) || any(tr != "raw")) stop("spreads must be differences of two rate (raw) variables")
  shown_in_levels <- as.character(unlist(e$shown_in_levels))
  tr <- series$transform[match(shown_in_levels, series$short_name)]
  if (any(is.na(tr)) || any(tr != "log")) {
    stop("export$shown_in_levels must list log series (rates are always shown in levels)")
  }
  list(run = e$run, horizon = as.integer(h), objects = objects, spreads = spreads,
       shown_in_levels = shown_in_levels)
}

#' The variables of the bundle: the model's series, then the spreads.
#'
#' Each variable has a display transform for the app: rates, indices kept
#' raw and spreads are shown as they are (`level`); log series listed in
#' export$shown_in_levels (prices and quantities read as levels, such as
#' the oil price in $/barrel or housing starts) in natural units
#' (`natural`); other log series follow the reader's growth setting
#' (`growth`). Every transform is exported for every log series anyway.
#'
#' @param series Series configuration.
#' @param spreads The spreads table of export_settings().
#' @param shown_in_levels Log series shown in natural units (export_settings()).
#' @return data.frame: variable, label, kind ("series" or "spread"),
#'   model_transform ("log" or "raw"; spreads are "raw"), display ("level",
#'   "natural" or "growth"), a and b (the spread's two variables; NA for a
#'   series).
export_variables <- function(series, spreads, shown_in_levels = character(0)) {
  display <- ifelse(series$transform == "raw", "level",
                    ifelse(series$short_name %in% shown_in_levels, "natural", "growth"))
  rbind(data.frame(variable = series$short_name, label = series$variable, kind = "series",
                   model_transform = series$transform, display = display,
                   a = NA_character_, b = NA_character_),
        data.frame(variable = spreads$id, label = spreads$label, kind = rep("spread", nrow(spreads)),
                   model_transform = rep("raw", nrow(spreads)), display = rep("level", nrow(spreads)),
                   a = spreads$a, b = spreads$b))
}

#' Display transforms of a variable.
#'
#' @param model_transform "log" or "raw".
#' @return Character vector of transforms.
display_transforms <- function(model_transform) {
  if (model_transform == "raw") "level" else c("level", "natural", "qoq_ann", "yoy")
}


# ---- Transforms ----

#' One variable (or spread) of an array of paths, as a matrix.
#'
#' @param x Array rows x H x n (paths or impulse-response draws), with the
#'   variable names as the third dimnames.
#' @param v Name of a variable of the bundle (export_variables()).
#' @param vars The table of export_variables().
#' @return Matrix rows x H in model units (a spread: a minus b).
variable_slice <- function(x, v, vars) {
  get <- function(name) matrix(x[, , name], dim(x)[1], dim(x)[2])
  k <- match(v, vars$variable)
  out <- if (vars$kind[k] == "spread") get(vars$a[k]) - get(vars$b[k]) else get(v)
  colnames(out) <- dimnames(x)[[2]]
  out
}

#' The same variable in the observed history (a vector over quarters).
#'
#' @param history Data matrix in model units (quarters x n, row names such
#'   as "2026Q2").
#' @param v,vars As in variable_slice().
#' @return Named numeric vector.
history_slice <- function(history, v, vars) {
  k <- match(v, vars$variable)
  out <- if (vars$kind[k] == "spread") history[, vars$a[k]] - history[, vars$b[k]] else history[, v]
  stats::setNames(out, rownames(history))
}

#' Display transform of paths in model units, with observed base quarters.
#'
#' Implements the transforms of the file header. Row j of the result is the
#' transform of path j: growth rates are computed path by path, so their
#' quantiles are quantiles of growth rates.
#'
#' @param y Matrix J x H of one variable in model units (paths, or per-draw
#'   means). For the history, a 1 x T matrix with `base` all NA.
#' @param base Values of the quarters just before y's first column, oldest
#'   first (the last four observed quarters for a forecast); NA where not
#'   available, which makes the growth rates that need them NA.
#' @param transform "level", "natural", "qoq_ann" or "yoy".
#' @return Matrix J x H, with y's column names.
transform_paths <- function(y, base, transform) {
  out <- switch(transform,
    level   = y,
    natural = exp(y / 100),
    qoq_ann = ,
    yoy = {
      lag <- if (transform == "yoy") 4L else 1L
      nb <- length(base)
      if (nb < lag) stop("the ", transform, " transform needs ", lag, " base quarter(s)")
      H <- ncol(y)
      full <- cbind(matrix(base, nrow(y), nb, byrow = TRUE), y)
      d <- (full[, nb + seq_len(H), drop = FALSE] - full[, nb + seq_len(H) - lag, drop = FALSE]) / 100
      if (transform == "yoy") 100 * (exp(d) - 1) else 100 * (exp(4 * d) - 1)
    },
    stop("unknown transform: ", transform))
  out <- matrix(out, nrow(y), ncol(y))
  colnames(out) <- colnames(y)
  out
}


# ---- Tables ----

#' First day of a quarter, for plotting ("2026Q3" -> 2026-07-01).
#'
#' @param q Quarter labels "YYYYQn".
#' @return Date vector.
quarter_start_date <- function(q) {
  as.Date(sprintf("%s-%02d-01", substr(q, 1, 4), 3L * as.integer(substr(q, 6, 6)) - 2L))
}

#' Quantiles of each column of a matrix, as rows of a table.
#'
#' @param values Matrix rows x H (paths or draws).
#' @param probs Probabilities (config: quantiles).
#' @param weights NULL, or one weight per row (tilting).
#' @return Matrix H x probs, columns q05, q16, ...
export_quantiles <- function(values, probs, weights = NULL) {
  q <- matrix(apply(values, 2, weighted_quantiles, weights, probs), ncol(values), length(probs), byrow = TRUE)
  colnames(q) <- sprintf("q%02d", round(100 * probs))
  q
}

#' Long table of quantiles for every variable and display transform.
#'
#' For each variable of the bundle and each of its transforms,
#' `values_of(v, tr)` gives a matrix (rows = paths, columns = quarters) and
#' its quantiles become rows of the table.
#'
#' @param object Id of the object (label column).
#' @param vars Table of export_variables().
#' @param values_of function(v, tr) returning a matrix rows x H.
#' @param quarters Labels of the H quarters.
#' @param probs Probabilities.
#' @param weights NULL or one weight per row.
#' @return data.frame: object, variable, transform, quarter, h, date, q05..q95.
quantile_table <- function(object, vars, values_of, quarters, probs, weights = NULL) {
  rows <- list()
  for (k in seq_len(nrow(vars))) {
    v <- vars$variable[k]
    for (tr in display_transforms(vars$model_transform[k])) {
      q <- export_quantiles(values_of(v, tr), probs, weights)
      rows[[length(rows) + 1]] <- data.frame(object = object, variable = v, transform = tr, quarter = quarters,
                                             h = seq_along(quarters), date = quarter_start_date(quarters), q)
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Forecast quantiles of one object, over its paths.
#'
#' @param paths Array J x H x n of paths in model units (the first H
#'   quarters are used).
#' @param history Observed data matrix (base quarters).
#' @param object,vars,quarters,probs As in quantile_table().
#' @param weights NULL, or the tilting weights of the paths.
#' @return Table of quantile_table().
export_forecast <- function(paths, history, object, vars, quarters, probs, weights = NULL) {
  H <- length(quarters)
  values_of <- function(v, tr) {
    y <- variable_slice(paths, v, vars)[, seq_len(H), drop = FALSE]
    transform_paths(y, utils::tail(history_slice(history, v, vars), 4), tr)
  }
  quantile_table(object, vars, values_of, quarters, probs, weights)
}

#' The observed history in every display transform.
#'
#' Quarters whose growth rate needs an earlier quarter than the sample's
#' first are left out.
#'
#' @param history Observed data matrix (quarters x n).
#' @param vars Table of export_variables().
#' @return data.frame: variable, transform, quarter, date, value.
export_history <- function(history, vars) {
  rows <- list()
  for (k in seq_len(nrow(vars))) {
    v <- vars$variable[k]
    y <- matrix(history_slice(history, v, vars), 1, dimnames = list(NULL, rownames(history)))
    for (tr in display_transforms(vars$model_transform[k])) {
      value <- as.vector(transform_paths(y, rep(NA_real_, 4), tr))
      keep <- !is.na(value)
      rows[[length(rows) + 1]] <- data.frame(variable = v, transform = tr, quarter = rownames(history)[keep],
                                             date = quarter_start_date(rownames(history)[keep]),
                                             value = value[keep])
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Impulse responses: quantiles over draws.
#'
#' Transforms: level (the response in model units per unit shock: log
#' points for log series, percentage points for rates); for log series also
#' yoy, the response at h minus the response at h - 4 (zero before the
#' shock), the response of the 4-quarter log change in percentage points.
#' Spreads respond by the difference of their two yields' responses.
#'
#' @param irf_draws Named list of arrays draws x (H + 1) x n (irf_draws(),
#'   names "girf_<shock>" or "recursive_<shock>").
#' @param vars Table of export_variables().
#' @param probs Probabilities.
#' @return data.frame: identification, scheme ("girf" or "recursive"),
#'   shock, variable, transform, horizon (0..H), q05..q95.
export_irf <- function(irf_draws, vars, probs) {
  rows <- list()
  for (id in names(irf_draws)) {
    x <- irf_draws[[id]]
    horizon <- seq_len(dim(x)[2]) - 1L
    scheme <- sub("_.*$", "", id)
    shock <- sub("^[a-z]+_", "", id)
    for (k in seq_len(nrow(vars))) {
      v <- vars$variable[k]
      r <- variable_slice(x, v, vars)
      trs <- if (vars$model_transform[k] == "log") c("level", "yoy") else "level"
      for (tr in trs) {
        val <- if (tr == "level") r else r - cbind(matrix(0, nrow(r), 4), r)[, seq_len(ncol(r)), drop = FALSE]
        rows[[length(rows) + 1]] <- data.frame(identification = id, scheme = scheme, shock = shock, variable = v,
                                               transform = tr, horizon = horizon, export_quantiles(val, probs))
      }
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


# ---- The bundle ----

#' Build the bundle from the engine's outputs (in memory).
#'
#' @param baseline The engine's baseline (engine_baseline(), with
#'   inputs$quarters).
#' @param scenarios Named list of the engine's scenarios (engine_scenario()),
#'   named by id, in the order to show.
#' @param irf The engine's impulse responses (engine_irfs()).
#' @param irf_extra NULL, or a table of export_irf() for further shocks
#'   (run_engine_irf_extra()), appended to the responses of `irf`.
#' @param history Observed data matrix (the engine's origin data).
#' @param series Series configuration.
#' @param index Scenario index (read_scenario_index()), for labels.
#' @param settings Output of export_settings().
#' @param probs Probabilities (config: quantiles).
#' @return List: variables, objects, history, forecasts, conditions (fixed
#'   cells), targets (soft targets and their check), irf.
build_bundle <- function(baseline, scenarios, irf, history, series, index, settings, probs, irf_extra = NULL) {
  vars <- export_variables(series, settings$spreads, settings$shown_in_levels)
  quarters <- baseline$inputs$quarters[seq_len(settings$horizon)]
  if (rownames(history)[nrow(history)] != baseline$inputs$last_observed) {
    stop("the history ends in ", rownames(history)[nrow(history)], "; the forecast starts after ",
         baseline$inputs$last_observed)
  }
  objects <- data.frame(id = "baseline", label = "Baseline", short_label = "Unconditional",
                        description = "The model's unconditional forecast.",
                        type = "baseline", fixed_cells = 0L, ess_share = NA_real_)
  forecasts <- list(export_forecast(baseline$paths, history, "baseline", vars, quarters, probs))
  conditions <- list(); targets <- list()
  for (id in names(scenarios)) {
    s <- scenarios[[id]]
    tilted <- !is.null(s$tilt)
    k <- match(id, index$id)
    if (is.na(k)) stop("scenario ", id, " is not in the scenario index")
    short <- if (is.null(index$short_label)) index$label[k] else index$short_label[k]
    objects <- rbind(objects, data.frame(id = id, label = index$label[k], short_label = short,
                                         description = index$description[k],
                                         type = if (tilted) "tilted" else "hard", fixed_cells = s$fixed_cells,
                                         ess_share = if (tilted) s$tilt$ess_share else NA_real_))
    w <- if (tilted) s$tilt$weights else NULL
    forecasts[[id]] <- export_forecast(s$paths, history, id, vars, quarters, probs, w)
    if (tilted) {
      ck <- s$tilt$check
      targets[[id]] <- data.frame(object = id, variable = ck$variable, quarter = ck$quarter, h = ck$h,
                                  measure = ck$measure, statistic = ck$statistic, target = ck$target,
                                  untilted_median = ck$untilted_median, goal = ck$goal, achieved = ck$achieved)
    }
    fixed <- which(!is.na(s$conditions[seq_len(settings$horizon), , drop = FALSE]), arr.ind = TRUE)
    if (nrow(fixed)) {
      v <- colnames(s$conditions)[fixed[, 2]]
      value <- s$conditions[fixed]
      log_v <- series$transform[match(v, series$short_name)] == "log"
      conditions[[id]] <- rbind(
        data.frame(object = id, variable = v, transform = "level", quarter = quarters[fixed[, 1]],
                   h = fixed[, 1], value = value),
        data.frame(object = id, variable = v[log_v], transform = rep("natural", sum(log_v)),
                   quarter = quarters[fixed[log_v, 1]], h = fixed[log_v, 1], value = exp(value[log_v] / 100)))
    }
  }
  bind <- function(x, empty) {
    if (!length(x)) return(empty)
    out <- do.call(rbind, unname(x))
    rownames(out) <- NULL
    out
  }
  list(variables = vars, objects = objects,
       history = export_history(history, vars),
       forecasts = bind(forecasts, NULL),
       conditions = bind(conditions, data.frame(object = character(0), variable = character(0),
                                                transform = character(0), quarter = character(0),
                                                h = integer(0), value = numeric(0))),
       targets = bind(targets, data.frame(object = character(0), variable = character(0))),
       irf = rbind(export_irf(irf$draws, vars, probs), irf_extra))
}

#' Metadata of the bundle: what it holds, from which run and code.
#'
#' @param cfg Parsed config/run.yml.
#' @param settings Output of export_settings().
#' @param engine_meta Parsed <run_dir>/engine/meta.json.
#' @param engine_meta_file Its path (for the MD5).
#' @param history Observed data matrix.
#' @param bundle Output of build_bundle().
#' @return List (written to meta.json and kept in the bundle).
export_meta <- function(cfg, settings, engine_meta, engine_meta_file, history, bundle) {
  lock <- "renv.lock"
  pkgs <- c("targets", "yaml", "jsonlite")
  list(
    run_id = engine_meta$run_id, created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    data = list(qd_vintage_file = cfg$data$qd_vintage_file, api_realtime_date = cfg$data$api_realtime_date,
                first_quarter = rownames(history)[1], last_observed = rownames(history)[nrow(history)]),
    model = list(lags = cfg$model$lags, covid_method = settings$run,
                 covid_window = if (settings$run == "exclude_window") unlist(cfg$covid$window) else NULL),
    forecast = list(quarters = unique(bundle$forecasts$quarter)[c(1, settings$horizon)],
                    horizon_exported = settings$horizon, horizon_computed = engine_meta$settings$horizon,
                    n_draws = engine_meta$n_draws, paths_per_draw = engine_meta$settings$paths_per_draw,
                    tilting_paths_per_draw = engine_meta$settings$tilting$paths_per_draw,
                    conditional_method = engine_meta$settings$method, seed = engine_meta$settings$seed,
                    quantiles = unlist(cfg$quantiles)),
    objects = bundle$objects[c("id", "label", "short_label", "type", "fixed_cells", "ess_share")],
    transforms = list(
      level = "model units: 100 ln(x) for log series, natural units for rates",
      natural = "x = exp(level / 100), log series only",
      qoq_ann = "100 [(x_t / x_{t-1})^4 - 1], percent change at an annual rate",
      yoy = "100 (x_t / x_{t-4} - 1), percent change from a year earlier",
      irf_yoy = "response at h minus response at h - 4, percentage points"),
    differences = paste("not exported: the app shows the difference between two objects' medians,",
                        "in the same transform and quarter"),
    spreads = settings$spreads,
    engine = list(meta_md5 = unname(tools::md5sum(engine_meta_file)), created = engine_meta$created,
                  draws_md5 = engine_meta$draws_md5),
    code = git_state(), r_version = R.version.string,
    packages = as.list(stats::setNames(vapply(pkgs, function(p) as.character(utils::packageVersion(p)), ""), pkgs)),
    renv_lock_md5 = if (file.exists(lock)) unname(tools::md5sum(lock)) else NA_character_)
}

#' Export the published run: read its engine files, build and write the
#' bundle.
#'
#' Checks that the engine files read are the ones engine/meta.json records
#' (same MD5), so the bundle never mixes outputs of two engine runs.
#'
#' @param run_dir Folder of the estimation run.
#' @param engine_meta_file Path of <run_dir>/engine/meta.json.
#' @param series Series configuration.
#' @param index Scenario index (read_scenario_index()).
#' @param cfg Parsed config/run.yml.
#' @param irf_extra_file NULL, or the path written by run_engine_irf_extra().
#' @return Paths of bundle.rds and meta.json.
run_export <- function(run_dir, engine_meta_file, series, index, cfg, irf_extra_file = NULL) {
  settings <- export_settings(cfg, series)
  eng <- file.path(run_dir, "engine")
  em <- jsonlite::read_json(engine_meta_file, simplifyVector = TRUE)
  files <- c("baseline.rds", "irf.rds", paste0("scenario_", settings$objects, ".rds"))
  for (f in files) {
    path <- file.path(eng, f)
    if (!file.exists(path)) stop("engine file not found: ", path)
    if (!identical(unname(tools::md5sum(path)), em$files[[f]])) {
      stop(f, " is not the file recorded in engine/meta.json; rerun the engine")
    }
  }
  baseline <- readRDS(file.path(eng, "baseline.rds"))
  scenarios <- lapply(stats::setNames(settings$objects, settings$objects),
                      function(id) readRDS(file.path(eng, paste0("scenario_", id, ".rds"))))
  irf <- readRDS(file.path(eng, "irf.rds"))
  history <- est_data_to_matrix(readRDS(file.path(run_dir, "covid.rds"))$origin_data, series)
  # Generalized IRFs of the other variables (quantiles only), if computed:
  # they must come from the same draws as the rest of the engine.
  extra <- NULL
  if (!is.null(irf_extra_file)) {
    x <- readRDS(irf_extra_file)
    if (!identical(x$draws_md5, em$draws_md5)) stop("irf_extra.rds was computed from other draws; rerun it")
    extra <- x$table
  }

  bundle <- build_bundle(baseline, scenarios, irf, history, series, index, settings, unlist(cfg$quantiles), extra)
  bundle$meta <- export_meta(cfg, settings, em, engine_meta_file, history, bundle)
  bundle$meta$irf <- list(identifications = unique(bundle$irf$identification),
                          draws_stored = names(irf$draws),
                          irf_extra_md5 = if (is.null(irf_extra_file)) NULL else unname(tools::md5sum(irf_extra_file)))
  out_dir <- file.path(run_dir, "bundle")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  paths <- file.path(out_dir, c("bundle.rds", "meta.json"))
  saveRDS(bundle, paths[1])
  jsonlite::write_json(bundle$meta, paths[2], auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
  paths
}
