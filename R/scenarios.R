# scenarios.R: read, check and build scenario files for the conditional
# forecast engine.
#
# Why it exists: a scenario is a conditional forecast. The user fixes future
# values of some variables in some quarters, and the model fills in every
# other cell (Crump et al. 2021, Section 2.3). The engine
# (forecast_conditional() in R/forecast_conditional.R) takes those fixed values as a
# conditions matrix, H x n, with NA for free cells. This file is the layer
# between the CSV a person writes and that matrix: it reads the scenario
# index and files, checks them, and turns them into conditions. It also
# provides helpers to write scenario paths in the model's units.
#
# Scenario files (docs/data.md for units):
#   scenarios/scenarios.yml   the index: one entry per scenario with id,
#                             label and description, and optionally `soft`,
#                             the columns that hold soft targets (entropic
#                             tilting, R/tilting.R) instead of fixed values;
#   scenarios/<id>.csv        one row per quarter: a column `quarter`
#                             (e.g. 2027Q1), then any subset of the model's
#                             variables (short names). Values are absolute
#                             levels in model units: 100 ln(x) for log
#                             series, natural units for rates. A blank cell
#                             (or NA) is free. Rows may skip quarters and
#                             columns may come in any order.
#
# Quarters are calendar quarters, so a scenario belongs to a data vintage:
# it must start after the last observed quarter of the data it is run on.
# A file whose quarters are no longer in the future is rejected, not
# shifted, because a scenario written against old data should be revisited.
#
# Inputs:  scenarios/scenarios.yml, scenarios/<id>.csv, the variable names
#          and transforms (config/series.csv), the last observed quarter and
#          the forecast horizon H.
# Outputs: conditions matrices H x n (NA = free) for forecast_conditional();
#          tables of soft targets for tilting;
#          scenario paths in model units; a template CSV of the baseline.


# ---- Quarters ----

#' Quarter label as a count of quarters (year * 4 + quarter - 1).
#'
#' Counting quarters makes "how many quarters after the last observed one"
#' a subtraction.
#'
#' @param x Character vector of labels such as "2027Q1".
#' @return Integer vector; NA where the label is not a valid quarter.
quarter_index <- function(x) {
  x  <- trimws(as.character(x))
  ok <- !is.na(x) & grepl("^[0-9]{4}Q[1-4]$", x)
  out <- rep(NA_integer_, length(x))
  out[ok] <- as.integer(substr(x[ok], 1, 4)) * 4L + as.integer(substr(x[ok], 6, 6)) - 1L
  out
}

#' Labels of the forecast quarters T+1, ..., T+H.
#'
#' @param last_observed Label of the last observed quarter T (e.g. "2026Q2").
#' @param horizon Number of forecast quarters H.
#' @return Character vector of H quarter labels.
forecast_quarters <- function(last_observed, horizon) {
  t0 <- quarter_index(last_observed)
  if (is.na(t0)) stop("last_observed is not a quarter label (expected e.g. 2026Q2): ", last_observed)
  k <- t0 + seq_len(horizon)
  paste0(k %/% 4L, "Q", k %% 4L + 1L)
}


# ---- Scenario index ----

#' Read the scenario index, scenarios/scenarios.yml.
#'
#' The file holds a list `scenarios`, one entry per scenario with the fields
#' id (also the CSV file name, scenarios/<id>.csv), label (a short title)
#' and description (what the scenario assumes and where its numbers come
#' from), and optionally `short_label` (the name shown in the website's
#' scenario list; the label if absent) and `soft`: the variables whose
#' columns are soft conditions (see parse_soft()). Any other field is an
#' error: a field the engine does not read would otherwise be ignored
#' without notice.
#'
#' @param path Path to scenarios.yml.
#' @return data.frame with columns id, label, short_label, description and
#'   soft (a list column: for each scenario the table of parse_soft(), zero
#'   rows if it has no soft conditions); zero rows if the index is empty.
#'   Stops if a field is missing or unknown, a soft entry is malformed, an
#'   id is repeated or an id is not a plain file name.
read_scenario_index <- function(path) {
  idx <- yaml::read_yaml(path)
  entries <- idx$scenarios
  if (is.null(entries)) entries <- list()
  fields <- c("id", "label", "description")
  soft <- vector("list", length(entries))
  for (k in seq_along(entries)) {
    e <- entries[[k]]
    miss <- fields[!vapply(fields, function(f) is.character(e[[f]]) && length(e[[f]]) == 1 && nzchar(e[[f]]),
                           logical(1))]
    if (length(miss) > 0) stop("scenario entry ", k, " in ", path, " lacks: ", paste(miss, collapse = ", "))
    extra <- setdiff(names(e), c(fields, "short_label", "soft"))
    if (length(extra) > 0) stop("scenario entry ", k, " in ", path, " has unknown field(s): ",
                                paste(extra, collapse = ", "))
    soft[[k]] <- tryCatch(parse_soft(e$soft),
                          error = function(err) stop("scenario entry ", k, " in ", path, ": ",
                                                     conditionMessage(err), call. = FALSE))
  }
  out <- data.frame(id = vapply(entries, function(e) e$id, ""),
                    label = vapply(entries, function(e) e$label, ""),
                    short_label = vapply(entries, function(e) {
                      s <- e$short_label
                      if (is.character(s) && length(s) == 1 && nzchar(s)) s else e$label
                    }, ""),
                    description = vapply(entries, function(e) trimws(e$description), ""),
                    stringsAsFactors = FALSE)
  out$soft <- soft
  # The id names the CSV file, so it must be a safe file name.
  bad <- out$id[!grepl("^[a-z0-9_]+$", out$id)]
  if (length(bad) > 0) stop("scenario ids must use lower-case letters, digits and _: ", paste(bad, collapse = ", "))
  dup <- unique(out$id[duplicated(out$id)])
  if (length(dup) > 0) stop("repeated scenario ids: ", paste(dup, collapse = ", "))
  out
}

#' Parse the `soft` field of a scenario entry.
#'
#' A soft column holds targets for the forecast distribution instead of
#' fixed values; they are imposed by entropic tilting (R/tilting.R; Crump
#' et al. 2021, Section 4.2). Each soft variable has a measure, the form in
#' which its column's values are written, and a statistic:
#'   short form   soft: {pce: yoy}                     (statistic: median)
#'   long form    soft: {pce: {measure: yoy, statistic: mean}}
#' Measures: level (model units), yoy, qoq_ann (exact percent changes, as
#' path_from_growth()), annual_avg (growth of the four-quarter average over
#' the previous four; in a Q4 row, calendar-year average growth; R/tilting.R).
#' Statistics: median (default) or mean.
#'
#' @param soft The `soft` field as read by yaml (NULL if absent).
#' @return data.frame variable, measure, statistic (zero rows if NULL).
parse_soft <- function(soft) {
  empty <- data.frame(variable = character(0), measure = character(0), statistic = character(0),
                      stringsAsFactors = FALSE)
  if (is.null(soft)) return(empty)
  if (!is.list(soft) || is.null(names(soft)) || any(!nzchar(names(soft)))) {
    stop("soft must map variable names to a measure, e.g. soft: {pce: yoy}")
  }
  rows <- lapply(names(soft), function(v) {
    s <- soft[[v]]
    if (is.character(s) && length(s) == 1) s <- list(measure = s)
    if (!is.list(s)) stop("soft entry for ", v, " must be a measure or {measure: ..., statistic: ...}")
    unknown <- setdiff(names(s), c("measure", "statistic"))
    if (length(unknown)) stop("soft entry for ", v, " has unknown field(s): ", paste(unknown, collapse = ", "))
    measure <- s$measure
    statistic <- if (is.null(s$statistic)) "median" else s$statistic
    if (!(is.character(measure) && length(measure) == 1 &&
          measure %in% c("level", "yoy", "qoq_ann", "annual_avg"))) {
      stop("soft entry for ", v, ": measure must be level, yoy, qoq_ann or annual_avg")
    }
    if (!(is.character(statistic) && length(statistic) == 1 && statistic %in% c("median", "mean"))) {
      stop("soft entry for ", v, ": statistic must be median or mean")
    }
    data.frame(variable = v, measure = measure, statistic = statistic, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  dup <- unique(out$variable[duplicated(out$variable)])
  if (length(dup)) stop("soft lists a variable more than once: ", paste(dup, collapse = ", "))
  out
}

#' The soft conditions of one scenario in the index.
#'
#' @param index Output of read_scenario_index().
#' @param id Scenario id.
#' @return data.frame variable, measure, statistic (zero rows if none).
scenario_soft <- function(index, id) {
  k <- match(id, index$id)
  if (is.na(k)) stop("scenario not in the index: ", id)
  index$soft[[k]]
}


# ---- Scenario files ----

#' Check a scenario table against the rules of a scenario file.
#'
#' All problems are collected and reported together, so a file can be fixed
#' in one pass. The rules: a `quarter` column; every other column a model
#' variable, at most once; quarters valid, distinct, after the last observed
#' quarter and within the horizon; every non-blank cell a finite number; in
#' no quarter are all variables fixed (the model would have nothing to fill
#' in, and the conditioning has no free cell to absorb the conditions).
#'
#' Soft columns hold targets, not fixed values, so they do not count
#' towards "every variable fixed".
#'
#' @param df data.frame read with all columns as character (see
#'   read_scenario()).
#' @param variables Model variable names (short names, series order).
#' @param last_observed Label of the last observed quarter T.
#' @param horizon Forecast horizon H.
#' @param soft_vars Names of the soft columns (parse_soft()).
#' @return TRUE invisibly; stops with the list of problems otherwise.
validate_scenario <- function(df, variables, last_observed, horizon, soft_vars = character(0)) {
  problems <- character(0)
  cols <- names(df)
  if (!"quarter" %in% cols) problems <- c(problems, "no `quarter` column")
  vars_in <- setdiff(cols, "quarter")
  unknown <- setdiff(vars_in, variables)
  if (length(unknown) > 0) problems <- c(problems, paste0("unknown column(s): ", paste(unknown, collapse = ", ")))
  dup_cols <- unique(cols[duplicated(cols)])
  if (length(dup_cols) > 0) problems <- c(problems, paste0("repeated column(s): ", paste(dup_cols, collapse = ", ")))

  if ("quarter" %in% cols) {
    q <- trimws(df$quarter)
    qi <- quarter_index(q)
    if (any(is.na(qi))) {
      problems <- c(problems, paste0("not a quarter label (expected e.g. 2027Q1): ",
                                     paste(q[is.na(qi)], collapse = ", ")))
    }
    dup_q <- unique(q[!is.na(qi) & duplicated(qi)])
    if (length(dup_q) > 0) problems <- c(problems, paste0("repeated quarter(s): ", paste(dup_q, collapse = ", ")))
    h <- qi - quarter_index(last_observed)
    out <- !is.na(h) & (h < 1 | h > horizon)
    if (any(out)) {
      rng <- forecast_quarters(last_observed, horizon)[c(1, horizon)]
      problems <- c(problems, paste0("quarter(s) outside the forecast range ", rng[1], " to ", rng[2],
                                     " (data end in ", last_observed, "): ", paste(q[out], collapse = ", ")))
    }
  }

  vals <- df[, intersect(vars_in, variables), drop = FALSE]
  for (v in names(vals)) {
    x <- trimws(vals[[v]])
    filled <- !is.na(x) & nzchar(x) & x != "NA"
    num <- suppressWarnings(as.numeric(x[filled]))
    bad <- !is.finite(num)
    if (any(bad)) {
      problems <- c(problems, paste0("non-numeric value(s) in ", v, ": ", paste(x[filled][bad], collapse = ", ")))
    }
  }
  absent_soft <- setdiff(soft_vars, vars_in)
  if (length(absent_soft) > 0) {
    problems <- c(problems, paste0("soft variable(s) without a column: ", paste(absent_soft, collapse = ", ")))
  }
  if (length(variables) > 0 && all(variables %in% names(vals)) && nrow(vals) > 0 &&
      length(intersect(soft_vars, variables)) == 0) {
    filled_all <- apply(vals[, variables, drop = FALSE], 1, function(r) {
      r <- trimws(r)
      all(!is.na(r) & nzchar(r) & r != "NA")
    })
    if (any(filled_all)) {
      problems <- c(problems, paste0("every variable is fixed in quarter(s): ",
                                     paste(trimws(df$quarter)[filled_all], collapse = ", ")))
    }
  }

  if (length(problems) > 0) stop("invalid scenario:\n  - ", paste(problems, collapse = "\n  - "), call. = FALSE)
  invisible(TRUE)
}

#' Turn a checked scenario table into the engine's conditions matrix.
#'
#' Row h of the result is quarter T + h; columns are the model variables in
#' series order, whatever the column order of the file, so the conditions
#' (and the forecast) do not depend on how the file is laid out.
#'
#' @param df Scenario table that passed validate_scenario().
#' @param variables Model variable names (series order).
#' @param last_observed Label of the last observed quarter T.
#' @param horizon Forecast horizon H.
#' @return Numeric matrix H x n, NA = free, with row names the forecast
#'   quarters and column names the variables.
scenario_conditions <- function(df, variables, last_observed, horizon) {
  out <- matrix(NA_real_, horizon, length(variables),
                dimnames = list(forecast_quarters(last_observed, horizon), variables))
  h <- quarter_index(df$quarter) - quarter_index(last_observed)
  for (v in intersect(names(df), variables)) {
    x <- trimws(df[[v]])
    x[is.na(x) | !nzchar(x) | x == "NA"] <- NA
    out[h, v] <- as.numeric(x)
  }
  out
}

#' Read and check a scenario CSV (shared by read_scenario() and
#' soft_targets()).
#'
#' @param path,variables,last_observed,horizon As in read_scenario().
#' @param soft_vars Names of the soft columns.
#' @return The checked table, all columns as text.
read_scenario_table <- function(path, variables, last_observed, horizon, soft_vars = character(0)) {
  # Every column is read as text so that a stray word is reported as a
  # non-numeric value rather than silently turning a column into text.
  df <- utils::read.csv(path, colClasses = "character", na.strings = character(0),
                        check.names = FALSE, strip.white = TRUE)
  tryCatch(validate_scenario(df, variables, last_observed, horizon, soft_vars),
           error = function(e) stop(basename(path), ": ", conditionMessage(e), call. = FALSE))
  df
}

#' Read a scenario CSV and return its conditions matrix (hard cells only).
#'
#' @param path Path to scenarios/<id>.csv.
#' @param variables Model variable names (series order).
#' @param last_observed Label of the last observed quarter T.
#' @param horizon Forecast horizon H.
#' @param soft Soft conditions of the scenario (scenario_soft()); their
#'   columns hold targets, not fixed values, and are left free here.
#' @return Numeric matrix H x n for forecast_conditional() (NA = free).
read_scenario <- function(path, variables, last_observed, horizon, soft = NULL) {
  soft_vars <- if (is.null(soft)) character(0) else soft$variable
  df <- read_scenario_table(path, variables, last_observed, horizon, soft_vars)
  conditions <- scenario_conditions(df, variables, last_observed, horizon)
  conditions[, intersect(soft_vars, variables)] <- NA
  conditions
}

#' Read the soft targets of a scenario.
#'
#' Every non-blank cell of a soft column is a target, written in the
#' column's measure (parse_soft()). A growth measure applies only to log
#' series, since a percent change of a rate is not what a user means.
#'
#' @param path,variables,last_observed,horizon As in read_scenario().
#' @param soft Soft conditions of the scenario (scenario_soft()).
#' @param series Series configuration (read_series_config()).
#' @return data.frame variable, quarter, h (1 = T + 1), measure, statistic,
#'   target; zero rows if the scenario has no soft conditions. Stops if a
#'   growth measure is given for a series that is not in logs, or a soft
#'   column has no target.
soft_targets <- function(path, variables, last_observed, horizon, soft, series) {
  out <- data.frame(variable = character(0), quarter = character(0), h = integer(0),
                    measure = character(0), statistic = character(0), target = numeric(0),
                    stringsAsFactors = FALSE)
  if (is.null(soft) || nrow(soft) == 0) return(out)
  unknown <- setdiff(soft$variable, variables)
  if (length(unknown)) stop(basename(path), ": soft variable(s) not in the model: ", paste(unknown, collapse = ", "))
  tr <- series$transform[match(soft$variable, series$short_name)]
  bad <- soft$variable[soft$measure != "level" & tr != "log"]
  if (length(bad)) stop(basename(path), ": growth measures apply only to log series, not to ",
                        paste(bad, collapse = ", "))
  df <- read_scenario_table(path, variables, last_observed, horizon, soft$variable)
  h_row <- quarter_index(df$quarter) - quarter_index(last_observed)
  rows <- list()
  for (k in seq_len(nrow(soft))) {
    x <- trimws(df[[soft$variable[k]]])
    filled <- !is.na(x) & nzchar(x) & x != "NA"
    if (!any(filled)) stop(basename(path), ": soft column ", soft$variable[k], " has no target")
    rows[[k]] <- data.frame(variable = soft$variable[k], quarter = trimws(df$quarter[filled]),
                            h = as.integer(h_row[filled]), measure = soft$measure[k],
                            statistic = soft$statistic[k], target = as.numeric(x[filled]),
                            stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  out[order(out$h, match(out$variable, variables)), , drop = FALSE]
}


# ---- Helpers to write scenario paths ----

#' Convert values in natural units to model units.
#'
#' Log series are stored as 100 ln(x) (a level in log points, docs/data.md
#' section 2); rates, ratios and indices marked `raw` are used as published.
#' So a real GDP level of 24,500 (billions of chained dollars) becomes
#' 100 ln(24500) = 1010.6, and an unemployment rate of 4.3 stays 4.3.
#'
#' @param var Short name of the variable.
#' @param x Numeric vector in the variable's natural units.
#' @param series Series configuration (read_series_config()).
#' @return Numeric vector in model units.
to_model_units <- function(var, x, series) {
  tr <- series$transform[match(var, series$short_name)]
  if (is.na(tr)) stop("unknown variable: ", var)
  if (tr == "raw") return(x)
  if (any(!is.na(x) & x <= 0)) stop(var, " is a log series: values must be positive")
  100 * log(x)
}

#' Build a path in model units from growth rates and the observed data.
#'
#' For log series only. Growth rates are percent changes as published:
#'   qoq_ann  percent change at an annual rate from the previous quarter,
#'            g = 100 [(x_t / x_{t-1})^4 - 1], the convention of the BEA's
#'            national accounts;
#'   yoy      percent change from the same quarter a year earlier,
#'            g = 100 (x_t / x_{t-4} - 1).
#' In model units y = 100 ln(x), these invert exactly to
#'   qoq_ann  y_t = y_{t-1} + 100 ln(1 + g/100) / 4,
#'   yoy      y_t = y_{t-4} + 100 ln(1 + g/100).
#' The base values come from the last observed quarters of `history`, then
#' from the path itself. An NA growth leaves that quarter free; a later
#' growth rate whose base quarter is free cannot be turned into a level and
#' is an error.
#'
#' @param var Short name of a log series.
#' @param growth Growth rates in percent for quarters T+1, T+2, ... (NA =
#'   free).
#' @param type "qoq_ann" or "yoy".
#' @param history Data matrix in model units, rows = quarters (row names
#'   such as "2026Q2"), columns = variables; its last rows are the base.
#' @param series Series configuration (read_series_config()).
#' @return Named numeric vector (forecast quarters) of levels in model units,
#'   NA where free.
path_from_growth <- function(var, growth, type, history, series) {
  tr <- series$transform[match(var, series$short_name)]
  if (is.na(tr)) stop("unknown variable: ", var)
  if (tr != "log") stop(var, " is not a log series: write its path in natural units instead")
  lag <- switch(type, qoq_ann = 1L, yoy = 4L, stop("type must be \"qoq_ann\" or \"yoy\""))
  step <- function(g) if (type == "qoq_ann") 100 * log(1 + g / 100) / 4 else 100 * log(1 + g / 100)

  T_ <- nrow(history)
  y <- c(history[(T_ - lag + 1):T_, var], rep(NA_real_, length(growth)))   # base quarters, then the path
  for (h in seq_along(growth)) {
    if (is.na(growth[h])) next
    base <- y[h]                                     # quarter T + h - lag
    if (is.na(base)) {
      stop(var, ": the growth rate for quarter ", h, " needs the level of quarter ", h - lag,
           ", which is free")
    }
    y[h + lag] <- base + step(growth[h])
  }
  stats::setNames(y[-seq_len(lag)], forecast_quarters(rownames(history)[T_], length(growth)))
}

#' Write a scenario template: the baseline median path of every variable.
#'
#' The baseline is the unconditional forecast (forecast_unconditional(),
#' with shocks) from the run's posterior draws and its origin data; the
#' template holds its median path over draws in model units. A user copies
#' the file, keeps the quarters and columns to fix, and edits the values.
#'
#' @param run_dir Folder of an estimation run (draws.rds, covid.rds).
#' @param series Series configuration (read_series_config()).
#' @param lags,horizon Lags p and forecast horizon H.
#' @param seed Seed of the random number generator (forecast shocks).
#' @param path File to write (NULL: do not write).
#' @param digits Decimal places written.
#' @return data.frame with `quarter` and one column per variable (invisibly).
make_scenario_template <- function(run_dir, series, lags, horizon, seed, path = NULL, digits = 4) {
  draws <- readRDS(file.path(run_dir, "draws.rds"))
  covid <- readRDS(file.path(run_dir, "covid.rds"))
  history <- est_data_to_matrix(covid$origin_data, series)
  if (!identical(draws$variables, colnames(history))) stop("variable order of draws and history differ")
  fc <- forecast_unconditional(draws$beta, draws$sigma, history, lags, horizon, seed = seed)
  med <- apply(fc, c(2, 3), stats::median)
  out <- data.frame(quarter = forecast_quarters(rownames(history)[nrow(history)], horizon),
                    round(med, digits), check.names = FALSE)
  names(out)[-1] <- colnames(history)
  if (!is.null(path)) utils::write.csv(out, path, row.names = FALSE)
  invisible(out)
}
