# data_assemble.R: build the 31 quarterly model series and the two data tables
# of pipeline stage 1 (see _targets.R).
#
# Why it exists: the loaders (data_load.R) and the constructed series (at the
# end of this file) each handle one source. This file applies them to the series listed in
# config/series.csv, lines the 31 series up on a common quarterly index, and
# applies the transformations of docs/data.md (section 1):
#   levels     : natural units (billions of chained dollars, index points,
#                percent, dollars per barrel, ...), one column per variable;
#   model_data : what the VAR sees, 100 * ln(x) for "log" series (a level in
#                log points, not a growth rate) and x unchanged for "raw" series.
#
# Inputs:  config/series.csv, the data block of config/run.yml, raw files.
# Outputs: data.frames returned to the caller (scripts/build_data.R saves them).


#' Path of a raw file named in the run config.
#'
#' @param cfg Parsed config/run.yml.
#' @param relative Path relative to cfg$data$raw_dir.
#' @return Full path. Stops if the file does not exist.
raw_path <- function(cfg, relative) {
  p <- file.path(cfg$data$raw_dir, relative)
  if (!file.exists(p)) stop("raw file not found: ", p)
  p
}

#' Build each model variable as a quarterly series in natural units.
#'
#' Source rules (docs/data.md, sections 2 and 3):
#'   QD     : column of the QD vintage file named in the config.
#'   API    : ALFRED vintage history read as of cfg$data$api_realtime_date;
#'            monthly series averaged to quarterly (whole quarters only).
#'   CALC   : rvol (3.4) from daily closes.
#'   MANUAL : gsci_nonenergy (3.3), World Bank monthly index averaged to quarterly.
#' A series listed under cfg$data$backfill is then extended back before its
#' first observation with the named substitute (docs/data.md 3.6).
#'
#' @param series Series config from read_series_config().
#' @param cfg Parsed config/run.yml.
#' @return Named list (names = short names, in config order) of data.frames
#'   with columns quarter (Date) and value (numeric, natural units). Attribute
#'   "backfill" lists every backfilled quarter (series, quarter, value, source)
#'   for audit.
build_series_list <- function(series, cfg) {
  qd <- read_qd_vintage(raw_path(cfg, cfg$data$qd_vintage_file), series$id[series$source == "QD"])
  as_of <- as.Date(cfg$data$api_realtime_date)

  out <- list()
  for (i in seq_len(nrow(series))) {
    s <- series[i, ]

    if (s$source == "QD") {
      x <- data.frame(quarter = qd$quarter, value = qd[[s$id]])

    } else if (s$source == "API") {
      file <- cfg$data$alfred_files[[s$id]]
      if (is.null(file)) stop("no ALFRED file in config for ", s$id)
      v <- alfred_as_of(read_alfred_vintages(raw_path(cfg, file)), as_of)
      if (s$api_aggregation == "avg") {
        x <- monthly_to_quarterly(v)
      } else {
        if (any(v$date != quarter_start(v$date))) stop(s$id, " is not quarterly")
        x <- data.frame(quarter = v$date, value = v$value)
      }

    } else if (s$short_name == "gsci_nonenergy") {
      x <- monthly_to_quarterly(read_worldbank_index(raw_path(cfg, cfg$data$worldbank_file)))

    } else if (s$short_name == "rvol") {
      rv <- construct_realized_volatility(read_yahoo_daily(raw_path(cfg, cfg$data$yahoo_file)),
                                          as.Date(cfg$data$yahoo_complete_through))
      x <- rv[, c("quarter", "value")]

    } else {
      stop("no construction rule for ", s$short_name, " (source ", s$source, ")")
    }

    x <- x[!is.na(x$value), , drop = FALSE]
    rownames(x) <- NULL
    out[[s$short_name]] <- x
  }

  # Backfills (docs/data.md 3.6): extend a series before its first observation.
  audit <- data.frame(series = character(), quarter = as.Date(character()), value = numeric(),
                      source = character())
  for (nm in names(cfg$data$backfill)) {
    b <- cfg$data$backfill[[nm]]
    if (b$source != "gsw") stop("unknown backfill source for ", nm, ": ", b$source)
    filler <- daily_to_quarterly(read_gsw_daily(raw_path(cfg, b$file), b$column),
                                 as.Date(b$complete_through))
    # Fill only inside the sample, so the audit lists exactly the quarters used.
    filler <- filler[filler$quarter >= parse_quarter(cfg$sample$start), , drop = FALSE]
    filled <- backfill_before_start(out[[nm]], filler)
    out[[nm]] <- filled$series
    if (nrow(filled$added) > 0) {
      audit <- rbind(audit, data.frame(series = nm, quarter = filled$added$quarter,
                                       value = filled$added$value,
                                       source = paste0("GSW ", b$column)))
    }
  }
  attr(out, "backfill") <- audit
  out
}

#' Extend a quarterly series back in time with a substitute series.
#'
#' Only quarters before the series' first observation are filled; where the
#' series exists, its own values are kept unchanged. The substitute is used
#' as is (no rescaling), so it must measure the same quantity in the same
#' units (docs/data.md 3.6 gives the evidence for GS2).
#'
#' @param x data.frame (quarter, value), the series to extend.
#' @param filler data.frame (quarter, value), the substitute.
#' @return list(series = extended data.frame, added = the rows taken from filler).
backfill_before_start <- function(x, filler) {
  added <- filler[filler$quarter < min(x$quarter), , drop = FALSE]
  series <- rbind(added, x)
  series <- series[order(series$quarter), , drop = FALSE]
  rownames(series) <- NULL
  list(series = series, added = added)
}

#' Line the series up on one quarterly index: the levels table.
#'
#' The sample starts at `start`. With end = "latest", it ends at the last
#' quarter in which all series have a value (docs/data.md, section 1), so no
#' series has a missing last quarter. Missing values inside the sample (for
#' example a series that starts after the sample start and has no backfill) are
#' kept as NA and reported by the validation checks; nothing is filled in.
#'
#' @param series_list Output of build_series_list().
#' @param start First quarter (Date or "1973Q1").
#' @param end Last quarter (Date or label), or "latest".
#' @return data.frame with column quarter and one column per short name, in
#'   config order, natural units.
assemble_levels <- function(series_list, start, end = "latest") {
  start <- if (inherits(start, "Date")) start else parse_quarter(start)
  last_obs  <- vapply(series_list, function(x) as.numeric(max(x$quarter)), numeric(1))
  if (identical(end, "latest")) {
    end <- as.Date(min(last_obs), origin = "1970-01-01")
  } else {
    end <- if (inherits(end, "Date")) end else parse_quarter(end)
  }
  quarters <- seq(start, end, by = "3 months")

  levels <- data.frame(quarter = quarters)
  for (nm in names(series_list)) {
    x <- series_list[[nm]]
    levels[[nm]] <- x$value[match(quarters, x$quarter)]
  }
  # With end = "latest", every series must have its last quarter.
  if (any(is.na(levels[nrow(levels), -1]))) {
    stop("some series have no value in the last quarter, ", format_quarter(end))
  }
  levels
}

#' Apply the transformations of docs/data.md: the model_data table.
#'
#' "log" series become 100 * ln(x), a level measured in log points, so that a
#' change of 1 is roughly a 1 percent change; "raw" series are unchanged
#' (rates and ratios in their natural units). No differencing: the VAR is in
#' levels (Crump et al. 2021, Section 2.1).
#'
#' @param levels Output of assemble_levels().
#' @param series Series config (gives the transform of each column).
#' @return data.frame with the same shape as levels. Stops if a log series has
#'   a zero or negative value (log undefined).
transform_to_model_units <- function(levels, series) {
  out <- levels
  for (i in seq_len(nrow(series))) {
    nm <- series$short_name[i]
    if (series$transform[i] == "log") {
      bad <- !is.na(levels[[nm]]) & levels[[nm]] <= 0
      if (any(bad)) {
        stop("non-positive values in log series ", nm, " at ",
             paste(format_quarter(levels$quarter[bad]), collapse = ", "))
      }
      out[[nm]] <- 100 * log(levels[[nm]])
    }
  }
  out
}


# ---- Constructed series (docs/data.md, section 3) ----
#
# Series that are not read directly from a source.
#
# Why it exists: one of the 31 model variables is computed by the project:
# row 31, S&P 500 realized volatility, from daily closes (3.4). Row 30 (World
# Bank non-energy index) and the monthly API series only need monthly-to-
# quarterly averaging, done with monthly_to_quarterly() (R/data_load.R).
#
# Inputs:  tables returned by the loaders in data_load.R.
# Outputs: quarterly data.frames with columns quarter and value, in natural
#          units (no logs; transformations come later, docs/data.md section 1).

#' Quarterly realized volatility of the S&P 500 (docs/data.md 3.4).
#'
#' Notation: P_d is the close on trading day d; r_d = ln(P_d / P_prev) is the
#' daily log return, using the previous available trading day's close. For
#' quarter q,
#'
#'   rvol_q = 100 * sqrt(252 * mean(r_d^2)),   mean over trading days d in q,
#'
#' the root mean square of daily returns (zero mean assumed, no demeaning),
#' annualized with 252 trading days and expressed in percent, so it is in the
#' same units as the VIX. The mean (not the sum) makes the measure independent
#' of the number of trading days in the quarter.
#'
#' Whole-quarter rule: a quarter is kept only if the daily data cover all of
#' it. The first quarter is kept only if the data start before it (so its
#' first return has a previous close); the last only if it has ended by
#' complete_through.
#'
#' @param daily data.frame with columns date (Date) and close (numeric, > 0),
#'   from read_yahoo_daily().
#' @param complete_through Date up to which the daily data are complete,
#'   normally the day before the download date.
#' @return data.frame with columns quarter (Date), value (rvol, annualized %)
#'   and n_days (trading days used).
construct_realized_volatility <- function(daily, complete_through) {
  daily <- daily[order(daily$date), , drop = FALSE]
  if (any(is.na(daily$close)) || any(daily$close <= 0)) {
    stop("daily closes must be positive and non-missing")
  }
  r <- data.frame(date = daily$date[-1], r = diff(log(daily$close)))
  r$quarter <- quarter_start(r$date)

  mean_sq <- tapply(r$r^2, r$quarter, mean)
  n_days  <- tapply(r$r, r$quarter, length)
  out <- data.frame(
    quarter = as.Date(names(mean_sq)),
    value   = 100 * sqrt(252 * as.numeric(mean_sq)),
    n_days  = as.integer(n_days)
  )

  first_close <- min(daily$date)
  keep <- out$quarter > first_close & quarter_end(out$quarter) <= as.Date(complete_through)
  out <- out[keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}
