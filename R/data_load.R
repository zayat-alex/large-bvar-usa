# data_load.R: read the series configuration and the raw files under
# data/raw/ into clean R tables.
#
# Why it exists: every raw source has its own layout (QD's two header rows and
# end-of-quarter date stamps, ALFRED's real-time periods, Yahoo's daily rows,
# the World Bank's Excel sheet). These functions turn each one into a plain
# table indexed by date, so the construction and assembly steps do not need to
# know anything about file formats.
#
# Conventions (docs/data.md, section 1):
#   - A quarter is stored as the R Date of its first day (1973-01-01 = 1973:Q1)
#     and printed as "1973Q1".
#   - Monthly data are averaged to quarterly only when all three months of the
#     quarter are available (the "whole-quarter rule").
#   - Nothing here transforms values (no logs, no QD transformation codes).
#
# Inputs:  config/series.csv; files under data/raw/ (see data/raw/SOURCE.md).
# Outputs: data.frames returned to the caller; nothing is written to disk.


# ---- The series configuration (config/series.csv) ----
#
# Reading and checking the table of the 31 model variables.
#
# Why it exists: every later data step (loading QD columns, pulling API series,
# building the constructed series, applying 100*ln or no transform, setting the
# Minnesota prior mean) is driven by one table with one row per model variable.
# Keeping that table in a file, instead of hard-coding names in functions, means
# the dataset definition lives in one place and can be checked against the
# series table in docs/data.md (section 2).
#
# Input:  config/series.csv (31 rows, one per variable, in the order of Table 1
#         of Crump et al. 2021).
# Output: a data.frame with the columns described in read_series_config().

# Allowed values of the coded columns (docs/data.md, section 2).
#   source:    QD = column of the FRED-QD CSV; API = FRED/ALFRED API;
#              CALC = constructed from other data; MANUAL = file downloaded by hand.
#   transform: log = 100 * ln(x), a level in log points (not a growth rate);
#              raw = used as published (rates, ratios, indices in natural units).
#   prior:     prior mean of the variable's own first-lag coefficient in the
#              Minnesota prior. RW = random walk (mean 1), WN = white noise (mean 0).
#   api_aggregation: only for API rows. none = series is already quarterly;
#              avg = monthly series averaged to quarterly by FRED.
SERIES_SOURCES      <- c("QD", "API", "CALC", "MANUAL")
SERIES_TRANSFORMS   <- c("log", "raw")
SERIES_PRIORS       <- c("RW", "WN")
SERIES_AGGREGATIONS <- c("none", "avg")


#' Read the series configuration and check its structure.
#'
#' @param path Path to config/series.csv.
#' @return A data.frame with one row per model variable and columns
#'   row (integer, 1 to 31, order of Table 1 in the paper), short_name,
#'   variable, source, id, transform, prior, api_aggregation (NA unless
#'   source is API) and notes.
#'
#' The function stops with an error if the file breaks a structural rule
#' (unknown codes, duplicated names, rows out of order, API rows without an
#' aggregation rule). It does not check the content against the docs/data.md;
#' that comparison is a test (tests/testthat/test-data-config.R).
read_series_config <- function(path) {
  series <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names      = FALSE,
    na.strings       = "",
    strip.white      = TRUE
  )

  expected_cols <- c("row", "short_name", "variable", "source", "id",
                     "transform", "prior", "api_aggregation", "notes")
  missing_cols <- setdiff(expected_cols, names(series))
  if (length(missing_cols) > 0) {
    stop("series config is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  # Rows must be numbered 1..n in order, because the order of the variables is
  # the order of the VAR and later the Cholesky ordering for the IRFs.
  if (!identical(as.integer(series$row), seq_len(nrow(series)))) {
    stop("series config: column 'row' must run 1, 2, ..., n in order")
  }
  series$row <- as.integer(series$row)

  if (anyDuplicated(series$short_name) > 0) {
    stop("series config: duplicated short_name: ",
         paste(series$short_name[duplicated(series$short_name)], collapse = ", "))
  }
  if (any(is.na(series$id))) {
    stop("series config: every row needs an id")
  }

  check_codes <- function(column, allowed) {
    bad <- setdiff(unique(series[[column]]), allowed)
    if (length(bad) > 0) {
      stop("series config: unknown value(s) in '", column, "': ", paste(bad, collapse = ", "))
    }
  }
  check_codes("source", SERIES_SOURCES)
  check_codes("transform", SERIES_TRANSFORMS)
  check_codes("prior", SERIES_PRIORS)

  # API rows need an aggregation rule; other rows must not have one, so that a
  # stray value cannot be silently ignored.
  is_api <- series$source == "API"
  if (any(is.na(series$api_aggregation[is_api])) ||
      any(!series$api_aggregation[is_api] %in% SERIES_AGGREGATIONS)) {
    stop("series config: API rows need api_aggregation in {",
         paste(SERIES_AGGREGATIONS, collapse = ", "), "}")
  }
  if (any(!is.na(series$api_aggregation[!is_api]))) {
    stop("series config: api_aggregation is only allowed on API rows")
  }

  series
}


# ---- Quarter helpers ----

#' First day of the quarter that contains each date.
#'
#' @param d Date vector.
#' @return Date vector: 1 January, 1 April, 1 July or 1 October.
quarter_start <- function(d) {
  d     <- as.Date(d)
  year  <- as.integer(format(d, "%Y"))
  month <- as.integer(format(d, "%m"))
  as.Date(sprintf("%04d-%02d-01", year, (month - 1) %/% 3 * 3 + 1))
}

#' Last calendar day of the quarter that contains each date.
#'
#' @param d Date vector.
#' @return Date vector: 31 March, 30 June, 30 September or 31 December.
quarter_end <- function(d) {
  # The day before the start of the next quarter. Adding 92 days to the first
  # day of any quarter lands inside the next quarter (quarters have 90 to 92 days).
  quarter_start(quarter_start(d) + 92) - 1
}

#' Print quarters as "1973Q1".
#'
#' @param q Date vector (first days of quarters).
#' @return Character vector of the same length as q: an empty vector for no
#'   quarters, and NA for a missing date. (paste0() alone would return "Q" for
#'   an empty input and "NAQNA" for a missing date, so both cases are handled
#'   explicitly.)
format_quarter <- function(q) {
  q   <- as.Date(q)
  out <- rep(NA_character_, length(q))
  ok  <- !is.na(q)
  out[ok] <- paste0(format(q[ok], "%Y"), "Q", (as.integer(format(q[ok], "%m")) - 1) %/% 3 + 1)
  out
}

#' Parse "1973Q1" (or "1973:Q1") into the first day of the quarter.
#'
#' @param x Character vector.
#' @return Date vector. Stops on anything that is not a valid quarter label.
parse_quarter <- function(x) {
  x  <- gsub(":", "", as.character(x))
  ok <- grepl("^[0-9]{4}Q[1-4]$", x)
  if (!all(ok)) stop("not a quarter label (expected e.g. 1973Q1): ", paste(x[!ok], collapse = ", "))
  year <- as.integer(substr(x, 1, 4))
  qtr  <- as.integer(substr(x, 6, 6))
  as.Date(sprintf("%04d-%02d-01", year, (qtr - 1) * 3 + 1))
}

#' Convert QD date stamps to quarters.
#'
#' FRED-QD stamps each quarter with the first day of its LAST month
#' ("3/1/1973" = 1973:Q1). The project stores the first day of the quarter.
#'
#' @param x Character vector in QD's m/d/yyyy format.
#' @return Date vector (first days of quarters). Stops if a stamp is not the
#'   first day of March, June, September or December, which would mean the
#'   file does not have the expected layout.
qd_date_to_quarter <- function(x) {
  d <- as.Date(x, format = "%m/%d/%Y")
  if (any(is.na(d))) stop("unreadable QD date stamp(s): ", paste(x[is.na(d)], collapse = ", "))
  month <- as.integer(format(d, "%m"))
  if (any(!month %in% c(3, 6, 9, 12)) || any(format(d, "%d") != "01")) {
    stop("QD date stamps must be the first day of the last month of a quarter")
  }
  quarter_start(d)
}


# ---- Monthly to quarterly ----

#' Average a monthly series to quarterly, applying the whole-quarter rule.
#'
#' A quarter gets a value only if all three of its months have a non-missing
#' value; otherwise it is dropped (docs/data.md, section 1). This is how FRED's
#' own quarterly "avg" aggregation and FRED-QD treat monthly series, and it
#' avoids a last quarter computed from one or two months.
#'
#' @param df data.frame with columns date (Date, any day within the month;
#'   normally the first) and value (numeric).
#' @return data.frame with columns quarter (Date, first day) and value,
#'   sorted by quarter, complete quarters only.
monthly_to_quarterly <- function(df) {
  if (anyDuplicated(format(df$date, "%Y-%m")) > 0) stop("monthly series has duplicated months")
  df <- df[!is.na(df$value), , drop = FALSE]
  df$quarter <- quarter_start(df$date)
  n_months <- tapply(df$value, df$quarter, length)
  means    <- tapply(df$value, df$quarter, mean)
  out <- data.frame(quarter = as.Date(names(means)), value = as.numeric(means))
  out <- out[as.integer(n_months) == 3, , drop = FALSE]
  rownames(out) <- NULL
  out[order(out$quarter), , drop = FALSE]
}


# ---- FRED-QD vintage ----

#' Read selected columns from a FRED-QD vintage CSV.
#'
#' Layout of the file (data/raw/SOURCE.md): row 1 holds the column names
#' (sasdate, then the series); two rows labelled "factors" and "transform"
#' follow; then one row per quarter. The "transform" codes are QD's suggested
#' stationarity transformations, which the project does not use (the VAR is in
#' levels; docs/data.md, section 1), so both rows are dropped.
#'
#' @param path Path to a FRED-QD_<yyyy>m<mm>.csv file.
#' @param columns Character vector of column names to keep (exact match,
#'   e.g. "S&P 500").
#' @return data.frame with column quarter (Date) and one numeric column per
#'   requested name. Empty cells become NA.
read_qd_vintage <- function(path, columns) {
  raw <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                         colClasses = "character", na.strings = "")
  missing_cols <- setdiff(columns, names(raw))
  if (length(missing_cols) > 0) {
    stop("columns not in ", basename(path), ": ", paste(missing_cols, collapse = ", "))
  }
  # Drop the two metadata rows by their label, and any blank trailing rows.
  keep <- !is.na(raw$sasdate) & !raw$sasdate %in% c("factors", "transform")
  raw  <- raw[keep, , drop = FALSE]

  out <- data.frame(quarter = qd_date_to_quarter(raw$sasdate))
  for (col in columns) {
    out[[col]] <- as.numeric(raw[[col]])
  }
  if (anyDuplicated(out$quarter) > 0) stop("duplicated quarters in ", basename(path))
  out
}


# ---- ALFRED vintage histories ----

#' Read an ALFRED vintage-history file written by scripts/pull_raw.R.
#'
#' @param path Path to <ID>_vintages_<start>_to_<end>.csv.
#' @return data.frame with columns series_id, date, value, realtime_start,
#'   realtime_end (dates as Date, value numeric, NA where nothing was published).
read_alfred_vintages <- function(path) {
  v <- utils::read.csv(path, stringsAsFactors = FALSE, na.strings = "")
  for (col in c("date", "realtime_start", "realtime_end")) v[[col]] <- as.Date(v[[col]])
  v$value <- as.numeric(v$value)
  v
}

#' The series as it was published on a given date.
#'
#' Each row of an ALFRED vintage history holds a value together with the
#' period over which it was the published figure (realtime_start to
#' realtime_end, inclusive). The series as of date d is the set of rows whose
#' period contains d: exactly one row per observation date.
#'
#' @param vintages data.frame from read_alfred_vintages().
#' @param as_of Date of the vintage wanted (config: data$api_realtime_date).
#' @return data.frame with columns date and value, sorted by date, rows with a
#'   missing value dropped. Stops if as_of is outside the downloaded real-time
#'   window, or if any observation has more than one value on that date.
alfred_as_of <- function(vintages, as_of) {
  as_of <- as.Date(as_of)
  if (as_of < min(vintages$realtime_start) || as_of > max(vintages$realtime_end)) {
    stop("as_of date ", as_of, " is outside the downloaded real-time window (",
         min(vintages$realtime_start), " to ", max(vintages$realtime_end), ")")
  }
  s <- vintages[vintages$realtime_start <= as_of & vintages$realtime_end >= as_of,
                c("date", "value"), drop = FALSE]
  if (anyDuplicated(s$date) > 0) stop("more than one value per date as of ", as_of)
  s <- s[!is.na(s$value), , drop = FALSE]
  rownames(s) <- NULL
  s[order(s$date), , drop = FALSE]
}


# ---- Yahoo daily closes ----

#' Read a Yahoo daily file written by scripts/pull_raw.R.
#'
#' @param path Path to GSPC_daily_<date>.csv.
#' @return data.frame with columns date (Date) and close (numeric), sorted by
#'   date. Only the close is used (docs/data.md 3.4: price index, daily closes).
read_yahoo_daily <- function(path) {
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  out <- data.frame(date = as.Date(d$date), close = as.numeric(d$close))
  if (anyDuplicated(out$date) > 0) stop("duplicated trading days in ", basename(path))
  out[order(out$date), , drop = FALSE]
}


# ---- Gürkaynak-Sack-Wright Treasury yield curve ----

#' Read one column of the Federal Reserve's GSW yield curve file.
#'
#' The file (feds200628.csv) starts with a few lines of notes; the data begin
#' at the line starting with "Date". Values are daily, in percent; "NA" marks
#' maturities not estimated on that day.
#'
#' @param path Path to feds200628_<date>.csv.
#' @param column Column name, e.g. "SVENPY02" (2-year par yield, coupon-equivalent).
#' @return data.frame with columns date (Date) and value (numeric), missing
#'   values dropped, sorted by date.
read_gsw_daily <- function(path, column) {
  head_lines <- readLines(path, n = 50)
  skip <- grep("^Date,", head_lines) - 1
  if (length(skip) != 1) stop("no 'Date' header line in ", basename(path))
  x <- utils::read.csv(path, skip = skip, na.strings = "NA", stringsAsFactors = FALSE)
  if (!column %in% names(x)) stop("column ", column, " not in ", basename(path))
  out <- data.frame(date = as.Date(x$Date), value = as.numeric(x[[column]]))
  out <- out[!is.na(out$value), , drop = FALSE]
  out[order(out$date), , drop = FALSE]
}

#' Average a daily series to quarterly, keeping complete quarters only.
#'
#' A quarter is kept only if the daily data cover all of it: the data start
#' before the quarter begins and run past its end (`complete_through`).
#'
#' @param df data.frame with columns date and value (trading days only).
#' @param complete_through Date up to which the daily data are complete.
#' @return data.frame with columns quarter and value.
daily_to_quarterly <- function(df, complete_through) {
  df$quarter <- quarter_start(df$date)
  means <- tapply(df$value, df$quarter, mean)
  out <- data.frame(quarter = as.Date(names(means)), value = as.numeric(means))
  keep <- out$quarter > min(df$date) & quarter_end(out$quarter) <= as.Date(complete_through)
  out <- out[keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}


# ---- World Bank Pink Sheet ----

#' Read the monthly non-energy price index from a World Bank Pink Sheet file.
#'
#' The header row sits at different positions in different releases
#' (docs/data.md 3.3), so the column is located by its header text, not by position.
#' Data rows are those whose first cell is a month label like "1960M01".
#'
#' @param path Path to CMO-Historical-Data-Monthly*.xlsx.
#' @param sheet Sheet name, "Monthly Indices".
#' @param header Exact header text of the column, "Non-energy **".
#' @return data.frame with columns date (Date, first day of the month) and
#'   value (index, nominal US dollars, 2010 = 100).
read_worldbank_index <- function(path, sheet = "Monthly Indices", header = "Non-energy **") {
  x <- readxl::read_excel(path, sheet = sheet, col_names = FALSE, col_types = "text",
                          .name_repair = "minimal")
  x <- as.data.frame(x, stringsAsFactors = FALSE)

  hits <- which(as.matrix(x) == header, arr.ind = TRUE)
  if (nrow(hits) != 1) {
    stop("expected exactly one cell '", header, "' in sheet '", sheet, "', found ", nrow(hits))
  }
  col <- hits[1, "col"]

  is_month <- grepl("^[0-9]{4}M[0-9]{2}$", x[[1]])
  labels   <- x[[1]][is_month]
  out <- data.frame(
    date  = as.Date(sprintf("%s-%s-01", substr(labels, 1, 4), substr(labels, 6, 7))),
    value = as.numeric(x[[col]][is_month])
  )
  if (anyDuplicated(out$date) > 0) stop("duplicated months in ", basename(path))
  out[order(out$date), , drop = FALSE]
}
