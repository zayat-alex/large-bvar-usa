# data_pull.R: download raw data from external sources and save them as
# immutable files under data/raw/.
#
# Why it exists: the pipeline itself never goes online.
# Downloads happen once, through scripts/pull_raw.R, and are saved as dated raw
# snapshots; every later step reads those files only. This keeps every result
# reproducible from the files on disk.
#
# Inputs:  FRED API (key from the FRED_API_KEY environment variable; never printed).
# Outputs: CSV files written by save_raw_csv(), which refuses to overwrite.

#' Download the full vintage history of one FRED series from ALFRED.
#'
#' ALFRED keeps every published version ("vintage") of a FRED series. Asking for
#' a range of real-time dates returns one row per observation and per period in
#' which that value was the published one: a value is valid from realtime_start
#' to realtime_end. The series as it stood on a date d is the set of rows with
#' realtime_start <= d <= realtime_end.
#'
#' The series is downloaded at its native frequency, with no aggregation by
#' FRED: the API only aggregates (e.g. monthly to quarterly) for a single
#' real-time date, not over a range of vintages. Monthly series are averaged
#' to quarterly later, in our own code (docs/data.md, section 1).
#'
#' @param series_id FRED series ID, e.g. "GDPC1".
#' @param realtime_start First real-time date (Date).
#' @param realtime_end Last real-time date (Date), normally the retrieval date.
#' @param page_size Rows per API request; the request is repeated with an
#'   offset until all rows are returned.
#' @return A data.frame with columns series_id, date (first day of the month
#'   or quarter the value refers to), value, realtime_start, realtime_end.
#'   value is NA where FRED reports a missing value (".").
pull_alfred_vintages <- function(series_id, realtime_start, realtime_end,
                                 page_size = 100000) {
  if (!nzchar(Sys.getenv("FRED_API_KEY"))) {
    stop("FRED_API_KEY is not set; add it to the user-level ~/.Renviron")
  }

  pages  <- list()
  offset <- 0
  repeat {
    page <- fredr::fredr_series_observations(
      series_id      = series_id,
      realtime_start = as.Date(realtime_start),
      realtime_end   = as.Date(realtime_end),
      limit          = page_size,
      offset         = offset
    )
    pages[[length(pages) + 1]] <- page
    if (nrow(page) < page_size) break
    offset <- offset + page_size
  }

  out <- do.call(rbind, pages)
  out <- as.data.frame(out[, c("series_id", "date", "value", "realtime_start", "realtime_end")])
  out[order(out$date, out$realtime_start), , drop = FALSE]
}

#' Download current FRED metadata for a set of series (docs/data.md, check 4.7).
#'
#' @param series_ids Character vector of FRED series IDs.
#' @return A data.frame with one row per series: id, title, units,
#'   frequency, seasonal_adjustment, observation_start, observation_end,
#'   last_updated.
pull_fred_metadata <- function(series_ids) {
  rows <- lapply(series_ids, function(id) {
    m <- fredr::fredr_series(id)
    data.frame(
      id                  = m$id,
      title               = m$title,
      units               = m$units,
      frequency           = m$frequency,
      seasonal_adjustment = m$seasonal_adjustment,
      observation_start   = as.character(m$observation_start),
      observation_end     = as.character(m$observation_end),
      last_updated        = m$last_updated,
      stringsAsFactors    = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Download daily price history of one Yahoo Finance symbol.
#'
#' Used for row 31 (realized volatility, docs/data.md 3.4), computed later from
#' daily closes of the S&P 500 price index (^GSPC, not total return). The raw
#' file keeps every column Yahoo returns; nothing is filtered or filled here.
#'
#' @param symbol Yahoo symbol, e.g. "^GSPC".
#' @param from First date requested (Date). Yahoo returns what it has from
#'   that date on; the actual first date must be checked on load.
#' @param to Last date requested (Date).
#' @return A data.frame with columns date, open, high, low, close, volume,
#'   adjusted, one row per trading day.
pull_yahoo_daily <- function(symbol, from, to) {
  x <- quantmod::getSymbols(symbol, src = "yahoo", from = as.Date(from),
                            to = as.Date(to), auto.assign = FALSE)
  out <- data.frame(date = zoo::index(x), zoo::coredata(x), row.names = NULL)
  # quantmod names columns like "GSPC.Close"; keep only the field name.
  names(out)[-1] <- tolower(sub("^.*\\.", "", names(out)[-1]))
  out
}

#' Write a raw snapshot to CSV, refusing to overwrite an existing file.
#'
#' Raw files are immutable: a new download is always
#' a new file. Writing to an existing path is an error, not a silent overwrite.
#'
#' @param df Data frame to write.
#' @param path Destination path.
#' @return path, invisibly.
save_raw_csv <- function(df, path) {
  if (file.exists(path)) {
    stop("refusing to overwrite raw file: ", path)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, path, row.names = FALSE, na = "")
  invisible(path)
}
