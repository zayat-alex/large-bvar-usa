# pull_raw.R: download raw data into data/raw/ (outside the targets pipeline).
#
# Run from the project root (so .Rprofile activates renv):
#   Rscript scripts/pull_raw.R            # all sources
#   Rscript scripts/pull_raw.R yahoo      # only the named source(s): alfred, yahoo
#   Rscript scripts/pull_raw.R alfred GS2 # only the named ALFRED series (and their metadata)
#
# Reads config/raw_pulls.yml. Every file is new and dated; existing files are
# never overwritten. After a run, record the files in data/raw/SOURCE.md.
#
# Sources: alfred (ALFRED vintage histories plus FRED metadata) and
# yahoo (^GSPC daily closes). An output already present for today's date is
# skipped, so a run that stopped partway can be restarted safely.

source("R/data_pull.R")

all_sources <- c("alfred", "yahoo")
args        <- commandArgs(trailingOnly = TRUE)
sources     <- intersect(args, all_sources)
series_only <- setdiff(args, all_sources)       # ALFRED series IDs to restrict to
if (length(sources) == 0) sources <- all_sources
cfg       <- yaml::read_yaml("config/raw_pulls.yml")
if (length(series_only) > 0) {
  known <- vapply(cfg$alfred$series, function(s) s$id, character(1))
  if (!all(series_only %in% known)) {
    stop("not a source or a series in config/raw_pulls.yml: ",
         paste(setdiff(series_only, known), collapse = ", "))
  }
  cfg$alfred$series <- Filter(function(s) s$id %in% series_only, cfg$alfred$series)
}
# FRED's "today" is the US Central date; a later local date (e.g. after midnight
# in South America) is rejected as a realtime_end, so use FRED's date.
retrieved <- as.Date(format(Sys.time(), tz = "America/Chicago"))
cat("Retrieval date:", format(retrieved), " sources:", paste(sources, collapse = ", "), "\n\n")

# ---- ALFRED vintage histories, one file per series, plus metadata ----
if ("alfred" %in% sources) {
  alfred   <- cfg$alfred
  rt_start <- as.Date(alfred$realtime_start)
  cat("FRED_API_KEY set:", nzchar(Sys.getenv("FRED_API_KEY")), "\n")

  for (s in alfred$series) {
    path <- file.path(alfred$out_dir, sprintf("%s_vintages_%s_to_%s.csv",
                                              s$id, format(rt_start), format(retrieved)))
    if (file.exists(path)) {
      cat(sprintf("%-11s already pulled, skipped: %s\n", s$id, path))
      next
    }
    df <- pull_alfred_vintages(s$id, rt_start, retrieved)
    save_raw_csv(df, path)
    cat(sprintf("%-11s rows=%6d  obs %s to %s  real-time periods=%d  -> %s\n",
                s$id, nrow(df), format(min(df$date)), format(max(df$date)),
                length(unique(df$realtime_start)), path))
  }

  ids       <- vapply(alfred$series, function(s) s$id, character(1))
  meta_name <- if (length(series_only) > 0) {
    sprintf("fred_metadata_%s_%s.csv", paste(ids, collapse = "_"), format(retrieved))
  } else {
    sprintf("fred_metadata_%s.csv", format(retrieved))
  }
  meta_path <- file.path(alfred$out_dir, meta_name)
  if (file.exists(meta_path)) {
    cat("\nMetadata already pulled, skipped:", meta_path, "\n")
  } else {
    meta <- pull_fred_metadata(ids)
    meta$retrieved <- format(retrieved)
    save_raw_csv(meta, meta_path)
    cat("\nMetadata ->", meta_path, "\n")
    print(meta[, c("id", "frequency", "units", "seasonal_adjustment", "last_updated")],
          row.names = FALSE)
  }
}

# ---- Yahoo daily closes ----
if ("yahoo" %in% sources) {
  yahoo      <- cfg$yahoo
  yahoo_path <- file.path(yahoo$out_dir, sprintf("%s_daily_%s.csv",
                                                 gsub("\\^", "", yahoo$symbol), format(retrieved)))
  if (file.exists(yahoo_path)) {
    cat("\nYahoo already pulled, skipped:", yahoo_path, "\n")
  } else {
    daily <- pull_yahoo_daily(yahoo$symbol, yahoo$from, retrieved)
    save_raw_csv(daily, yahoo_path)
    cat(sprintf("\n%s daily: rows=%d  %s to %s  missing closes=%d  -> %s\n",
                yahoo$symbol, nrow(daily), format(min(daily$date)), format(max(daily$date)),
                sum(is.na(daily$close)), yahoo_path))
  }
}
