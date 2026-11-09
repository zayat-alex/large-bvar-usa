# data_outputs.R: write the data stage's outputs to disk (pipeline stage 1).
#
# Why it exists: the targets pipeline (_targets.R) keeps the computations
# (assembly, checks) separate from the files they produce, so each output is a
# tracked file target that is rebuilt only when its inputs change. These
# functions do the writing; each returns the paths it wrote, which is what
# targets needs for a file target.
#
# Inputs:  tables and check results from data_assemble.R and data_validate.R.
# Outputs: files under cfg$data$output_dir (e.g. outputs/data_2026m09/):
#   levels.csv, model_data.csv, data.rds, validation.rds, checks/*.png,
#   exploratory_series.pdf, data_validation.html.


#' Every raw file the run config names, as full paths.
#'
#' Used by the pipeline to track the raw inputs: if any file changes (or the
#' config points to a different one), everything built from it is rebuilt.
#'
#' @param cfg Parsed config/run.yml.
#' @return Character vector of paths under cfg$data$raw_dir. Stops if one is missing.
raw_file_paths <- function(cfg) {
  d <- cfg$data
  rel <- c(d$qd_vintage_file, d$worldbank_file, d$yahoo_file,
           unlist(d$alfred_files), unlist(d$alfred_metadata_files),
           d$qd_twexmmth_file, d$worldbank_previous_file,
           vapply(d$backfill, function(b) b$file, character(1)))
  vapply(unique(rel), function(p) raw_path(cfg, p), character(1), USE.NAMES = FALSE)
}

#' One-line description of the data behind a build (printed on every output).
sample_description <- function(levels, cfg) {
  sprintf("QD vintage %s, ALFRED as of %s, sample %s to %s (%d quarters)",
          basename(cfg$data$qd_vintage_file), cfg$data$api_realtime_date,
          format_quarter(min(levels$quarter)), format_quarter(max(levels$quarter)), nrow(levels))
}

#' Write the two data tables (CSV) and data.rds.
#'
#' data.rds holds both tables, the series config, the run config passed in
#' (the pipeline passes its data and sample sections), the backfill audit and
#' the R session information, so a build can be traced.
#'
#' @return Paths written.
write_data_tables <- function(levels, model_data, series, cfg, backfill, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  write_table <- function(df, path) {
    df$quarter <- format_quarter(df$quarter)
    utils::write.csv(df, path, row.names = FALSE, na = "")
    path
  }
  rds <- file.path(dir, "data.rds")
  saveRDS(list(levels = levels, model_data = model_data, series = series, config = cfg,
               backfill = backfill, sample = sample_description(levels, cfg),
               built = Sys.time(), session = utils::sessionInfo()), rds)
  c(write_table(levels, file.path(dir, "levels.csv")),
    write_table(model_data, file.path(dir, "model_data.csv")),
    rds)
}

#' Save the check results read by the validation report.
#'
#' @return Path written.
write_validation_results <- function(checks, levels, cfg, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, "validation.rds")
  saveRDS(list(checks = checks, sample = sample_description(levels, cfg), built = Sys.time()), path)
  path
}

#' Render the Quarto validation report into the output folder.
#'
#' The report only reads validation.rds and the check plots; it runs no checks.
#'
#' @param qmd Report source (reports/data_validation.qmd).
#' @param dir Output folder holding validation.rds and checks/.
#' @param depends_on Paths of the files the report reads (passed only so the
#'   pipeline rebuilds the report when they change).
#' @return Path of the HTML report.
render_validation_report <- function(qmd, dir, depends_on = NULL) {
  out_abs <- normalizePath(dir, winslash = "/")
  quarto::quarto_render(qmd, execute_params = list(output_dir = out_abs), quiet = TRUE)
  rendered <- sub("\\.qmd$", ".html", qmd)
  target   <- file.path(dir, "data_validation.html")
  if (file.exists(target)) file.remove(target)
  if (!file.rename(rendered, target)) stop("could not move the rendered report to ", target)
  target
}
