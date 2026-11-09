# render_reports_md.R: the check reports as Markdown pages that GitHub
# displays, written to reports/results/.
#
# Why it exists: each check report (reports/*.qmd) is rendered to a
# self-contained HTML file in its output folder (outputs/..., not part of
# the repository), and GitHub shows an HTML file as source code rather than
# as a page. This script renders the same reports again as GitHub Markdown,
# from the results their scripts saved (nothing is recomputed):
#   reports/results/<report>.md          the report
#   reports/results/<report>_figures/    its figures
#
# Run from the project root (so .Rprofile activates renv), after the
# reports' own scripts have saved their results:
#   Rscript scripts/render_reports_md.R                  all reports
#   Rscript scripts/render_reports_md.R forecast_review  only the named ones

cfg <- yaml::read_yaml("config/run.yml")
results <- "reports/results"
dir.create(results, showWarnings = FALSE, recursive = TRUE)
abs_path <- function(p) normalizePath(p, winslash = "/")

# The reports, with the parameters their scripts pass to quarto (see each
# script's `report` stage, or the pipeline's render functions for the data
# validation and MCMC reports).
check <- function(name, dir) list(qmd = name, out = name, params = list(output_dir = abs_path(dir)))
reports <- list(
  check("data_validation", cfg$data$output_dir),
  check("conditional_checks", "outputs/checks_conditional"),
  check("forecast_checks", "outputs/checks_forecast"),
  check("forecast_review", "outputs/checks_forecast_review"),
  check("irf_review", "outputs/checks_irf_review"),
  check("estimation_diagnostics", "outputs/checks_estimation_diagnostics"),
  check("hyper_stability", "outputs/checks_hyper_stability"),
  check("prior_checks", "outputs/checks_priors"),
  check("prior_sensitivity", "outputs/checks_prior_sensitivity"),
  check("recovery_checks", "outputs/checks_recovery"),
  check("tilting_checks", "outputs/checks_tilting"),
  list(qmd = "covid_comparison", out = "covid_comparison",
       params = list(output_dir = abs_path("outputs/checks_covid_comparison"),
                     forecast_dir = abs_path("outputs/checks_forecast_review"),
                     irf_dir = abs_path("outputs/checks_irf_review"))))
# The MCMC report of each estimation run (one per COVID method).
for (run_dir in dirname(Sys.glob(file.path(cfg$output_root, "*_full_*", "mcmc_diagnostics.html")))) {
  method <- sub("^[0-9]{4}m[0-9]{2}_(.*)_full_[0-9a-f]+$", "\\1", basename(run_dir))
  reports[[length(reports) + 1]] <- list(qmd = "mcmc_diagnostics", out = paste0("mcmc_diagnostics_", method),
                                         params = list(run_dir = abs_path(run_dir)))
}

#' Render one report to GitHub Markdown and gather its figures.
#'
#' The report is rendered in reports/ (its code may read files relative to
#' it); its images, wherever they are (the output folder, or the figures
#' Quarto writes beside the page), are copied to <out>_figures/ and the
#' links rewritten to point there.
#'
#' @param r One entry of `reports`.
#' @return Number of figures, invisibly.
render_md <- function(r) {
  quarto::quarto_render(file.path("reports", paste0(r$qmd, ".qmd")), output_format = "gfm",
                        execute_params = r$params, pandoc_args = "--wrap=none", quiet = TRUE)
  md_file <- file.path("reports", paste0(r$qmd, ".md"))
  quarto_figs <- file.path("reports", paste0(r$qmd, "_files"))
  md <- readLines(md_file, warn = FALSE, encoding = "UTF-8")
  fig_dir <- file.path(results, paste0(r$out, "_figures"))
  unlink(fig_dir, recursive = TRUE)
  dir.create(fig_dir, showWarnings = FALSE)
  # Image files linked from the page: ](path) and src="path".
  pattern <- "(\\]\\(|src=\")([^)\"]+[.](png|jpe?g|svg|gif))"
  text <- paste(md, collapse = "\n")
  paths <- unique(sub(pattern, "\\2", regmatches(text, gregexpr(pattern, text, perl = TRUE))[[1]], perl = TRUE))
  paths <- paths[!grepl("^(https?:|data:)", paths)]
  for (k in seq_along(paths)) {
    p <- utils::URLdecode(paths[k])                  # Quarto encodes spaces as %20
    src <- if (file.exists(p)) p else file.path("reports", p)
    if (!file.exists(src)) stop(r$out, ": figure not found: ", paths[k])
    name <- sprintf("%02d_%s", k, basename(src))
    file.copy(src, file.path(fig_dir, name), overwrite = TRUE)
    md <- gsub(paths[k], paste0(basename(fig_dir), "/", name), md, fixed = TRUE)
  }
  writeLines(md, file.path(results, paste0(r$out, ".md")), useBytes = TRUE)
  file.remove(md_file)
  unlink(quarto_figs, recursive = TRUE)
  cat(sprintf("%-36s %3d figures\n", r$out, length(paths)))
  invisible(length(paths))
}

only <- commandArgs(trailingOnly = TRUE)
if (length(only)) reports <- Filter(function(r) r$out %in% only, reports)
for (r in reports) render_md(r)
