# check_covid_comparison.R: the COVID treatments compared (forecasts, impulse
# responses, hyperparameters, and the diagnostics of R/forecast_drift_checks.R).
#
# Run from the project root (so .Rprofile activates renv), after
# scripts/check_forecast_review.R run and scripts/check_irf_review.R run
# (the report reads their results):
#   Rscript scripts/check_covid_comparison.R run
#   Rscript scripts/check_covid_comparison.R report
#
# For each production run (one per COVID method): the hyperparameter
# posterior and MCMC headline numbers (from the run's files), the drift of
# the no-shock forecast to 80 quarters, forecasts from earlier origins with
# the same estimates against what happened, and the effect of replacing
# unusual recent data (consumer sentiment, the oil price). No-shock paths
# only: nothing is simulated or re-estimated.
#
# Writes outputs/checks_covid_comparison/ (covid_comparison.rds,
# covid_comparison.html).

for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)

args  <- commandArgs(trailingOnly = TRUE)
stage <- args[1]
if (is.na(stage)) stop("give a stage: run or report")

cfg     <- yaml::read_yaml("config/run.yml")
series  <- read_series_config("config/series.csv")
vars    <- series$short_name
is_log  <- setNames(series$transform == "log", vars)
out_dir <- "outputs/checks_covid_comparison"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

show_vars    <- c("rgdp", "payems", "urate", "pce", "corepce", "gs2", "gs10")
origins      <- c("2017Q4", "2018Q4", "2022Q4", "2023Q4", "2024Q2", "2024Q4")
method_order <- c("none", "precovid_params", "exclude_window")

# ---- Stage: run ----
if (stage == "run") {
  runs <- list.files(cfg$output_root, pattern = "^[0-9]{4}m[0-9]{2}_.+_full_[0-9a-f]{8}$", full.names = TRUE)
  methods <- vapply(runs, function(d) jsonlite::read_json(file.path(d, "meta.json"))$covid$method, character(1))
  res <- list(time = Sys.time(), methods = list())
  for (k in seq_along(runs)) {
    m <- methods[[k]]; d <- runs[k]
    t0 <- proc.time()[["elapsed"]]
    dr <- readRDS(file.path(d, "draws.rds"))
    cv <- readRDS(file.path(d, "covid.rds"))
    hist <- est_data_to_matrix(cv$origin_data, series)
    summ <- readRDS(file.path(d, "mcmc_summary.rds"))
    est <- readRDS(file.path(d, "estimation.rds"))
    n <- nrow(hist)
    pre <- parse_quarter(rownames(hist)) <= parse_quarter("2019Q4")

    med80 <- no_shock_median(dr, hist, 80, is_log)
    cf <- rbind(
      cbind(change = "umich at its 1973-2019 mean (last 5 quarters)",
            counterfactual_effect(dr, hist, is_log, "umich", mean(hist[pre, "umich"]), 5, c(4, 8, 12), show_vars)),
      cbind(change = "oil flat at its 2025Q1 level (last 5 quarters)",
            counterfactual_effect(dr, hist, is_log, "wti", hist[n - 5, "wti"], 5, c(4, 8, 12), show_vars)))

    res$methods[[m]] <- list(
      run_id = basename(d),
      est_quarters = if (is.null(cv$est_weights)) nrow(cv$est_data) - dr$lags else sum(cv$est_weights[-seq_len(dr$lags)]),
      est_sample = paste(format_quarter(min(cv$est_data$quarter)), "to", format_quarter(max(cv$est_data$quarter))),
      excluded = cv$notes$quarter_label[grepl("excluded|dropped", cv$notes$action)],
      mode = exp(est$mode$theta_mode),
      hyper = summ$hyper, convergence = summ$convergence, rhat_multivariate = summ$rhat_multivariate,
      max_modulus = summ$max_modulus_summary,
      drift = drift_table(med80, show_vars),
      no_shock = med80[seq_len(40), show_vars],
      origin_errors = origin_errors(dr, hist, is_log, origins, c(4, 8), show_vars),
      counterfactual = cf,
      seconds = proc.time()[["elapsed"]] - t0)
    cat(m, "done in", round(res$methods[[m]]$seconds), "s\n")
  }
  res$history_change <- history_change_per_year(est_data_to_matrix(cv$origin_data, series), show_vars, is_log)
  res$methods <- res$methods[method_order[method_order %in% names(res$methods)]]
  saveRDS(res, file.path(out_dir, "covid_comparison.rds"))
}

# ---- Stage: report ----
if (stage == "report") {
  out_abs <- normalizePath(out_dir, winslash = "/")
  quarto::quarto_render("reports/covid_comparison.qmd",
                        execute_params = list(output_dir = out_abs,
                                              forecast_dir = normalizePath("outputs/checks_forecast_review", winslash = "/"),
                                              irf_dir = normalizePath("outputs/checks_irf_review", winslash = "/")),
                        quiet = TRUE)
  target <- file.path(out_dir, "covid_comparison.html")
  if (file.exists(target)) file.remove(target)
  if (!file.rename("reports/covid_comparison.html", target)) stop("could not move the report to ", target)
  cat("Wrote", target, "\n")
}
