# check_forecast_review.R: forecast reasonableness at the latest origin
# (spec 9.4; R/forecast_review.R).
#
# Run from the project root (so .Rprofile activates renv):
#   Rscript scripts/check_forecast_review.R run
#   Rscript scripts/check_forecast_review.R report
#
# For each production run (one per COVID method), from its saved draws and
# engine outputs (nothing is re-estimated):
#   1. the baseline against the latest SPF median (data/raw/spf/), in the
#      SPF's measures: exact annualized quarterly growth, levels for rates;
#   2. history-based plausible ranges (5th-95th percentiles, 1973Q1-2019Q4);
#   3. the share of draws with negative yields, to 40 quarters;
#   4. the baseline simulated to 40 quarters (in this check only; the
#      engine's horizon is 20), against the history means and ranges;
#   5. the inflation-anchored scenario (entropic tilting) against the
#      untilted forecast;
#   6. COVID sensitivity: baselines and two scenarios' effects per method;
#   7. simulation noise: the baseline recomputed with a second seed;
#   8. a readable table of the baseline (both methods, 12 quarters):
#      baseline_table.csv.
#
# Seeds: the engine's (mcmc$seed + forecast$seed_offset) reproduces its
# baseline; the noise run uses mcmc$seed + 5000.
#
# Writes outputs/checks_forecast_review/ (forecast_review.rds, PNG figures,
# forecast_review.html).

for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)

args  <- commandArgs(trailingOnly = TRUE)
stage <- args[1]
if (is.na(stage)) stop("give a stage: run or report")

cfg     <- yaml::read_yaml("config/run.yml")
series  <- read_series_config("config/series.csv")
vars    <- series$short_name
is_log  <- setNames(series$transform == "log", vars)
lags    <- cfg$model$lags
out_dir <- "outputs/checks_forecast_review"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

engine_seed <- cfg$mcmc$seed + cfg$forecast$seed_offset
noise_seed  <- cfg$mcmc$seed + 5000
long_h      <- 40
main_vars   <- c("rgdp", "rpce", "payems", "urate", "pce", "corepce", "cpi", "gs2", "gs10", "houst")
yields      <- c("gs2", "gs5", "gs10")
spf_survey  <- list(year = 2026, quarter = 3,
                    level = "data/raw/spf/medianLevel_2026-11-01.xlsx",
                    growth = "data/raw/spf/medianGrowth_2026-11-01.xlsx")
# SPF sheet and column prefix for each of our variables (SPF documentation):
# growth sheets give annualized quarterly growth of the median level
# forecast; the inflation sheets give median annualized quarterly inflation;
# UNEMP, TBOND and BAABOND are levels.
spf_map <- data.frame(
  variable = c("rgdp", "rpce", "rnresinv", "rresinv", "rfedgov", "indpro", "payems", "houst",
               "cpi", "corecpi", "pce", "corepce", "urate", "gs10", "baa"),
  file = c(rep("growth", 8), rep("level", 7)),
  sheet = c("RGDP", "RCONSUM", "RNRESIN", "RRESINV", "RFEDGOV", "INDPROD", "EMP_PCG", "HOUSING",
            "CPI", "CORECPI", "PCE", "COREPCE", "UNEMP", "TBOND", "BAABOND"),
  prefix = c("drgdp", "drconsum", "drnresin", "drresinv", "drfedgov", "dindprod", "dempb", "dhousing",
             "CPI", "CORECPI", "PCE", "COREPCE", "UNEMP", "TBOND", "BAABOND"))
scenario_ids <- c("jpm_midyear_2026", "gs_inflation_oil_2026")
table_vars <- c("rgdp", "rpce", "rnresinv", "payems", "urate", "pce", "corepce", "cpi", "corecpi",
                "gs2", "gs10", "baa", "houst", "sp500", "wti")
table_h <- 12

#' Our forecast in the SPF's measures: exact annualized QoQ growth for log
#' series, levels for rates.
spf_measure <- function(paths, history) {
  m <- qoq_annualized_exact(paths, history, is_log)
  m[, , !is_log] <- paths[, , !is_log]
  m
}

#' Write the figures from saved results.
write_figures <- function(res) {
  units <- setNames(ifelse(is_log[res$spf_map$variable], "QoQ ann. %", "level"), res$spf_map$variable)
  for (m in names(res$methods)) {
    png(file.path(out_dir, paste0("spf_", m, ".png")), 1700, 1250, res = 130)
    plot_spf_panels(res$methods[[m]]$spf, res$spf_map$variable, units,
                    paste0("Model baseline (", m, ") against the SPF median, 2026Q3 survey"))
    dev.off()
  }
  png(file.path(out_dir, "long_horizon.png"), 1700, 1350, res = 130)
  plot_long_horizon(lapply(res$methods, `[[`, "long"), main_vars[main_vars != "houst"], res$ranges,
                    "Baseline to 40 quarters: YoY growth (log series) or level (rates)")
  dev.off()
}

# ---- Stage: run ----
if (stage == "run") {
  runs <- list.files(cfg$output_root, pattern = "^[0-9]{4}m[0-9]{2}_.+_full_[0-9a-f]{8}$", full.names = TRUE)
  methods <- vapply(runs, function(d) jsonlite::read_json(file.path(d, "meta.json"))$covid$method, character(1))
  spf <- read_spf_survey(spf_survey$level, spf_survey$growth, spf_survey$year, spf_survey$quarter, spf_map)
  res <- list(time = Sys.time(), spf_map = spf_map, spf_survey = spf_survey, methods = list())
  for (k in seq_along(runs)) {
    m <- methods[[k]]; d <- runs[k]
    t0 <- proc.time()[["elapsed"]]
    dr <- readRDS(file.path(d, "draws.rds"))
    hist <- est_data_to_matrix(readRDS(file.path(d, "covid.rds"))$origin_data, series)
    base <- readRDS(file.path(d, "engine", "baseline.rds"))
    H <- dim(base$paths)[2]

    # The engine's seed must reproduce its baseline (the noise run then
    # differs only in the seed).
    again <- forecast_unconditional(dr$beta, dr$sigma, hist, lags, H, engine_seed)
    reproduced <- max(abs(again - base$paths))
    alt   <- forecast_unconditional(dr$beta, dr$sigma, hist, lags, H, noise_seed)
    long  <- forecast_unconditional(dr$beta, dr$sigma, hist, lags, long_h, engine_seed)

    annual      <- growth_measure(base$paths, hist, is_log, "annual")
    annual_long <- growth_measure(long, hist, is_log, "annual")
    ranges <- historical_ranges(hist, is_log, "2019Q4")

    # 5. Tilted against untilted: the anchor scenario's paths and weights.
    anc <- readRDS(file.path(d, "engine", "scenario_anchor_pce_2pct.rds"))
    anc_annual <- growth_measure(anc$paths, hist, is_log, "annual")
    tilt <- rbind(cbind(forecast = "untilted", measure_quantiles(anc_annual, main_vars, c(4, 8, 12, 20))),
                  cbind(forecast = "tilted", measure_quantiles(anc_annual, main_vars, c(4, 8, 12, 20),
                                                               anc$tilt$weights)))

    # 6. Scenario effects (scenario minus baseline) and the scenario
    # forecasts themselves, in the annual measure. A scenario fixes levels,
    # so its effect depends on the baseline it is measured against; the
    # scenario forecast does not.
    scen <- do.call(rbind, lapply(scenario_ids, function(id) {
      s <- readRDS(file.path(d, "engine", paste0("scenario_", id, ".rds")))
      eff <- irf_annual_measure(s$effect, is_log)   # effect is 0 before the forecast starts
      rbind(cbind(scenario = id, quantity = "effect", measure_quantiles(eff, main_vars, c(2, 4, 6, 8, 12))),
            cbind(scenario = id, quantity = "forecast",
                  measure_quantiles(growth_measure(s$paths, hist, is_log, "annual"), main_vars, c(2, 4, 6, 8, 12))))
    }))

    res$methods[[m]] <- list(
      run_id = basename(d), reproduced = reproduced,
      spf = compare_spf(spf_measure(base$paths, hist), spf),
      ranges_check = range_check(annual_long, ranges, c(4, 8, 12, 20, 40)),
      negative = negative_share(long, yields),
      long = measure_quantiles(annual_long, main_vars, seq_len(long_h)),
      baseline_annual = measure_quantiles(annual, main_vars, c(4, 8, 12, 20)),
      tilt = tilt, tilt_ess_share = anc$tilt$ess_share, tilt_targets = anc$tilt$check,
      scenarios = scen,
      noise_annual = seed_difference(annual, growth_measure(alt, hist, is_log, "annual"), main_vars),
      noise_qoq = seed_difference(qoq_annualized_exact(base$paths, hist, is_log),
                                  qoq_annualized_exact(alt, hist, is_log), main_vars[is_log[main_vars]]),
      table = forecast_table(base$paths[, seq_len(table_h), , drop = FALSE], hist, is_log, table_vars,
                             base$inputs$quarters[seq_len(table_h)]),
      seconds = proc.time()[["elapsed"]] - t0)
    res$ranges <- ranges
    cat(m, ": baseline reproduced to", signif(reproduced, 3), "; done in",
        round(res$methods[[m]]$seconds), "s\n")
  }
  method_order <- c("none", "precovid_params", "exclude_window")
  res$methods <- res$methods[method_order[method_order %in% names(res$methods)]]
  saveRDS(res, file.path(out_dir, "forecast_review.rds"))
  tbl <- do.call(rbind, lapply(names(res$methods), function(m) cbind(method = m, res$methods[[m]]$table)))
  write.csv(tbl, file.path(out_dir, "baseline_table.csv"), row.names = FALSE)
  write_figures(res)
}

# ---- Stage: report ----
if (stage == "report") {
  write_figures(readRDS(file.path(out_dir, "forecast_review.rds")))
  out_abs <- normalizePath(out_dir, winslash = "/")
  quarto::quarto_render("reports/forecast_review.qmd", execute_params = list(output_dir = out_abs), quiet = TRUE)
  target <- file.path(out_dir, "forecast_review.html")
  if (file.exists(target)) file.remove(target)
  if (!file.rename("reports/forecast_review.html", target)) stop("could not move the report to ", target)
  cat("Wrote", target, "\n")
}
