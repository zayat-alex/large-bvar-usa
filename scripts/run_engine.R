# run_engine.R: run stage 4 of the pipeline (_targets.R), the forecasting
# engine, and print a summary of each run.
#
# Run from the project root (so .Rprofile activates renv):
#   Rscript scripts/run_engine.R
#
# For each COVID method in config/run.yml (covid$methods), the engine reads
# the saved estimation run (outputs/<run_id>/draws.rds and covid.rds) and
# computes the baseline forecast, every scenario listed in config/run.yml
# (scenarios) and the impulse responses (irf). It does not depend on the
# estimation targets, so it never re-estimates: run scripts/estimate.R
# first if the run does not exist yet, or after changing the model.
# Steps whose inputs have not changed are skipped (cache in _targets/):
# editing one scenario file reruns only that scenario.
#
# Time with 1,000 draws (Gaussian method): about 2-3 minutes per COVID method
# with two hard scenarios, plus about 1 minute per tilted scenario.
#
# Writes <output_root>/<run_id>/engine/ for each run (R/engine_run.R):
# baseline.rds, scenario_<id>.rds, irf.rds and meta.json.

library(targets)
source("R/tilting.R")                            # weighted_quantiles() for the summary
cfg <- yaml::read_yaml("config/run.yml")
methods <- unlist(cfg$covid$methods)

tar_make(names = starts_with("engine_meta_"))

# A few numbers to read at a glance, in the fourth quarter of the second
# forecast year (2027Q4 for data ending in 2026Q2): 4-quarter change of real
# GDP and core PCE prices (log points, the 100 ln levels four quarters apart,
# both in the forecast) and the unemployment rate (percent). Medians over
# paths for the baseline and scenarios (weighted for a tilted scenario); for
# an effect, the median over draws of the per-draw effect (the change of the
# effect over the four quarters for the growth columns). A tilted
# scenario's effect is a single path (weighted mean minus baseline mean).
#
# @param x Array draws x H x n (paths or effects), or a matrix H x n (a
#   tilted effect).
# @param quarters Labels of the forecast quarters.
# @param label Row label.
# @param weights NULL, or tilting weights (one per path).
# @return One-row data.frame.
summary_row <- function(x, quarters, label, weights = NULL) {
  q_end <- grep("Q4$", quarters, value = TRUE)[2]
  h <- match(q_end, quarters)
  if (is.matrix(x)) {
    g4 <- function(v) x[h, v] - x[h - 4, v]
    lev <- x[h, "urate"]
  } else {
    med <- function(v) weighted_quantiles(v, weights, 0.5)
    g4 <- function(v) med(x[, h, v] - x[, h - 4, v])
    lev <- med(x[, h, "urate"])
  }
  data.frame(object = label, quarter = q_end, rgdp_4q = g4("rgdp"), corepce_4q = g4("corepce"), urate = lev)
}

for (m in methods) {
  meta_path <- tar_read_raw(paste0("engine_meta_", m))
  meta <- jsonlite::read_json(meta_path)
  dir <- dirname(meta_path)
  cat("\n== Engine for run", meta$run_id, "==\n")
  cat("Origin", meta$origin_last_quarter, "| forecast", paste(unlist(meta$forecast_quarters), collapse = " to "),
      "|", meta$n_draws, "draws |", meta$settings$paths_per_draw, "path(s) per draw | method",
      meta$settings$method, "| seed", meta$settings$seed, "\n")
  base <- readRDS(file.path(dir, "baseline.rds"))
  rows <- summary_row(base$paths, base$inputs$quarters, "baseline")
  for (s in meta$scenarios) {
    sc <- readRDS(file.path(dir, paste0("scenario_", s$id, ".rds")))
    cat("Scenario", s$id, ":", s$fixed_cells, "fixed cells, largest |path - condition|",
        signif(s$fixed_max_error, 3), ",", round(s$seconds, 1), "s\n")
    if (is.null(sc$tilt)) {
      rows <- rbind(rows, summary_row(sc$paths, base$inputs$quarters, s$id),
                    summary_row(sc$effect, base$inputs$quarters, paste0(s$id, " effect")))
    } else {
      t <- s$tilting
      cat("  tilted:", t$targets, "target(s) on", t$paths, "paths, largest miss", signif(t$largest_miss, 3),
          "| effective sample size", round(t$ess), paste0("(", signif(100 * t$ess_share, 3), "%)"),
          if (isTRUE(t$ess_warning)) "BELOW THE WARNING SHARE" else "", "\n")
      print(sc$tilt$check[, c("variable", "quarter", "measure", "statistic", "target", "untilted_median")],
            row.names = FALSE)
      rows <- rbind(rows, summary_row(sc$paths, base$inputs$quarters, s$id, sc$tilt$weights),
                    summary_row(sc$tilt$effect, base$inputs$quarters, paste0(s$id, " effect (tilted)")))
    }
  }
  cat("Impulse responses:", paste(unlist(meta$irf$identifications), collapse = ", "), "to",
      meta$irf$horizon, "quarters,", round(meta$timings_seconds$irf, 1), "s\n")
  cat("Medians (log points for 4-quarter growth, percent for urate):\n")
  print(rows, digits = 3, row.names = FALSE)
  cat("Files in", dir, "\n")
}
