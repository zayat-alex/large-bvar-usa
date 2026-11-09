# check_conditional.R: compare the two conditional forecast methods on the
# model's posterior draws and write the results and the report.
#
# Run from the project root (so .Rprofile activates renv):
#   Rscript scripts/check_conditional.R checks         (all comparisons)
#   Rscript scripts/check_conditional.R checks quick   (a few draws; tests the script)
#   Rscript scripts/check_conditional.R report
#
# Draws used: the 1,000 engine draws of each COVID method's estimation run
# (outputs/<run_id>/draws.rds, with the run_id of config/run.yml), forecast
# from origin_data (the latest data, 2026Q2). Test scenarios
# (make_test_scenarios(), R/conditional_checks.R) are built around each
# run's median no-shock path.
#
# Comparisons (R/conditional_checks.R):
#   moments    1-2: conditional means and per-quarter variances of both
#              methods, and against the textbook formula on some draws;
#   paths      3 and 5: both methods with the same seed (all draws, one path
#              each), and the Gaussian method with a second seed;
#   simulated  4: many paths of one draw against the exact moments.
#
# Writes outputs/checks_conditional/ (conditional_checks.rds, PNG figures,
# conditional_checks.html).

for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)

args  <- commandArgs(trailingOnly = TRUE)
stage <- args[1]
quick <- isTRUE(args[2] == "quick")
if (is.na(stage)) stop("give a stage: checks or report")

# ---- Setup ----
cfg     <- yaml::read_yaml("config/run.yml")
series  <- read_series_config("config/series.csv")
lags    <- cfg$model$lags
seed    <- cfg$mcmc$seed
horizon <- cfg$forecast$horizon                     # 20 quarters
out_dir <- "outputs/checks_conditional"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Sizes of the comparisons (quick: only to test the script).
sizes <- if (quick) {
  list(moment_draws = 2, textbook_draws = 1, path_draws_main = 20, path_draws_other = 10,
       sim_paths = c(kalman = 200, gaussian = 400))
} else {
  list(moment_draws = 10, textbook_draws = 3, path_draws_main = 1000, path_draws_other = 100,
       sim_paths = c(kalman = 5000, gaussian = 20000))
}
main_scenario <- "mixed"            # compared on all draws; the others on path_draws_other
key_vars <- c("rgdp", "corepce", "urate", "gs2", "gs10", "baa", "sp500", "rvol")
key_horizons <- c(1, 2, 4, 8, 12, 20)
# Held at the baseline at T+1 in the many_cells_h1 test scenario: the
# variables listed as slow in the recursive identification proposal
# (used here only as a set of many cells in one quarter).
hold_vars <- c("rgdp", "rpce", "rresinv", "rnresinv", "rexp", "rimp", "rgov", "rfedgov", "rdpi",
               "gdpdefl", "pce", "corepce", "cpi", "corecpi", "comphr", "payems", "urate",
               "indpro", "caputil", "houst")

#' Draws, forecast history and test scenarios of one COVID method's run.
run_inputs <- function(method) {
  cfg_m <- cfg
  cfg_m$covid$method <- method
  id <- make_run_id(cfg_m, method, mcmc_settings(cfg_m), series)$run_id
  run_dir <- file.path(cfg$output_root, id)
  if (!dir.exists(run_dir)) stop("run folder not found: ", run_dir, " (run scripts/estimate.R)")
  draws <- readRDS(file.path(run_dir, "draws.rds"))
  covid <- readRDS(file.path(run_dir, "covid.rds"))
  history <- est_data_to_matrix(covid$origin_data, series)
  if (!identical(draws$variables, colnames(history))) stop("variable order of draws and history differ")
  no_shock <- forecast_unconditional(draws$beta, draws$sigma, history, lags, horizon, seed = 1,
                                     shocks = FALSE)
  base <- apply(no_shock, c(2, 3), median)
  list(run_id = id, draws = draws, history = history, base = base,
       scenarios = make_test_scenarios(base, hold_vars))
}

#' Evenly spaced draw numbers.
spaced <- function(n_total, n) unique(round(seq(1, n_total, length.out = n)))

# ---- Stage: checks ----
if (stage == "checks") {
  t_start <- proc.time()[["elapsed"]]
  res <- list(run = list(time = Sys.time(), seed = seed, horizon = horizon, lags = lags,
                         quick = quick, sizes = sizes, main_scenario = main_scenario,
                         key_vars = key_vars, key_horizons = key_horizons, hold_vars = hold_vars),
              methods = list())
  for (method in unlist(cfg$covid$methods)) {
    cat("==", method, "\n")
    inp <- run_inputs(method)
    d <- inp$draws
    n_total <- dim(d$beta)[1]
    r <- list(run_id = inp$run_id, origin = rownames(inp$history)[nrow(inp$history)],
              scenarios = lapply(inp$scenarios, function(s) {
                data.frame(fixed_cells = sum(!is.na(s)),
                           variables = paste(colnames(s)[colSums(!is.na(s)) > 0], collapse = ", "),
                           quarters = paste(which(rowSums(!is.na(s)) > 0), collapse = ", "))
              }))

    # 1-2: exact moments.
    idx <- spaced(n_total, sizes$moment_draws)
    r$moments <- compare_moments(d$beta, d$sigma, inp$history, lags, horizon, inp$scenarios,
                                 draw_index = idx, textbook_index = idx[seq_len(sizes$textbook_draws)])
    cat("moments done:", round(proc.time()[["elapsed"]] - t_start), "s\n")

    # 3 and 5: paths with the same seed; quantiles; timing.
    r$paths <- list()
    for (nm in names(inp$scenarios)) {
      n_use <- if (nm == main_scenario) sizes$path_draws_main else sizes$path_draws_other
      use <- spaced(n_total, n_use)
      cp <- compare_paths(d$beta[use, , , drop = FALSE], d$sigma[use, , , drop = FALSE], inp$history,
                          lags, horizon, inp$scenarios[[nm]], seed = seed + 3100, seed_alt = seed + 3101,
                          variables = key_vars, horizons = key_horizons)
      if (nm == main_scenario) {
        baseline <- forecast_unconditional(d$beta[use, , , drop = FALSE], d$sigma[use, , , drop = FALSE],
                                           inp$history, lags, horizon, seed = seed + 3102)
        units <- setNames(ifelse(series$transform[match(key_vars, series$short_name)] == "log",
                                 "log points from 2026Q2", "pp from 2026Q2"), key_vars)
        png(file.path(out_dir, paste0("fans_", method, ".png")), width = 2400, height = 1300, res = 200)
        op <- par(mfrow = c(2, 4), mar = c(4, 4, 2, 1), oma = c(2.5, 0, 0, 0))
        plot_conditional_fans(baseline, cp$paths$kalman, cp$paths$gaussian, inp$scenarios[[nm]],
                              inp$history, key_vars, units)
        legend_conditional_fans()
        par(op)
        dev.off()
      }
      cp$paths <- NULL
      cp$draws <- length(use)
      r$paths[[nm]] <- cp
      cat("paths", nm, "done:", round(proc.time()[["elapsed"]] - t_start), "s\n")
    }

    # 4: simulated against exact, one draw (the middle one), main scenario.
    i_mid <- spaced(n_total, 3)[2]
    r$simulated <- check_simulated_vs_exact(d$beta[i_mid, , ], d$sigma[i_mid, , ], inp$history, lags,
                                            horizon, inp$scenarios[[main_scenario]], sizes$sim_paths,
                                            seed = seed + 3103)
    r$simulated_draw <- i_mid
    cat("simulated done:", round(proc.time()[["elapsed"]] - t_start), "s\n")
    res$methods[[method]] <- r
  }
  res$run$seconds <- proc.time()[["elapsed"]] - t_start
  saveRDS(res, file.path(out_dir, "conditional_checks.rds"))

  # Summary.
  for (method in names(res$methods)) {
    r <- res$methods[[method]]
    cat("\n==", method, "(", r$run_id, ")\n")
    m <- r$moments
    cat("moments: max |mean diff|", signif(max(m$mean_max_abs_diff), 3),
        "; max rel var diff", signif(max(m$var_max_rel_diff), 3),
        "; textbook mean", signif(max(m$textbook_mean_max_abs_diff, na.rm = TRUE), 3),
        "; textbook var", signif(max(m$textbook_var_max_rel_diff, na.rm = TRUE), 3),
        "; fixed var", signif(max(m$fixed_var_max_rel), 3),
        "; var above unc", signif(max(m$var_above_unc_max_rel), 3), "\n")
    for (nm in names(r$paths)) {
      p <- r$paths[[nm]]
      cat(sprintf("paths %-14s draws %4d  max|k-g| %.2e  fixed err %.2e  ms/draw kalman %.0f gaussian %.0f\n",
                  nm, p$draws, p$max_abs_diff, p$fixed_cells_max_error, p$ms_per_draw[["kalman"]],
                  p$ms_per_draw[["gaussian"]]))
    }
    print(r$simulated, digits = 3, row.names = FALSE)
  }
  cat("total", round(res$run$seconds), "s\n")
}

# ---- Stage: report ----
if (stage == "report") {
  out_abs <- normalizePath(out_dir, winslash = "/")
  quarto::quarto_render("reports/conditional_checks.qmd", execute_params = list(output_dir = out_abs),
                        quiet = TRUE)
  target <- file.path(out_dir, "conditional_checks.html")
  if (file.exists(target)) file.remove(target)
  if (!file.rename("reports/conditional_checks.html", target)) stop("could not move the report to ", target)
  cat("Wrote", target, "\n")
}
