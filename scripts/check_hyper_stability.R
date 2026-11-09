# check_hyper_stability.R: hyperparameter stability at an earlier origin
# (R/hyper_stability_checks.R). Estimates the hyperparameter mode on data up
# to the origin, compares it with the current mode, and compares forecasts
# and impulse responses drawn under the two.
#
# Run from the project root (so .Rprofile activates renv):
#   Rscript scripts/check_hyper_stability.R run 2016Q2   (one origin)
#   Rscript scripts/check_hyper_stability.R run 2016Q2 quick
#     (tests the script: 5 optimizer iterations, 50 draws, the current
#     Hessian in place of the origin's; the numbers are not the check)
#   Rscript scripts/check_hyper_stability.R compare      (redo the comparisons of every
#     origin from its saved mode, adding growth rates; no mode search)
#   Rscript scripts/check_hyper_stability.R report       (every origin run so far)
#
# Model: COVID method none, 31 variables, 5 lags, the prior of
# config/run.yml. The current run is the full run of method none in
# outputs/ (its estimation.rds and draws.rds). Each origin's sample is
# 1973Q1 to the origin.
#
# Cost per origin: the mode search (three starts in parallel R processes,
# so as long as the slowest start), the Hessian at the mode (several
# minutes at full size), then 3 x 1,000 coefficient draws, forecasts and
# impulse responses (about a minute).
#
# Seeds: coefficient draws use mcmc$seed + 4000, forecast shocks
# mcmc$seed + 4001; the same for both sets of hyperparameters, so that they
# differ only through the hyperparameters. The simulation-noise comparison
# uses mcmc$seed + 4100 and + 4101.
#
# Writes outputs/checks_hyper_stability/ (origin_<origin>.rds, PNG figures,
# hyper_stability.html).

for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)

args  <- commandArgs(trailingOnly = TRUE)
stage <- args[1]
quick <- isTRUE(args[3] == "quick")
if (is.na(stage)) stop("give a stage: run <origin> or report")

# ---- Setup ----
cfg     <- yaml::read_yaml("config/run.yml")
cfg$covid$method <- "none"
series  <- read_series_config("config/series.csv")
lags    <- cfg$model$lags
out_dir <- "outputs/checks_hyper_stability"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_draws   <- if (quick) 50 else 1000
horizon   <- 20               # forecast and IRF horizon compared (5 years)
windows   <- list(`1-12` = 1:12, `13-20` = 13:20)   # 1-3 years is the main window
main_vars <- c("rgdp", "pce", "corepce", "urate", "gs2", "gs10")
irf_resp  <- c("rgdp", "urate", "corepce", "gs10", "baa", "sp500")
seeds     <- cfg$mcmc$seed + c(coef = 4000, shocks = 4001, coef_alt = 4100, shocks_alt = 4101)

#' The current full run of method none (the one folder matching the pattern).
current_run_dir <- function() {
  d <- list.files(cfg$output_root, pattern = "^[0-9]{4}m[0-9]{2}_none_full_[0-9a-f]{8}$",
                  full.names = TRUE)
  if (length(d) != 1) stop("expected one full run of method none in ", cfg$output_root,
                           "; found ", length(d))
  d
}

#' Write the figures of one origin from its saved results.
write_figures <- function(res) {
  o <- res$run$origin
  png(file.path(out_dir, paste0("shift_", o, ".png")), 1100, 1300, res = 130)
  plot_hyper_shift(res$shift, o)
  dev.off()
  png(file.path(out_dir, paste0("forecast_", o, ".png")), 1500, 950, res = 130)
  par(mfrow = c(2, 3), mar = c(4, 3, 2.5, 1))
  for (v in main_vars) plot_band_pair(res$paths$ref, res$paths$alt, v, seq_len(res$run$horizon),
                                      res$paths$last[[v]], paste0(v, " (change since ", o, ")"), TRUE)
  dev.off()
  png(file.path(out_dir, paste0("irf_", o, ".png")), 1500, 950, res = 130)
  par(mfrow = c(2, 3), mar = c(4, 3, 2.5, 1))
  for (v in irf_resp) plot_band_pair(res$paths$irf_ref, res$paths$irf_alt, v, 0:res$run$horizon,
                                     0, paste(v, "to recursive gs2 shock"), FALSE)
  dev.off()
}

#' Draw coefficients under the current and the origin's mode on the origin's
#' sample, and compare forecasts and impulse responses: against the band
#' (compare_bands()) and in growth rates and percentage points against each
#' variable's history (compare_medians(), historical_scale()).
#'
#' @param theta_now,theta_past Modes on the log scale.
#' @param setup estimation_setup() of the origin's sample.
#' @param y The origin's sample (quarters x n).
#' @return List: forecast, irf_recursive, irf_girf_rgdp (compare_bands()),
#'   growth (scale, annual, quarterly, irf_recursive, irf_girf_rgdp), max_root, paths.
draw_and_compare <- function(theta_now, theta_past, setup, y) {
  d_ref  <- draws_at_mode(theta_now, setup, n_draws, seeds[["coef"]])
  d_alt  <- draws_at_mode(theta_past, setup, n_draws, seeds[["coef"]])
  d_seed <- draws_at_mode(theta_now, setup, n_draws, seeds[["coef_alt"]])
  fc <- function(d, s) forecast_unconditional(d$beta, d$sigma, y, lags, horizon, s)
  f_ref <- fc(d_ref, seeds[["shocks"]]); f_alt <- fc(d_alt, seeds[["shocks"]])
  f_seed <- fc(d_seed, seeds[["shocks_alt"]])
  rec <- cfg$irf$recursive
  ir <- function(d) irf_draws(d$beta, d$sigma, lags, horizon, series$short_name, "rgdp",
                              list(slow = unlist(rec$slow), policy = rec$policy,
                                   fast = unlist(rec$fast)))
  i_ref <- ir(d_ref); i_alt <- ir(d_alt); i_seed <- ir(d_seed)
  last <- y[nrow(y), ]
  vars <- series$short_name
  is_log <- series$transform == "log"
  roots <- function(d) apply(d$beta, 1, companion_max_modulus, lags = lags)

  # Growth rates and percentage points, path by path.
  gm <- function(p, m) growth_measure(p, y, is_log, m)
  ia <- function(i) irf_annual_measure(i, is_log)
  growth <- list(
    scale = historical_scale(y_full, is_log, scale_end),
    annual = compare_medians(gm(f_ref, "annual"), gm(f_alt, "annual"), gm(f_seed, "annual"),
                             vars, seq_len(horizon)),
    quarterly = compare_medians(gm(f_ref, "quarterly"), gm(f_alt, "quarterly"),
                                gm(f_seed, "quarterly"), vars, seq_len(horizon)),
    irf_recursive = compare_medians(ia(i_ref$recursive_gs2), ia(i_alt$recursive_gs2),
                                    ia(i_seed$recursive_gs2), vars, 0:horizon),
    irf_girf_rgdp = compare_medians(ia(i_ref$girf_rgdp), ia(i_alt$girf_rgdp),
                                    ia(i_seed$girf_rgdp), vars, 0:horizon))

  list(
    forecast = compare_bands(f_ref, f_alt, f_seed, vars, seq_len(horizon), base = last),
    irf_recursive = compare_bands(i_ref$recursive_gs2, i_alt$recursive_gs2, i_seed$recursive_gs2,
                                  vars, 0:horizon),
    irf_girf_rgdp = compare_bands(i_ref$girf_rgdp, i_alt$girf_rgdp, i_seed$girf_rgdp, vars, 0:horizon),
    growth = growth,
    max_root = list(ref = roots(d_ref), alt = roots(d_alt)),
    paths = list(ref = f_ref[, , main_vars], alt = f_alt[, , main_vars], last = last[main_vars],
                 irf_ref = i_ref$recursive_gs2[, , irf_resp], irf_alt = i_alt$recursive_gs2[, , irf_resp]))
}

# Data and the historical yardstick (full sample up to scale_end, the same
# for every origin so that origins are compared on one scale).
model_data <- readRDS(file.path(cfg$data$output_dir, "data.rds"))$model_data
covid <- covid_treatment(model_data, cfg)
y_full <- est_data_to_matrix(covid$est_data, series)
scale_end <- "2019Q4"

# ---- Stage: run ----
if (stage == "run") {
  origin <- args[2]
  if (is.na(origin)) stop("give the origin quarter, e.g. run 2016Q2")
  t_start <- proc.time()[["elapsed"]]

  # Step 1: sample up to the origin and the current run.
  y <- est_data_to_matrix(sample_to_origin(covid$est_data, origin), series)
  setup <- estimation_setup(y, series, cfg$prior, lags)
  run_dir <- current_run_dir()
  est_now <- readRDS(file.path(run_dir, "estimation.rds"))
  draws_now <- readRDS(file.path(run_dir, "draws.rds"))
  cat("origin", origin, ": sample", rownames(y)[1], "to", rownames(y)[nrow(y)],
      "(", nrow(y), "quarters ); current run", basename(run_dir), "\n")

  # Step 2: mode (starts in parallel) and Hessian at the origin.
  mode_past <- find_hyper_mode_parallel(setup, unlist(cfg$mcmc$mode_search$start_factors),
                                        if (quick) 5 else cfg$mcmc$mode_search$maxit)
  cat("mode search:", round(mode_past$seconds / 60, 1), "min; all starts converged:",
      all(mode_past$table$convergence == 0), "\n")
  print(mode_past$table[, c("start_factor", "convergence", "log_post_theta", "lambda", "soc")])
  prop_past <- if (quick) modifyList(est_now$proposal, list(seconds = 0)) else
    timed_call(hyper_proposal, mode_past$theta_mode, setup)
  cat("Hessian:", round(prop_past$seconds / 60, 1), "min\n")

  # Step 3: the modes compared.
  shift <- hyper_shift_table(mode_past, prop_past, est_now$mode, est_now$proposal,
                             draws_now$hyper, setup, series$short_name)

  # Step 4: coefficients on the origin's sample under each mode; forecasts
  # from the origin and impulse responses, compared.
  res <- c(list(
    run = list(time = Sys.time(), quick = quick, origin = origin, current_run = basename(run_dir),
               sample = paste(rownames(y)[1], "to", rownames(y)[nrow(y)]), n_sample = nrow(y),
               n_draws = n_draws, horizon = horizon, windows = windows, seeds = seeds,
               mode_seconds = mode_past$seconds, hessian_seconds = prop_past$seconds),
    mode_table = mode_past$table, mode_past = mode_past, proposal_past = prop_past,
    neg_hessian_min_eigen = min(prop_past$neg_hessian_eigenvalues),
    shift = shift),
    draw_and_compare(est_now$mode$theta_mode, mode_past$theta_mode, setup, y))
  res$run$seconds <- proc.time()[["elapsed"]] - t_start
  if (quick) res$run$origin <- paste0(origin, "_quick")
  saveRDS(res, file.path(out_dir, paste0("origin_", res$run$origin, ".rds")))
  write_figures(res)

  # Console summary.
  cat("\nlargest |shift| in current posterior s.d.:",
      round(max(abs(shift$shift_in_sd_now)), 2), "(", shift$hyperparameter[which.max(abs(shift$shift_in_sd_now))], ")\n")
  print(shift[1:2, c("hyperparameter", "mode_now", "mode_past", "log_ratio", "shift_in_sd_now",
                     "shift_in_sd_past")], digits = 3)
  cat("\nforecasts, main variables:\n")
  print(summarize_bands(res$forecast[res$forecast$variable %in% main_vars, ], windows), digits = 2)
  cat("\nrecursive gs2 IRF:\n")
  print(summarize_bands(res$irf_recursive[res$irf_recursive$variable %in% irf_resp, ],
                        list(`0-12` = 0:12, `13-20` = 13:20)), digits = 2)
  cat("done in", round(res$run$seconds / 60, 1), "minutes\n")
}

# ---- Stage: compare ----
# Redo step 4 for every origin already run, from its saved mode: no new
# mode search. The draws are the same as in the run (same seeds), so the
# band comparison is reproduced and the growth comparison added.
if (stage == "compare") {
  run_dir <- current_run_dir()
  est_now <- readRDS(file.path(run_dir, "estimation.rds"))
  files <- list.files(out_dir, pattern = "^origin_.*[.]rds$", full.names = TRUE)
  files <- files[!grepl("_quick[.]rds$", files)]
  for (f in files) {
    t0 <- proc.time()[["elapsed"]]
    res <- readRDS(f)
    y <- est_data_to_matrix(sample_to_origin(covid$est_data, res$run$origin), series)
    setup <- estimation_setup(y, series, cfg$prior, lags)
    new <- draw_and_compare(est_now$mode$theta_mode, res$mode_past$theta_mode, setup, y)
    # The band comparison must come out as in the run (same draws).
    same <- isTRUE(all.equal(new$forecast, res$forecast)) &&
      isTRUE(all.equal(new$irf_recursive, res$irf_recursive))
    cat(res$run$origin, ": band comparison reproduced:", same, "\n")
    if (!same) stop("the redrawn comparison differs from the saved one for ", res$run$origin)
    res[names(new)] <- new
    res$run$compare_time <- Sys.time()
    saveRDS(res, f)
    s <- summarize_medians(res$growth$annual[res$growth$annual$variable %in% main_vars, ],
                           res$growth$scale, windows["1-12"])
    print(s, digits = 2)
    cat("in", round(proc.time()[["elapsed"]] - t0), "s\n")
  }
}

# ---- Stage: report ----
if (stage == "report") {
  files <- list.files(out_dir, pattern = "^origin_.*[.]rds$", full.names = TRUE)
  files <- files[!grepl("_quick[.]rds$", files)]   # quick runs test the script only
  if (!length(files)) stop("no origin has been run yet (quick runs are not reported)")
  for (f in files) write_figures(readRDS(f))
  out_abs <- normalizePath(out_dir, winslash = "/")
  quarto::quarto_render("reports/hyper_stability.qmd", execute_params = list(output_dir = out_abs),
                        quiet = TRUE)
  target <- file.path(out_dir, "hyper_stability.html")
  if (file.exists(target)) file.remove(target)
  if (!file.rename("reports/hyper_stability.html", target)) stop("could not move the report to ", target)
  cat("Wrote", target, "\n")
}
