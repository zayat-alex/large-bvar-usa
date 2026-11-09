# check_prior_sensitivity.R: sensitivity of forecasts and impulse responses
# to the hyperparameters (R/prior_sensitivity_checks.R).
#
# Run from the project root (so .Rprofile activates renv):
#   Rscript scripts/check_prior_sensitivity.R run
#   Rscript scripts/check_prior_sensitivity.R run quick   (50 draws; tests the script)
#   Rscript scripts/check_prior_sensitivity.R report
#
# Model: COVID method none, the current full run (its mode and posterior
# draws), the full sample 1973Q1 to the last quarter. For the mode and for
# each variant, 1,000 draws of (B, Sigma) given the hyperparameters, then
# forecasts from the last quarter and impulse responses (recursive gs2
# shock, GIRF of real GDP), up to 20 quarters.
#
# Seeds: coefficient draws mcmc$seed + 4000, forecast shocks mcmc$seed +
# 4001, the same for the mode and every variant, so that they differ only
# through the hyperparameters; the noise rerun at the mode uses + 4100 and
# + 4101 (as in scripts/check_hyper_stability.R).
#
# Writes outputs/checks_prior_sensitivity/ (prior_sensitivity.rds, PNG
# figures, prior_sensitivity.html).

for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)

args  <- commandArgs(trailingOnly = TRUE)
stage <- args[1]
quick <- isTRUE(args[2] == "quick")
if (is.na(stage)) stop("give a stage: run or report")

# ---- Setup ----
cfg     <- yaml::read_yaml("config/run.yml")
cfg$covid$method <- "none"
series  <- read_series_config("config/series.csv")
lags    <- cfg$model$lags
out_dir <- "outputs/checks_prior_sensitivity"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_draws   <- if (quick) 50 else 1000
horizon   <- 20
windows   <- list(`1-12` = 1:12, `13-20` = 13:20)
main_vars <- c("rgdp", "pce", "corepce", "urate", "gs2", "gs10")
irf_resp  <- c("rgdp", "urate", "corepce", "gs10", "baa", "sp500")
seeds     <- cfg$mcmc$seed + c(coef = 4000, shocks = 4001, coef_alt = 4100, shocks_alt = 4101)

#' Write the figures from saved results.
write_figures <- function(res) {
  png(file.path(out_dir, "sensitivity_forecast.png"), 1500, 950, res = 130)
  par(mfrow = c(2, 3), mar = c(4, 3, 2.5, 1))
  for (v in main_vars) plot_sensitivity_panel(res$forecast, v, paste(v, "forecast"))
  dev.off()
  png(file.path(out_dir, "sensitivity_irf.png"), 1500, 950, res = 130)
  par(mfrow = c(2, 3), mar = c(4, 3, 2.5, 1))
  for (v in irf_resp) plot_sensitivity_panel(res$irf_recursive, v, paste(v, "to recursive gs2 shock"))
  dev.off()
}

# ---- Stage: run ----
if (stage == "run") {
  t_start <- proc.time()[["elapsed"]]

  # Step 1: data, setup and the current run.
  model_data <- readRDS(file.path(cfg$data$output_dir, "data.rds"))$model_data
  covid <- covid_treatment(model_data, cfg)
  y <- est_data_to_matrix(covid$est_data, series)
  setup <- estimation_setup(y, series, cfg$prior, lags)
  d <- list.files(cfg$output_root, pattern = "^[0-9]{4}m[0-9]{2}_none_full_[0-9a-f]{8}$",
                  full.names = TRUE)
  if (length(d) != 1) stop("expected one full run of method none; found ", length(d))
  est_now <- readRDS(file.path(d, "estimation.rds"))
  draws_now <- readRDS(file.path(d, "draws.rds"))
  theta_mode <- est_now$mode$theta_mode
  sd_theta <- apply(log(draws_now$hyper[, names(theta_mode)]), 2, sd)
  variants <- sensitivity_variants(theta_mode, sd_theta)
  info <- attr(variants, "info")
  # How far each variant is from the mode in posterior density.
  lp_mode <- log_post_theta(theta_mode, setup)
  info$log_post_vs_mode <- sapply(variants, log_post_theta, setup = setup) - lp_mode
  info$hyper_value <- sapply(seq_along(variants), function(k) {
    nm <- c(lambda = "lambda", mu = "soc", `psi (all)` = NA)[[info$hyperparameter[k]]]
    if (is.na(nm)) NA_real_ else exp(variants[[k]][[nm]])
  })
  print(info, digits = 3)

  # Step 2: draws, forecasts and IRFs at the mode, the noise rerun, each variant.
  rec <- cfg$irf$recursive
  one <- function(theta, seed_coef, seed_shocks) {
    dr <- draws_at_mode(theta, setup, n_draws, seed_coef)
    list(fc = forecast_unconditional(dr$beta, dr$sigma, y, lags, horizon, seed_shocks),
         irf = irf_draws(dr$beta, dr$sigma, lags, horizon, series$short_name, "rgdp",
                         list(slow = unlist(rec$slow), policy = rec$policy, fast = unlist(rec$fast))))
  }
  ref  <- one(theta_mode, seeds[["coef"]], seeds[["shocks"]])
  seed <- one(theta_mode, seeds[["coef_alt"]], seeds[["shocks_alt"]])
  last <- y[nrow(y), ]
  vars <- series$short_name
  # Growth rates and percentage points (R/hyper_stability_checks.R): YoY and
  # QoQ annualized growth for log series, levels for raw series; IRFs as the
  # effect on YoY growth. Yardstick: each variable's history to 2019Q4.
  is_log <- series$transform == "log"
  scale <- historical_scale(y, is_log, "2019Q4")
  gm <- function(p, m) growth_measure(p, y, is_log, m)
  ia <- function(i) irf_annual_measure(i, is_log)
  ref_g  <- list(a = gm(ref$fc, "annual"), q = gm(ref$fc, "quarterly"),
                 ir = ia(ref$irf$recursive_gs2), ig = ia(ref$irf$girf_rgdp))
  seed_g <- list(a = gm(seed$fc, "annual"), q = gm(seed$fc, "quarterly"),
                 ir = ia(seed$irf$recursive_gs2), ig = ia(seed$irf$girf_rgdp))

  forecast <- irf_recursive <- irf_girf <- growth <- list()
  for (nm in names(variants)) {
    alt <- one(variants[[nm]], seeds[["coef"]], seeds[["shocks"]])
    forecast[[nm]] <- compare_bands(ref$fc, alt$fc, seed$fc, vars, seq_len(horizon), base = last)
    irf_recursive[[nm]] <- compare_bands(ref$irf$recursive_gs2, alt$irf$recursive_gs2,
                                         seed$irf$recursive_gs2, vars, 0:horizon)
    irf_girf[[nm]] <- compare_bands(ref$irf$girf_rgdp, alt$irf$girf_rgdp, seed$irf$girf_rgdp,
                                    vars, 0:horizon)
    growth[[nm]] <- list(
      annual = compare_medians(ref_g$a, gm(alt$fc, "annual"), seed_g$a, vars, seq_len(horizon)),
      quarterly = compare_medians(ref_g$q, gm(alt$fc, "quarterly"), seed_g$q, vars, seq_len(horizon)),
      irf_recursive = compare_medians(ref_g$ir, ia(alt$irf$recursive_gs2), seed_g$ir, vars, 0:horizon),
      irf_girf_rgdp = compare_medians(ref_g$ig, ia(alt$irf$girf_rgdp), seed_g$ig, vars, 0:horizon))
    cat(nm, "done\n")
  }

  res <- list(
    run = list(time = Sys.time(), quick = quick, current_run = basename(d),
               sample = paste(rownames(y)[1], "to", rownames(y)[nrow(y)]), n_draws = n_draws,
               horizon = horizon, windows = windows, seeds = seeds,
               seconds = proc.time()[["elapsed"]] - t_start),
    info = info, theta_mode = theta_mode, sd_theta = sd_theta,
    forecast = forecast, irf_recursive = irf_recursive, irf_girf_rgdp = irf_girf,
    growth = growth, scale = scale)
  saveRDS(res, file.path(out_dir, if (quick) "prior_sensitivity_quick.rds" else "prior_sensitivity.rds"))
  if (!quick) write_figures(res)

  # Console summary: forecasts of the main variables, quarters 1-12.
  for (nm in names(variants)) {
    s <- summarize_bands(forecast[[nm]][forecast[[nm]]$variable %in% main_vars, ], windows["1-12"])
    cat(sprintf("%-14s forecast max|med diff| %.3f (noise %.3f), width ratio %.2f-%.2f\n", nm,
                max(s$max_abs_median_diff), max(s$seed_max_abs_median_diff),
                min(s$width_ratio_min), max(s$width_ratio_max)))
  }
  cat("done in", round(res$run$seconds / 60, 1), "minutes\n")
}

# ---- Stage: report ----
if (stage == "report") {
  write_figures(readRDS(file.path(out_dir, "prior_sensitivity.rds")))
  out_abs <- normalizePath(out_dir, winslash = "/")
  quarto::quarto_render("reports/prior_sensitivity.qmd", execute_params = list(output_dir = out_abs),
                        quiet = TRUE)
  target <- file.path(out_dir, "prior_sensitivity.html")
  if (file.exists(target)) file.remove(target)
  if (!file.rename("reports/prior_sensitivity.html", target)) stop("could not move the report to ", target)
  cat("Wrote", target, "\n")
}
