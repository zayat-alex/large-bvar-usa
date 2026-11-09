# check_recovery.R: simulation recovery check (R/recovery_checks.R): simulate
# data from a known VAR many times, estimate the model on each sample, and
# record how often the posterior bands contain the truth.
#
# Run from the project root (so .Rprofile activates renv):
#   Rscript scripts/check_recovery.R run          (100 replications)
#   Rscript scripts/check_recovery.R run quick    (3 replications; tests the script)
#   Rscript scripts/check_recovery.R report
#
# Model: 6 of the 31 variables (as in scripts/check_priors.R), 5 lags,
# the prior of config/run.yml, short chains. The true model is the
# posterior mean of B and Sigma from this small model estimated on the data
# (full sample, COVID method none); the simulations start from the data's
# first 5 quarters and have as many quarters as the data.
#
# Replications run in parallel (mcmc$n_cores R processes); results do not
# depend on the number of processes. Seeds: the true model uses seed + 3500;
# replication k uses seed + 3500 + 10 k and the next two integers.
#
# Writes outputs/checks_recovery/ (recovery_checks.rds, PNG figures,
# recovery_checks.html).

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
seed    <- cfg$mcmc$seed + 3500
out_dir <- "outputs/checks_recovery"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

small_vars <- c("rgdp", "corepce", "urate", "gs2", "sp500", "rvol")
ser_small  <- series[match(small_vars, series$short_name), ]
# Recursive identification in the small model, as in the paper's
# classification (config irf$recursive): economic variables, then gs2, then
# financial variables.
recursive  <- list(slow = c("rgdp", "corepce", "urate"), policy = "gs2", fast = c("sp500", "rvol"))
n_reps      <- if (quick) 3 else 100
horizon     <- 8           # forecast quarters checked (the site shows 8 to 12)
irf_horizon <- 12
mcmc        <- list(n_draw = 6000, n_burn = 1000, n_thin = 5)   # 1,000 retained draws per fit

#' Write the figures from saved results (also redrawn by the report stage).
write_figures <- function(res) {
  vars <- res$run$variables
  png(file.path(out_dir, "recovery_forecast_coverage.png"), 1500, 900, res = 130)
  plot_coverage(res$forecast, res$run$n_ok, seq_len(res$run$horizon), "quarters ahead", vars)
  legend_coverage()
  dev.off()
  png(file.path(out_dir, "recovery_irf_coverage.png"), 1500, 900, res = 130)
  plot_coverage(res$irf, res$run$n_ok, 0:res$run$irf_horizon, "quarters after impact", vars)
  legend_coverage()
  dev.off()
  png(file.path(out_dir, "recovery_example_fan.png"), 1500, 900, res = 130)
  plot_recovery_fan(res$example, vars)
  dev.off()
}

# ---- Stage: run ----
if (stage == "run") {
  t_start <- proc.time()[["elapsed"]]

  # The true model: posterior mean of the small model on the data.
  model_data <- readRDS(file.path(cfg$data$output_dir, "data.rds"))$model_data
  covid <- covid_treatment(model_data, cfg)
  y <- est_data_to_matrix(covid$est_data[, c("quarter", small_vars)], ser_small)
  fit0 <- estimate_bvar(estimation_setup(y, ser_small, cfg$prior, lags),
                        mcmc$n_draw, mcmc$n_burn, mcmc$n_thin, seed = seed)
  truth <- list(B = apply(fit0$beta, c(2, 3), mean), Sigma = apply(fit0$sigma, c(2, 3), mean),
                init = y[seq_len(lags), , drop = FALSE])
  dimnames(truth$Sigma) <- list(small_vars, small_vars)
  colnames(truth$B) <- small_vars
  truth$irf <- recursive_irf(truth$B, truth$Sigma, lags, irf_horizon, small_vars,
                             recursive$slow, recursive$policy, recursive$fast)
  truth$max_root <- companion_max_modulus(truth$B, lags)
  cat("true model: largest root", signif(truth$max_root, 6), "\n")
  if (truth$max_root >= 1) stop("the true model is not stable; the simulated data would explode")

  # Replications, in parallel.
  n_workers <- min(cfg$mcmc$n_cores, n_reps)
  cl <- parallel::makePSOCKcluster(n_workers)
  r_files <- normalizePath(list.files("R", pattern = "[.]R$", full.names = TRUE))
  parallel::clusterCall(cl, function(lib, files) {
    .libPaths(lib)
    for (f in files) source(f)
    NULL
  }, .libPaths(), r_files)
  reps <- parallel::clusterApplyLB(cl, seq_len(n_reps), recovery_replication, truth = truth,
                                   series = ser_small, prior_cfg = cfg$prior, lags = lags,
                                   n_sample = nrow(y), horizon = horizon, mcmc = mcmc,
                                   recursive = recursive, irf_horizon = irf_horizon, seed = seed)
  parallel::stopCluster(cl)

  ok <- vapply(reps, function(r) r$ok, logical(1))
  cat("replications:", sum(ok), "of", n_reps, "fitted\n")
  good <- reps[ok]
  res <- list(
    run = list(time = Sys.time(), quick = quick, n_reps = n_reps, n_ok = sum(ok), seed = seed,
               variables = small_vars, lags = lags, n_sample = nrow(y),
               sample = paste(rownames(y)[1], "to", rownames(y)[nrow(y)]), horizon = horizon,
               irf_horizon = irf_horizon, mcmc = mcmc, recursive = recursive, n_workers = n_workers,
               seconds = proc.time()[["elapsed"]] - t_start),
    truth = truth,
    failures = lapply(reps[!ok], function(r) r[c("rep", "error")]),
    forecast = coverage_rates(good, "forecast"),
    irf = coverage_rates(good, "irf"),
    parameters = parameter_coverage_table(good, truth, lags),
    per_rep = data.frame(rep = sapply(good, `[[`, "rep"), seconds = sapply(good, `[[`, "seconds"),
                         acceptance = sapply(good, `[[`, "acceptance"),
                         lambda_median = sapply(good, `[[`, "lambda_median"),
                         max_root_mean = sapply(good, `[[`, "max_root_mean")),
    example = good[[1]])
  dimnames(res$forecast$c68) <- dimnames(res$forecast$c90) <- list(paste0("h", seq_len(horizon)), small_vars)
  dimnames(res$irf$c68) <- dimnames(res$irf$c90) <- list(paste0("h", 0:irf_horizon), small_vars)
  saveRDS(res, file.path(out_dir, "recovery_checks.rds"))

  write_figures(res)

  cat("forecast coverage, all cells: 68%:", round(mean(res$forecast$c68), 3),
      " 90%:", round(mean(res$forecast$c90), 3), "\n")
  cat("IRF coverage, horizons 1..", irf_horizon, ": 68%:", round(mean(res$irf$c68[-1, ]), 3),
      " 90%:", round(mean(res$irf$c90[-1, ]), 3), "\n")
  print(res$parameters, digits = 3)
  cat("done in", round(res$run$seconds / 60, 1), "minutes\n")
}

# ---- Stage: report ----
if (stage == "report") {
  write_figures(readRDS(file.path(out_dir, "recovery_checks.rds")))
  out_abs <- normalizePath(out_dir, winslash = "/")
  quarto::quarto_render("reports/recovery_checks.qmd", execute_params = list(output_dir = out_abs),
                        quiet = TRUE)
  target <- file.path(out_dir, "recovery_checks.html")
  if (file.exists(target)) file.remove(target)
  if (!file.rename("reports/recovery_checks.html", target)) stop("could not move the report to ", target)
  cat("Wrote", target, "\n")
}
