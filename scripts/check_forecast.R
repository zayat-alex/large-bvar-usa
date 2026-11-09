# check_forecast.R: answer the forecasting questions with the full model and
# write the results and the report.
#
# Run from the project root (so .Rprofile activates renv):
#   Rscript scripts/check_forecast.R draws     (coefficient draws for the checks)
#   Rscript scripts/check_forecast.R checks    (all checks; reads the draws)
#   Rscript scripts/check_forecast.R report
#
# Draws used: for each COVID method, 500 draws of (Sigma, B) from their
# posterior given the hyperparameters at the posterior mode found by
# scripts/check_full_size.R (outputs/checks_full_size/mode_<method>.rds).
# Hyperparameter uncertainty is therefore not included; the forecasting code
# does not depend on where the draws come from.
# That mode search was a check of computing cost at full size, which is not
# included in the public version of this repository; to rerun the `draws`
# stage without it, give it the hyperparameter mode of an estimation run
# (find_hyper_mode(), R/estimate_bvar.R) in place of the mode file.
#
# Checks (R/forecast_checks.R):
#   known     the simulator against exact forecast moments of a known VAR;
#   predict   the simulator against BVAR::predict() with the same draws;
#   history   forecasting from another history: BVAR's newdata, and the
#             pre-2020 coefficients applied to the latest data;
#   cond      cost of BVAR's conditional forecasts at full size;
#   longrun   a first look at 20- and 40-quarter forecasts (preliminary).
#
# Writes outputs/checks_forecast/ (draws_<method>.rds, forecast_checks.rds,
# PNG figures, forecast_checks.html).

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

stage <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(stage)) stop("give a stage: draws, checks or report")

# ---- Setup ----
cfg     <- yaml::read_yaml("config/run.yml")
series  <- read_series_config("config/series.csv")
lags    <- cfg$model$lags
seed    <- cfg$mcmc$seed
horizon <- cfg$forecast$horizon                     # 20 quarters
out_dir <- "outputs/checks_forecast"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
model_data <- readRDS(file.path(cfg$data$output_dir, "data.rds"))$model_data
n_check_draws <- 500

#' COVID treatment output and estimation setup for one method.
method_data <- function(method) {
  cfg_m <- cfg
  cfg_m$covid$method <- method
  covid <- covid_treatment(model_data, cfg_m)
  y <- est_data_to_matrix(covid$est_data, series)
  list(covid = covid, y = y, setup = estimation_setup(y, series, cfg$prior, lags),
       origin = est_data_to_matrix(covid$origin_data, series))
}

# ---- Stage: draws ----
if (stage == "draws") {
  for (method in c("none", "precovid_params")) {
    md   <- method_data(method)
    mode <- readRDS(file.path("outputs/checks_full_size", paste0("mode_", method, ".rds")))
    h    <- complete_hyper(exp(mode$mode$theta_mode), md$setup)
    set.seed(seed + 10)
    d <- draw_coefficients(niw_posterior(h, md$setup), n_check_draws)
    saveRDS(list(method = method, hyper = h, beta = d$beta, sigma = d$sigma,
                 variables = colnames(md$y)),
            file.path(out_dir, paste0("draws_", method, ".rds")))
    cat("Saved", n_check_draws, "draws for", method, "\n")
  }
}

# ---- Stage: checks ----
if (stage == "checks") {
  res <- list(run = list(time = Sys.time(), seed = seed, horizon = horizon, lags = lags,
                         n_draws = n_check_draws, bvar_version = as.character(packageVersion("BVAR"))))
  d_none <- readRDS(file.path(out_dir, "draws_none.rds"))
  d_pre  <- readRDS(file.path(out_dir, "draws_precovid_params.rds"))
  md_none <- method_data("none")
  md_pre  <- method_data("precovid_params")
  key_vars <- c("rgdp", "corepce", "urate", "gs2", "gs10", "rvol")

  # known: a VAR(2) with known coefficients (as in tests/testthat/test-forecast.R).
  B1 <- matrix(c(0.6, 0.1, 0.2, 0.5), 2); B2 <- matrix(c(0.2, 0.0, -0.1, 0.3), 2)
  B  <- rbind(c(0.5, -0.2), t(B1), t(B2))
  S  <- matrix(c(1, 0.4, 0.4, 0.5), 2)
  hist_known <- matrix(c(1, 2, 3, 0.5, 1.5, 1.0), 3, 2, dimnames = list(NULL, c("a", "b")))
  res$known <- check_known_var(B, S, hist_known, lags = 2, horizon = 8, n_sim = 20000, seed = seed + 20)

  # predict: same draws, same history (method none: the full sample).
  hist_none <- md_none$origin
  own_sim   <- forecast_unconditional(d_none$beta, d_none$sigma, hist_none, lags, horizon, seed + 21)
  own_mean  <- forecast_unconditional(d_none$beta, d_none$sigma, hist_none, lags, horizon, seed + 21,
                                      shocks = FALSE)
  shell     <- bvar_shell(d_none$beta, d_none$sigma, hist_none, lags)
  pkg_sim   <- bvar_predict_draws(shell, horizon, seed + 22)
  dimnames(pkg_sim) <- dimnames(own_sim)
  # BVAR's draws minus each draw's no-shock path: the noise BVAR adds.
  pkg_noise <- pkg_sim - own_mean
  # Expected s.d. of that noise if it is t(Sigma) u, u ~ N(0, I): sqrt of diag(Sigma Sigma'),
  # averaged over draws; the project's one-quarter-ahead shock s.d. is sqrt(diag(Sigma)).
  sig2_diag <- apply(d_none$sigma, 1, function(s) diag(s %*% t(s)))      # n x draws
  sig_diag  <- apply(d_none$sigma, 1, diag)
  res$predict <- list(
    spread = compare_spread(own_sim, pkg_sim, key_vars, c(1, 4, 8, 20)),
    noise = data.frame(variable = key_vars,
                       noise_sd_h1 = apply(pkg_noise[, 1, key_vars], 2, sd),
                       noise_sd_h20 = apply(pkg_noise[, horizon, key_vars], 2, sd),
                       sqrt_sigma_squared = sqrt(rowMeans(sig2_diag))[match(key_vars, colnames(hist_none))],
                       sqrt_sigma = sqrt(rowMeans(sig_diag))[match(key_vars, colnames(hist_none))],
                       noise_mean_h20 = apply(pkg_noise[, horizon, key_vars], 2, mean)))
  res$predict$mean_path_max_abs_diff <- max(abs(apply(pkg_sim, c(2, 3), mean) - apply(own_mean, c(2, 3), mean) -
                                                  apply(pkg_noise, c(2, 3), mean)))

  # history: (1) BVAR's newdata equal to the stored data reproduces the default exactly.
  a <- bvar_predict_draws(shell, horizon, seed + 23)
  b <- bvar_predict_draws(shell, horizon, seed + 23, newdata = hist_none)
  res$history$newdata_identical <- identical(a, b)
  # (2) The pre-2020 coefficients forecast from the latest data (2026Q2), and
  # BVAR's newdata gives the same mean paths as the project's simulator.
  hist_pre_origin <- md_pre$origin                 # full sample to 2026Q2
  own_pre_mean <- forecast_unconditional(d_pre$beta, d_pre$sigma, hist_pre_origin, lags, horizon,
                                         seed + 24, shocks = FALSE)
  shell_pre <- bvar_shell(d_pre$beta, d_pre$sigma, md_pre$y, lags)   # stored data end 2019Q4
  pkg_pre   <- bvar_predict_draws(shell_pre, horizon, seed + 25, newdata = hist_pre_origin)
  dimnames(pkg_pre) <- dimnames(own_pre_mean)
  noise_pre <- pkg_pre - own_pre_mean
  res$history$precovid <- data.frame(
    variable = key_vars,
    last_estimation_quarter = md_pre$y[nrow(md_pre$y), key_vars],
    last_observed_quarter = hist_pre_origin[nrow(hist_pre_origin), key_vars],
    own_median_h1 = apply(own_pre_mean[, 1, key_vars], 2, median),
    bvar_newdata_mean_h1_minus_noise = apply(pkg_pre[, 1, key_vars] - noise_pre[, 1, key_vars], 2, mean))
  res$history$labels <- c(rownames(md_pre$y)[nrow(md_pre$y)], rownames(hist_pre_origin)[nrow(hist_pre_origin)])
  # A cleaner check that newdata is used: BVAR's newdata forecast minus the
  # project's no-shock path has mean close to 0 (noise only) at every horizon.
  res$history$precovid_noise_mean_max <- max(abs(apply(noise_pre, c(2, 3), mean)) /
                                               apply(noise_pre, c(2, 3), sd))

  # cond: cost of BVAR's conditional forecasts, 50 draws, horizon 20.
  shell50 <- bvar_shell(d_none$beta[1:50, , , drop = FALSE], d_none$sigma[1:50, , , drop = FALSE],
                        hist_none, lags)
  res$cond <- list(
    light = time_conditional(shell50, c("gs2", "wti", "sp500"), n_cond_q = 8, horizon = horizon, seed = seed + 26),
    heavy = time_conditional(shell50, c("rgdp", "pce", "corepce", "urate", "gs2", "gs10", "baa", "wti",
                                        "sp500", "usd"), n_cond_q = 20, horizon = horizon, seed = seed + 27))

  # longrun: 40 quarters, both methods, own simulator from the latest data.
  long_h <- 40
  fc_long <- list(
    none = forecast_unconditional(d_none$beta, d_none$sigma, hist_none, lags, long_h, seed + 28),
    precovid_params = forecast_unconditional(d_pre$beta, d_pre$sigma, hist_pre_origin, lags, long_h, seed + 29))
  panels <- list(rgdp = TRUE, corepce = TRUE, urate = FALSE, gs2 = FALSE, gs10 = FALSE, rvol = FALSE)
  rows <- list()
  for (m in names(fc_long)) for (v in names(panels)) {
    q <- forecast_quantiles(fc_long[[m]], hist_none, v, growth = panels[[v]])
    for (h in c(4, 20, 40)) {
      rows[[length(rows) + 1]] <- data.frame(method = m, variable = v,
                                             units = if (panels[[v]]) "q/q % annualized" else "level",
                                             horizon = h, q05 = q[h, 1], median = q[h, 3], q95 = q[h, 5])
    }
  }
  res$longrun <- list(table = do.call(rbind, rows), horizon = long_h, panels = panels,
                      origin = rownames(hist_none)[nrow(hist_none)])
  saveRDS(res, file.path(out_dir, "forecast_checks.rds"))

  # Figures: fan charts, 40 quarters, one file per method.
  for (m in names(fc_long)) {
    png(file.path(out_dir, paste0("longrun_", m, ".png")), width = 2000, height = 1200, res = 200)
    op <- par(mfrow = c(2, 3), mar = c(4, 3, 2, 1))
    for (v in names(panels)) {
      q <- forecast_quantiles(fc_long[[m]], hist_none, v, growth = panels[[v]])
      j <- match(v, colnames(hist_none))
      recent <- hist_none[(nrow(hist_none) - 40):nrow(hist_none), j]
      if (panels[[v]]) recent <- 4 * diff(recent) else recent <- recent[-1]
      plot_fan(q, recent, paste0(v, if (panels[[v]]) " (q/q % annualized)" else " (level)"))
    }
    par(op)
    dev.off()
  }

  # Print a summary.
  cat("known VAR: max |z| of means", round(max(abs(res$known$table$z_mean)), 2),
      "; variance ratios", paste(round(range(res$known$table$var_ratio), 3), collapse = " to "), "\n")
  print(res$predict$spread, digits = 3, row.names = FALSE)
  print(res$predict$noise, digits = 3, row.names = FALSE)
  cat("newdata identical:", res$history$newdata_identical, "\n")
  print(res$history$precovid, digits = 4, row.names = FALSE)
  cat("precovid newdata: max |mean of noise| / sd:", round(res$history$precovid_noise_mean_max, 3), "\n")
  str(res$cond)
  print(res$longrun$table, digits = 3, row.names = FALSE)
}

# ---- Stage: report ----
if (stage == "report") {
  out_abs <- normalizePath(out_dir, winslash = "/")
  quarto::quarto_render("reports/forecast_checks.qmd", execute_params = list(output_dir = out_abs),
                        quiet = TRUE)
  target <- file.path(out_dir, "forecast_checks.html")
  if (file.exists(target)) file.remove(target)
  if (!file.rename("reports/forecast_checks.html", target)) stop("could not move the report to ", target)
  cat("Wrote", target, "\n")
}
