# check_estimation_diagnostics.R: residual diagnostics of both production
# runs, with a summary of their MCMC diagnostics (R/estimation_diagnostics.R).
#
# Run from the project root (so .Rprofile activates renv):
#   Rscript scripts/check_estimation_diagnostics.R run
#   Rscript scripts/check_estimation_diagnostics.R report
#
# Reads, for each full run in outputs/ (one per COVID method): covid.rds
# (estimation data), draws.rds, mcmc_summary.rds; and the companion roots of
# the hyperparameter stability check (outputs/checks_hyper_stability/), if
# it has been run. Estimates nothing.
#
# Writes outputs/checks_estimation_diagnostics/ (estimation_diagnostics.rds,
# PNG figures, estimation_diagnostics.html).

for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)

args  <- commandArgs(trailingOnly = TRUE)
stage <- args[1]
if (is.na(stage)) stop("give a stage: run or report")

cfg     <- yaml::read_yaml("config/run.yml")
series  <- read_series_config("config/series.csv")
out_dir <- "outputs/checks_estimation_diagnostics"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
covid_start <- "2020Q1"     # residual statistics are also shown without 2020 onward

#' Write the figures of one run from its saved results.
write_figures <- function(r) {
  m <- r$method
  png(file.path(out_dir, paste0("residual_grid_", m, ".png")), 1500, 900, res = 130)
  plot_residual_grid(r$z, paste("Standardized residuals, method", m))
  dev.off()
  png(file.path(out_dir, paste0("residual_correlation_", m, ".png")), 1100, 1100, res = 130)
  plot_residual_correlation(r$e_pre, paste("Residual correlations to 2019Q4, method", m))
  dev.off()
}

# ---- Stage: run ----
if (stage == "run") {
  runs <- list.files(cfg$output_root, pattern = "^[0-9]{4}m[0-9]{2}_.+_full_[0-9a-f]{8}$",
                     full.names = TRUE)
  res <- list(runs = list())
  for (d in runs) {
    meta <- jsonlite::read_json(file.path(d, "meta.json"))
    method <- meta$covid$method
    covid <- readRDS(file.path(d, "covid.rds"))
    draws <- readRDS(file.path(d, "draws.rds"))
    summ  <- readRDS(file.path(d, "mcmc_summary.rds"))
    y <- est_data_to_matrix(covid$est_data, series)
    rr <- var_residuals(y, draws)
    e_pre <- rows_before(rr$e, covid_start)
    z_pre <- rows_before(rr$z, covid_start)
    res$runs[[method]] <- list(
      method = method, run_id = basename(d),
      sample = paste(rownames(y)[1], "to", rownames(y)[nrow(y)]),
      residual_sample = paste(rr$quarters[1], "to", rr$quarters[length(rr$quarters)]),
      z = rr$z, e_pre = e_pre,
      outliers = residual_outliers(rr$z, 5),
      count_beyond_3 = c(all = sum(abs(rr$z) > 3), pre_2020 = sum(abs(z_pre) > 3)),
      cells = c(all = length(rr$z), pre_2020 = length(z_pre)),
      acf_all = residual_autocorrelation(rr$e, 8),
      acf_pre = if (nrow(e_pre) < nrow(rr$e)) residual_autocorrelation(e_pre, 8) else NULL,
      cor_top_all = top_residual_correlations(rr$e, 15),
      cor_top_pre = top_residual_correlations(e_pre, 15),
      cor_pairs = do.call(rbind, lapply(list(c("cpi", "corecpi"), c("pce", "corepce"),
                                             c("rgdp", "rpce"), c("gs5", "gs10"), c("aaa", "baa"),
                                             c("rgdp", "payems"), c("urate", "payems")),
        function(p) data.frame(pair = paste(p, collapse = " / "),
                               all = stats::cor(rr$e[, p[1]], rr$e[, p[2]]),
                               pre_2020 = stats::cor(e_pre[, p[1]], e_pre[, p[2]])))),
      mcmc = list(hyper = summ$hyper, chains = summ$chains, rhat_multivariate = summ$rhat_multivariate,
                  sigma_pd = summ$sigma_pd, max_modulus = summ$max_modulus,
                  max_modulus_summary = summ$max_modulus_summary))
    write_figures(res$runs[[method]])
    cat(method, ": residuals", res$runs[[method]]$residual_sample, "; |z| > 5:",
        nrow(res$runs[[method]]$outliers), "; |z| > 3:", res$runs[[method]]$count_beyond_3[["all"]],
        "of", res$runs[[method]]$cells[["all"]], "\n")
  }
  # Companion roots at the stability-check origins (method none, mode-based draws).
  hs <- list.files("outputs/checks_hyper_stability", pattern = "^origin_[0-9]{4}Q[1-4][.]rds$",
                   full.names = TRUE)
  res$roots_origins <- do.call(rbind, lapply(hs, function(f) {
    o <- readRDS(f)
    data.frame(origin = o$run$origin, median = median(o$max_root$alt), q95 = quantile(o$max_root$alt, 0.95),
               max = max(o$max_root$alt), share_above_1 = mean(o$max_root$alt > 1), row.names = NULL)
  }))
  res$time <- Sys.time()
  saveRDS(res, file.path(out_dir, "estimation_diagnostics.rds"))
}

# ---- Stage: report ----
if (stage == "report") {
  res <- readRDS(file.path(out_dir, "estimation_diagnostics.rds"))
  for (r in res$runs) write_figures(r)
  out_abs <- normalizePath(out_dir, winslash = "/")
  quarto::quarto_render("reports/estimation_diagnostics.qmd", execute_params = list(output_dir = out_abs),
                        quiet = TRUE)
  target <- file.path(out_dir, "estimation_diagnostics.html")
  if (file.exists(target)) file.remove(target)
  if (!file.rename("reports/estimation_diagnostics.html", target)) stop("could not move the report to ", target)
  cat("Wrote", target, "\n")
}
