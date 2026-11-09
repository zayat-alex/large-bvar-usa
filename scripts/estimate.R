# estimate.R: run the COVID treatment and estimation stages of the pipeline
# (_targets.R, stages 2-3) and print a summary of each run.
#
# Run from the project root (so .Rprofile activates renv):
#   Rscript scripts/estimate.R
#
# One run per COVID method in config/run.yml (covid$methods), with the chain
# lengths of the MCMC profile in mcmc$profile (smoke: a short check of the
# pipeline; full: the production run). The data stage is built first if it
# is out of date. Steps whose inputs have not changed are skipped (cache in
# _targets/): changing only the chain lengths reruns the chains but not the
# mode search or the Hessian.
#
# Time at full size (31 variables, 5 lags), per COVID method: mode search
# about 4-8 minutes, Hessian about 3 minutes, then about 35 ms per Metropolis
# iteration per chain (chains run in parallel on mcmc$n_cores processes) and
# about 20 ms per coefficient draw.
#
# Writes <output_root>/<run_id>/ for each run (R/estimate_run.R,
# write_run_files()): covid.rds, draws.rds, estimation.rds, mcmc_summary.rds,
# meta.json; the diagnostic figures in plots/ and the report
# mcmc_diagnostics.html (R/mcmc_diagnostics.R).

library(targets)
cfg <- yaml::read_yaml("config/run.yml")
methods <- unlist(cfg$covid$methods)

# Every mcmc_report_<method> target, and everything it depends on (the run
# files and figures included).
tar_make(names = starts_with("mcmc_report_"))

for (m in methods) {
  id <- tar_read_raw(paste0("run_id_", m))
  s  <- tar_read_raw(paste0("mcmc_summary_", m))
  mode <- tar_read_raw(paste0("mode_", m))
  cat("\n== Run", id$run_id, "==\n")
  cat("Mode search:", round(mode$seconds / 60, 1), "min; log p(theta) at the mode:",
      format(max(mode$table$log_post_theta), nsmall = 2), "; all starts converged:",
      all(mode$table$convergence == 0), "\n")
  print(s$chains, digits = 3, row.names = FALSE)
  h <- s$hyper
  cat("ESS (all chains) min", round(min(h$ess)), "median", round(median(h$ess)),
      "| R-hat max", signif(max(h$rhat), 3), "| max |Geweke z|", signif(max(h$geweke_max_abs_z), 3),
      "| draws near a bound (share, max)", signif(max(h$share_near_bound, na.rm = TRUE), 3), "\n")
  cat("Sigma positive definite in", s$sigma_pd * 100, "% of kept draws; largest root:",
      paste(names(s$max_modulus_summary), signif(s$max_modulus_summary, 4), collapse = ", "), "\n")
  print(h[h$hyperparameter %in% c("lambda", "soc"), c("hyperparameter", "q05", "median", "q95")],
        digits = 3, row.names = FALSE)
  if (!is.null(s$convergence)) {
    cat("Convergence targets:\n")
    for (i in seq_len(nrow(s$convergence))) {
      r <- s$convergence[i, ]
      cat(" ", if (r$pass) "pass" else "FLAG", "|", r$check, "| worst", r$worst, "| target", r$threshold, "\n")
    }
  }
  cat("Files and report (mcmc_diagnostics.html) in", file.path(cfg$output_root, id$run_id), "\n")
}
