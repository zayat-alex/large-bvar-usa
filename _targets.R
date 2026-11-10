# _targets.R: the pipeline definition.
#
# Run from the project root (so .Rprofile activates renv):
#   Rscript scripts/build_data.R          (stage 1: data)
#   Rscript scripts/estimate.R            (stages 2-3: COVID treatment, estimation)
#   Rscript scripts/run_engine.R          (stage 4: baseline, scenarios, IRFs)
#   Rscript scripts/export.R              (stage 5: the website's bundle)
#
# What targets does: each target below is one step with explicit inputs. On a
# run, targets rebuilds a step only if its code or one of its inputs changed
# since the last run, and caches the results in _targets/ (not in git).
# Config files and raw data files are tracked as files, so editing
# config/run.yml or adding a new raw file and pointing the config to it
# rebuilds exactly what depends on it. The pipeline never goes online.
#
# Stages: 1 data (assembly, validation, report); 2 COVID treatment; 3
# estimation; 4 engine; 5 export. Stages 2-4 are built once per COVID
# method in config/run.yml (covid$methods), so the models are estimated side
# by side. Stage 4 reads each run's saved draws (outputs/<run_id>/draws.rds)
# as a file, not the estimation targets, so building it never re-estimates.
# Stage 5 exports one run, the published one (export$run).

library(targets)

tar_option_set(packages = c("yaml", "readxl", "quarto"))
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

# The COVID methods are read here, when the pipeline is defined, because each
# method gets its own targets (named <step>_<method>; "covid_data_" rather
# than "covid_" because R/covid.R has functions named covid_<method>).
covid_methods <- unlist(yaml::read_yaml("config/run.yml")$covid$methods)
# Likewise each scenario in the config gets its own engine target, so editing
# one scenario file reruns only that scenario.
scenario_ids <- as.character(unlist(yaml::read_yaml("config/run.yml")$scenarios))
# A scenario built on a tilted base (scenarios.yml, `base`; docs/methodology.md
# section 9.4) reuses the base's weights, so its target reads the base's
# engine file: a rerun of the base reruns the scenarios built on it.
scenario_bases <- vapply(scenario_ids, function(id) {
  scenario_base(read_scenario_index("scenarios/scenarios.yml"), id)
}, "")
missing_base <- setdiff(stats::na.omit(scenario_bases), scenario_ids)
if (length(missing_base)) stop("config/run.yml must list the base scenario(s): ", paste(missing_base, collapse = ", "))
# The run whose bundle is exported (stage 5).
export_method <- yaml::read_yaml("config/run.yml")$export$run

#' Targets of stages 2-3 for one COVID method.
#'
#' Each step is its own target, so a change reruns only the steps it affects:
#' changing the chain lengths reruns the chains and what follows, but not the
#' mode search or the Hessian, which do not depend on them.
#'
#' @param method A COVID method (for example "none").
#' @return List of targets: covid_data_, setup_, mode_, proposal_, chains_,
#'   draws_, mcmc_summary_, run_id_, run_files_, mcmc_plots_ and
#'   mcmc_report_<method>.
estimation_targets <- function(method) {
  nm  <- function(step) paste0(step, "_", method)
  sym <- function(step) as.symbol(nm(step))
  list(
    # Stage 2: COVID treatment (R/covid.R): estimation sample and forecast history.
    tar_target_raw(nm("covid_data"), bquote(covid_for_method(model_data, cfg, .(method)))),

    # Stage 3: estimation (R/estimate_bvar.R, R/estimate_run.R).
    # covid_keep(): the equations used (NULL = all; method exclude_window
    # leaves out the window's quarters).
    tar_target_raw(nm("setup"), bquote(
      estimation_setup(est_data_to_matrix(.(sym("covid_data"))$est_data, series), series,
                       cfg$prior, cfg$model$lags, keep = covid_keep(.(sym("covid_data")))))),
    # Step 1: mode of the hyperparameter posterior.
    tar_target_raw(nm("mode"), bquote(
      timed_call(find_hyper_mode, .(sym("setup")), unlist(mode_search$start_factors),
                 mode_search$maxit))),
    # Step 2: proposal covariance from the Hessian at the mode.
    tar_target_raw(nm("proposal"), bquote(
      timed_call(hyper_proposal, .(sym("mode"))$theta_mode, .(sym("setup"))))),
    # Step 3: Metropolis chains.
    tar_target_raw(nm("chains"), bquote(
      run_chains(.(sym("setup")), .(sym("mode"))$theta_mode, .(sym("proposal"))$W, mcmc,
                 n_cores))),
    # Step 4: pooled draws of the hyperparameters, Sigma and B.
    tar_target_raw(nm("draws"), bquote(pool_and_draw(.(sym("chains")), .(sym("setup")), mcmc))),
    # Convergence summary; derived quantities are forecast from origin_data,
    # the history the engine uses.
    tar_target_raw(nm("mcmc_summary"), bquote(
      summarize_mcmc(.(sym("chains")), .(sym("draws")), .(sym("setup")),
                     history = est_data_to_matrix(.(sym("covid_data"))$origin_data, series),
                     convergence = convergence))),

    # Run identifier and the run's files in <output_root>/<run_id>/.
    tar_target_raw(nm("run_id"), bquote(make_run_id(cfg, .(method), mcmc, series))),
    tar_target_raw(nm("run_files"), bquote(
      write_run_files(.(sym("run_id")), cfg$output_root, .(method), .(sym("covid_data")),
                      .(sym("mode")), .(sym("proposal")), .(sym("chains")), .(sym("draws")),
                      .(sym("mcmc_summary")), mcmc, cfg)),
      format = "file"),

    # MCMC diagnostics: figures and report in the run's folder (R/mcmc_diagnostics.R).
    tar_target_raw(nm("mcmc_plots"), bquote(
      unname(save_mcmc_plots(.(sym("chains")), .(sym("mcmc_summary")), .(sym("setup")),
                             file.path(dirname(.(sym("run_files"))[1]), "plots")))),
      format = "file"),
    tar_target_raw(nm("mcmc_report"), bquote(
      render_mcmc_report(mcmc_report_source, dirname(.(sym("run_files"))[1]),
                         depends_on = c(.(sym("run_files")), .(sym("mcmc_plots"))))),
      format = "file")
  )
}

#' Targets of stage 4 (engine) for one COVID method.
#'
#' The engine reads the run's saved files, tracked as files: a new draws.rds
#' (after a real re-estimation) reruns the engine, but building the engine
#' does not depend on, and so never triggers, the estimation targets.
#'
#' @param method A COVID method.
#' @return List of targets: engine_run_dir_, engine_run_files_,
#'   engine_baseline_, engine_scenario_<method>_<id> for each scenario,
#'   engine_irf_ and engine_meta_<method>.
engine_targets <- function(method) {
  nm  <- function(step) paste0(step, "_", method)
  sym <- function(step) as.symbol(nm(step))
  scen_nm <- function(id) paste0("engine_scenario_", method, "_", id)
  scenario_targets <- lapply(scenario_ids, function(id) {
    b <- scenario_bases[[id]]
    base_file <- if (is.na(b)) NULL else as.symbol(scen_nm(b))     # the base's engine file
    tar_target_raw(scen_nm(id), bquote({
      .(sym("engine_run_files"))
      run_engine_scenario(.(id), .(as.symbol(paste0("scenario_file_", id))), .(sym("engine_baseline")),
                          .(sym("engine_run_dir")), series, engine_lags, engine_cfg,
                          .(as.symbol(paste0("scenario_soft_", id))), base_file = .(base_file))
    }), format = "file")
  })
  scenario_syms <- lapply(scenario_ids, function(id) as.symbol(scen_nm(id)))
  c(list(
    # The run's folder, from its identifier (cheap: the config, not the draws).
    tar_target_raw(nm("engine_run_dir"), bquote(file.path(cfg$output_root, .(sym("run_id"))$run_id))),
    tar_target_raw(nm("engine_run_files"), bquote(
      file.path(.(sym("engine_run_dir")), c("draws.rds", "covid.rds"))), format = "file"),
    # Baseline: unconditional forecast (paths and per-draw means).
    tar_target_raw(nm("engine_baseline"), bquote({
      .(sym("engine_run_files"))
      run_engine_baseline(.(sym("engine_run_dir")), series, engine_lags, engine_cfg)
    }), format = "file"),
    # Impulse responses of the configured identifications.
    tar_target_raw(nm("engine_irf"), bquote({
      .(sym("engine_run_files"))
      run_engine_irfs(.(sym("engine_run_dir")), series, engine_lags, engine_cfg)
    }), format = "file")),
    scenario_targets,
    list(
      # engine/meta.json: settings, inputs' MD5s, scenario checks, timings.
      tar_target_raw(nm("engine_meta"), bquote(
        write_engine_meta(.(sym("engine_run_dir")), .(sym("run_id")), engine_cfg,
                          c(.(sym("engine_baseline")), .(sym("engine_irf")), ..(scenario_syms))),
        splice = TRUE),
        format = "file")))
}

list(
  # ---- Configuration and raw inputs (tracked as files) ----
  tar_target(run_config_file, "config/run.yml", format = "file"),
  tar_target(series_config_file, "config/series.csv", format = "file"),
  tar_target(report_source, "reports/data_validation.qmd", format = "file"),
  tar_target(mcmc_report_source, "reports/mcmc_diagnostics.qmd", format = "file"),
  tar_target(cfg, yaml::read_yaml(run_config_file)),
  tar_target(series, read_series_config(series_config_file)),
  # Sections of the config, as separate targets: a step that uses only one
  # section is not rebuilt when another section changes (for example, the
  # data are not rebuilt when the MCMC settings change).
  tar_target(cfg_data, cfg[c("data", "sample")]),
  tar_target(mode_search, cfg$mcmc$mode_search),
  tar_target(mcmc, mcmc_settings(cfg)),
  tar_target(n_cores, cfg$mcmc$n_cores),
  # Convergence targets: they judge a run, they do not change it, so they are
  # kept out of the MCMC settings (and out of the run identifier).
  tar_target(convergence, cfg$mcmc$convergence),
  tar_target(raw_files, raw_file_paths(cfg_data), format = "file"),

  # ---- Stage 1: data (docs/data.md, sections 1-3) ----
  # raw_files is passed so that a changed raw file invalidates these steps.
  tar_target(series_list, { raw_files; build_series_list(series, cfg_data) }),
  tar_target(levels, assemble_levels(series_list, cfg_data$sample$start, cfg_data$sample$end)),
  tar_target(model_data, transform_to_model_units(levels, series)),
  tar_target(data_files,
             write_data_tables(levels, model_data, series, cfg_data, attr(series_list, "backfill"),
                               cfg_data$data$output_dir),
             format = "file"),

  # ---- Stage 1: validation (docs/data.md, section 4) ----
  tar_target(checks, { raw_files; run_data_checks(levels, model_data, series_list, series, cfg_data) }),
  tar_target(check_plots,
             unname(save_check_plots(checks, file.path(cfg_data$data$output_dir, "checks"),
                                     cfg_data$data$api_realtime_date)),
             format = "file"),
  tar_target(validation_file,
             write_validation_results(checks, levels, cfg_data, cfg_data$data$output_dir),
             format = "file"),
  tar_target(exploratory_pdf,
             plot_exploratory_pdf(model_data, series,
                                  file.path(cfg_data$data$output_dir, "exploratory_series.pdf"),
                                  sample_description(levels, cfg_data)),
             format = "file"),
  tar_target(validation_report,
             render_validation_report(report_source, cfg_data$data$output_dir,
                                      depends_on = c(validation_file, check_plots)),
             format = "file"),

  # ---- Stages 2-3: COVID treatment and estimation, one set per method ----
  lapply(covid_methods, estimation_targets),

  # ---- Stage 4: engine (R/engine_run.R), one set per method ----
  # Settings shared by the methods; the scenario index and files are
  # tracked as files.
  tar_target(engine_cfg, engine_settings(cfg, series$short_name)),
  tar_target(engine_lags, cfg$model$lags),
  tar_target(scenario_index_file, "scenarios/scenarios.yml", format = "file"),
  tar_target(scenario_index, read_scenario_index(scenario_index_file)),
  lapply(scenario_ids, function(id) {
    tar_target_raw(paste0("scenario_file_", id),
                   bquote(scenario_files(.(id), scenario_index, "scenarios")[[1]]), format = "file")
  }),
  # Each scenario's soft conditions (entropic tilting), as its own target:
  # editing one entry of scenarios.yml reruns only that scenario.
  lapply(scenario_ids, function(id) {
    tar_target_raw(paste0("scenario_soft_", id), bquote(scenario_soft(scenario_index, .(id))))
  }),
  lapply(covid_methods, engine_targets),

  # ---- Stage 5: export (R/export.R), the published run only ----
  # Generalized IRFs of every variable not in irf$girf_shocks, as quantile
  # tables only (R/engine_run.R, run_engine_irf_extra()): the website offers
  # them, and their draws are not kept. Built from the run's draws.rds.
  tar_target_raw("engine_irf_extra", bquote({
    .(as.symbol(paste0("engine_run_files_", export_method)))
    run_engine_irf_extra(.(as.symbol(paste0("engine_run_dir_", export_method))), series, engine_lags, engine_cfg,
                         export_settings(cfg, series)$spreads, unlist(cfg$quantiles))
  }), format = "file"),
  # The bundle the website reads, from the engine files of the run named in
  # config/run.yml (export$run); engine_meta_<method> is a file target, so a
  # rebuilt engine reruns the export.
  tar_target_raw("export_bundle", bquote(
    run_export(.(as.symbol(paste0("engine_run_dir_", export_method))),
               .(as.symbol(paste0("engine_meta_", export_method))), series, scenario_index, cfg,
               engine_irf_extra)),
    format = "file")
)
