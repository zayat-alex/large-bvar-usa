# engine_run.R: stage 4 of the pipeline, the forecasting engine. From each
# estimation run's saved posterior draws it computes the baseline forecast,
# the configured scenarios and their effects, and the impulse responses,
# and writes them to <output_root>/<run_id>/engine/.
#
# Why it exists: the estimation stage produces draws of the parameters
# (B, Sigma); everything a reader sees (fan charts, scenario effects, IRFs)
# is computed from those draws here. Crump et al. (2021), Section 2.3:
#   - baseline: the unconditional predictive density, simulated draw by draw
#     (docs/methodology.md, M18);
#   - scenario: the predictive density given fixed future values of some
#     variables (a conditional forecast; M21-M23 for the Gaussian method,
#     M24-M26 for the Kalman simulation smoother);
#   - scenario effect: for each draw theta, D = E[y | C, theta] - E[y | theta],
#     the conditional mean minus the unconditional mean. Both means are
#     exact (no simulated shocks), so the spread of D over draws is
#     parameter uncertainty only, as for the IRFs. The scenario's own fan
#     chart comes from the simulated paths, which add future shocks;
#   - impulse responses: closed-form Cholesky responses per draw (M27-M29,
#     R/irf.R);
#   - soft conditions: a scenario's soft targets are imposed by reweighting
#     its paths (entropic tilting, Section 4.2; M30-M33, R/tilting.R).
#
# Common random numbers: the baseline and every scenario use the same seed,
# so their simulated paths share the same underlying random numbers and
# differences between their fans are not sampling noise in the shocks.
#
# The engine reads the run's saved files (draws.rds, covid.rds) rather than
# the estimation targets of the pipeline, so running it never triggers a
# re-estimation; a new draws.rds (a real re-estimation) reruns it.
#
# Each piece is computed by a core function on in-memory inputs
# (engine_baseline(), engine_scenario(), engine_irfs()), which the tests
# and the stored reference outputs use, and run on a run folder by a
# wrapper that reads the files and writes the result (run_engine_*()).
#
# Inputs:  <run_dir>/draws.rds and covid.rds, config/run.yml (sections
#          forecast (including tilting), irf, quantiles, scenarios), config/series.csv,
#          scenarios/scenarios.yml and scenarios/<id>.csv.
# Outputs: <run_dir>/engine/baseline.rds, scenario_<id>.rds, irf.rds,
#          irf_extra.rds (published run only: generalized IRFs of the other
#          variables, quantiles only; run_engine_irf_extra()) and
#          meta.json. All values in model units (100 ln for log series,
#          natural units for rates).


# ---- Settings and inputs ----

#' Engine settings from the run configuration, checked.
#'
#' @param cfg Parsed config/run.yml.
#' @param variables Model variable names (series order).
#' @return List: horizon, method, seed, seed_offset, paths_per_draw,
#'   quantiles, irf (irf_settings()). Stops on an invalid setting. The
#'   scenario ids (cfg$scenarios) are not included: the pipeline gives each
#'   scenario its own target, so adding one does not rerun the others.
engine_settings <- function(cfg, variables) {
  f <- cfg$forecast
  if (is.null(f$horizon) || f$horizon < 1) stop("forecast$horizon must be a positive number of quarters")
  if (is.null(f$conditional_method) || !f$conditional_method %in% c("gaussian", "kalman")) {
    stop("forecast$conditional_method must be \"gaussian\" or \"kalman\"")
  }
  if (is.null(f$seed_offset) || is.null(cfg$mcmc$seed)) stop("forecast$seed_offset and mcmc$seed must be set")
  n_paths <- f$paths_per_draw
  if (is.null(n_paths) || n_paths < 1 || n_paths != round(n_paths)) {
    stop("forecast$paths_per_draw must be a positive whole number")
  }
  probs <- unlist(cfg$quantiles)
  if (is.null(probs) || any(probs <= 0 | probs >= 1)) stop("quantiles must lie strictly between 0 and 1")
  list(horizon = f$horizon, method = f$conditional_method,
       seed = cfg$mcmc$seed + f$seed_offset, seed_offset = f$seed_offset,
       paths_per_draw = as.integer(n_paths), quantiles = probs,
       tilting = tilting_settings(f$tilting),
       irf = irf_settings(cfg, variables))
}

#' Tilting settings (config forecast$tilting), checked.
#'
#' @param t The `tilting` section: paths_per_draw (paths per posterior
#'   draw for scenarios with soft conditions), tol and maxit (solver,
#'   tilt_weights()), ess_warn_share (warn when the effective sample size is
#'   below this share of the paths).
#' @return The section with paths_per_draw and maxit as integers.
tilting_settings <- function(t) {
  if (is.null(t)) stop("forecast$tilting must be set")
  ok <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x)
  if (!ok(t$paths_per_draw) || t$paths_per_draw < 1 || t$paths_per_draw != round(t$paths_per_draw)) {
    stop("forecast$tilting$paths_per_draw must be a positive whole number")
  }
  if (!ok(t$tol) || t$tol <= 0) stop("forecast$tilting$tol must be positive")
  if (!ok(t$maxit) || t$maxit < 1) stop("forecast$tilting$maxit must be at least 1")
  if (!ok(t$ess_warn_share) || t$ess_warn_share <= 0 || t$ess_warn_share >= 1) {
    stop("forecast$tilting$ess_warn_share must lie strictly between 0 and 1")
  }
  list(paths_per_draw = as.integer(t$paths_per_draw), tol = t$tol, maxit = as.integer(t$maxit),
       ess_warn_share = t$ess_warn_share)
}

#' Check that every configured scenario is in the index and has a file.
#'
#' @param ids Scenario ids from the config (cfg$scenarios).
#' @param index Output of read_scenario_index().
#' @param dir Folder of the scenario files.
#' @return The file paths, named by id; stops otherwise.
scenario_files <- function(ids, index, dir) {
  missing_id <- setdiff(ids, index$id)
  if (length(missing_id)) stop("scenario(s) not in scenarios.yml: ", paste(missing_id, collapse = ", "))
  paths <- stats::setNames(file.path(dir, paste0(ids, ".csv")), ids)
  absent <- paths[!file.exists(paths)]
  if (length(absent)) stop("scenario file(s) not found: ", paste(absent, collapse = ", "))
  paths
}

#' Read the inputs of the engine from an estimation run's folder.
#'
#' The forecast starts from origin_data (the latest data, whatever the
#' estimation sample; docs/methodology.md section 11).
#'
#' @param run_dir Folder of the run (draws.rds, covid.rds).
#' @param series Series configuration (read_series_config()).
#' @param lags Lags p of the configured model.
#' @param horizon Forecast horizon H.
#' @return List: beta, sigma, history (matrix, row names = quarters),
#'   lags, variables, last_observed, quarters (labels of T+1..T+H),
#'   draws_md5, covid_md5.
engine_inputs <- function(run_dir, series, lags, horizon) {
  files <- file.path(run_dir, c("draws.rds", "covid.rds"))
  absent <- files[!file.exists(files)]
  if (length(absent)) stop("run files not found (estimate the run first): ", paste(absent, collapse = ", "))
  draws <- readRDS(files[1])
  covid <- readRDS(files[2])
  history <- est_data_to_matrix(covid$origin_data, series)
  if (!identical(draws$variables, colnames(history))) stop("variable order of draws and history differ")
  if (!identical(as.integer(draws$lags), as.integer(lags))) {
    stop("the draws have ", draws$lags, " lags; the config has ", lags)
  }
  last_observed <- rownames(history)[nrow(history)]
  list(beta = draws$beta, sigma = draws$sigma, history = history, lags = as.integer(lags),
       variables = draws$variables, last_observed = last_observed,
       quarters = forecast_quarters(last_observed, horizon),
       draws_md5 = unname(tools::md5sum(files[1])), covid_md5 = unname(tools::md5sum(files[2])))
}


# ---- Quantiles ----

#' Quantiles over draws of a forecast array, per quarter and variable.
#'
#' Used for every forecast table of the engine. With tilting weights the
#' quantiles are weighted (weighted_quantiles(), R/tilting.R), which equal
#' the unweighted ones when the weights are equal.
#'
#' @param x Array draws x H x n.
#' @param probs Probabilities (config: quantiles).
#' @param quarters Labels of the H forecast quarters.
#' @param weights NULL, or one weight per row of x (tilting).
#' @return data.frame: quarter, variable, one column per probability (q05,
#'   q16, ...), in the units of x.
engine_quantiles <- function(x, probs, quarters, weights = NULL) {
  variables <- dimnames(x)[[3]]
  rows <- list()
  for (j in seq_along(variables)) {
    q <- t(apply(x[, , j, drop = FALSE], 2, weighted_quantiles, weights, probs))   # H x probs
    colnames(q) <- sprintf("q%02d", round(100 * probs))
    rows[[j]] <- data.frame(quarter = quarters, variable = variables[j], q)
  }
  do.call(rbind, rows)
}


# ---- Core computations ----

#' Baseline: the unconditional forecast.
#'
#' Computed with the conditional forecast function and no fixed cell, so the
#' baseline and the scenarios go through the same code: an empty scenario
#' then has an effect of exactly zero. With no fixed cell this is the
#' unconditional predictive density (M18); with one path per draw it uses
#' the same random numbers as forecast_unconditional().
#'
#' @param inputs Output of engine_inputs() (or a list with the same fields).
#' @param settings Output of engine_settings().
#' @return List: paths (draws x paths_per_draw rows, H x n; with shocks),
#'   mean (draws x H x n; each draw's forecast with future shocks at zero,
#'   E[y | theta]), quantiles (of the paths), seconds.
engine_baseline <- function(inputs, settings) {
  t0 <- proc.time()[["elapsed"]]
  free <- matrix(NA_real_, settings$horizon, length(inputs$variables),
                 dimnames = list(inputs$quarters, inputs$variables))
  out <- engine_forecast(free, inputs, settings)
  out$seconds <- proc.time()[["elapsed"]] - t0
  out
}

#' Scenario: the conditional forecast given the scenario's fixed cells,
#' tilted to its soft targets if it has any.
#'
#' Hard cells first: paths are drawn given the fixed cells (none for a
#' scenario with soft conditions only). Then, if there are soft targets,
#' the paths are reweighted by entropic tilting (R/tilting.R). Tilting
#' changes only weights, so the hard cells stay exactly met. A tilted
#' scenario uses forecast$tilting$paths_per_draw paths per draw, since the
#' weights fall on paths and need enough of them.
#'
#' Effects: `effect` is, per draw, the conditional mean given the hard
#' cells minus the unconditional mean (the paper's D; zero for a scenario
#' without hard cells). For a tilted scenario the effect of the whole
#' scenario is one path, tilt$effect: the weighted mean of the tilted paths
#' minus the baseline mean averaged over draws. The weights belong to paths,
#' not to draws, so it has no distribution over draws.
#'
#' Built on a base (a tilted forecast, docs/methodology.md section 9.4,
#' "Scenarios built on a tilted forecast"): the paths are drawn given the
#' scenario's fixed cells as usual, and its quantiles are weighted with the
#' base's tilting weights. This is valid because every scenario of a run
#' uses the same seed, so path j of the scenario and path j of the base come
#' from the same random shocks; the base reweights those shock histories,
#' and the scenario runs the reweighted histories through its conditions.
#' It needs the same number of paths per draw (checked).
#'
#' @param conditions H x n matrix in model units, NA = free (read_scenario()).
#' @param inputs Output of engine_inputs().
#' @param settings Output of engine_settings().
#' @param baseline_mean The baseline's `mean` (draws x H x n).
#' @param targets NULL or a table of soft targets (soft_targets()).
#' @param base NULL, or list(id, weights, ess, ess_share) of the base
#'   (run_engine_scenario()); not together with targets.
#' @return List: conditions, paths, mean (E[y | C, theta] per draw),
#'   effect (mean - baseline_mean per draw), quantiles (of the paths,
#'   weighted if tilted), effect_quantiles, fixed_cells, fixed_max_error
#'   (largest |path - condition| over fixed cells and paths), seconds; if
#'   tilted also quantiles_untilted and tilt (targets, weights, gamma,
#'   iterations, ess, ess_share, ess_warning, kl, check, effect); if built
#'   on a base also quantiles_unweighted and base (id, weights, ess,
#'   ess_share).
engine_scenario <- function(conditions, inputs, settings, baseline_mean, targets = NULL, base = NULL) {
  t0 <- proc.time()[["elapsed"]]
  tilted <- !is.null(targets) && nrow(targets) > 0
  if (tilted && !is.null(base)) stop("a scenario with a base cannot have soft targets of its own")
  # A based scenario draws as many paths as its (tilted) base, so path j of
  # both comes from the same random shocks.
  n_paths <- if (tilted || !is.null(base)) settings$tilting$paths_per_draw else settings$paths_per_draw
  out <- engine_forecast(conditions, inputs, settings, n_paths)
  if (!is.null(base)) {
    if (length(base$weights) != dim(out$paths)[1]) stop("base ", base$id, " has another number of paths")
    out$quantiles_unweighted <- out$quantiles
    out$quantiles <- engine_quantiles(out$paths, settings$quantiles, inputs$quarters, base$weights)
    out$base <- base
  }
  if (!identical(dim(baseline_mean), dim(out$mean))) stop("baseline mean does not match the scenario's draws")
  out$conditions <- conditions
  out$effect <- out$mean - baseline_mean
  out$effect_quantiles <- engine_quantiles(out$effect, settings$quantiles, inputs$quarters)
  if (tilted) {
    s <- settings$tilting
    tl <- tilt_paths(out$paths, targets, inputs$history, s$tol, s$maxit)
    J <- length(tl$weights)
    out$quantiles_untilted <- out$quantiles
    out$quantiles <- engine_quantiles(out$paths, settings$quantiles, inputs$quarters, tl$weights)
    # Weighted mean of the tilted paths minus the baseline mean (H x n).
    tilted_mean <- apply(out$paths, c(2, 3), function(v) sum(tl$weights * v))
    tl$effect <- tilted_mean - apply(baseline_mean, c(2, 3), mean)
    tl$targets <- targets
    tl$ess_share <- tl$ess / J
    tl$ess_warning <- tl$ess_share < s$ess_warn_share
    if (tl$ess_warning) {
      warning("tilting: effective sample size ", round(tl$ess), " of ", J, " paths (",
              signif(100 * tl$ess_share, 2), "%), below the warning share of ",
              100 * s$ess_warn_share, "%", call. = FALSE)
    }
    out$tilt <- tl
  }
  # Fixed cells must hold in every path (spec check; exact up to rounding).
  fixed <- which(!is.na(conditions), arr.ind = TRUE)
  out$fixed_cells <- nrow(fixed)
  out$fixed_max_error <- if (nrow(fixed) == 0) 0 else max(vapply(seq_len(nrow(fixed)), function(k) {
    h <- fixed[k, 1]; j <- fixed[k, 2]
    max(abs(out$paths[, h, j] - conditions[h, j]), abs(out$mean[, h, j] - conditions[h, j]))
  }, numeric(1)))
  out$seconds <- proc.time()[["elapsed"]] - t0
  out
}

#' Simulated paths and per-draw means for one conditions matrix.
#'
#' @param conditions H x n matrix, NA = free.
#' @param inputs,settings As in engine_scenario().
#' @param n_paths Paths per posterior draw.
#' @return List: paths ((draws x n_paths) x H x n, paths of draw 1 first),
#'   mean (draws x H x n), quantiles.
engine_forecast <- function(conditions, inputs, settings, n_paths = settings$paths_per_draw) {
  args <- list(beta = inputs$beta, sigma = inputs$sigma, history = inputs$history, lags = inputs$lags,
               horizon = settings$horizon, conditions = conditions, seed = settings$seed,
               method = settings$method)
  paths <- do.call(forecast_conditional, c(args, list(shocks = TRUE, n_paths = n_paths)))
  mean  <- do.call(forecast_conditional, c(args, list(shocks = FALSE)))
  list(paths = paths, mean = mean,
       quantiles = engine_quantiles(paths, settings$quantiles, inputs$quarters))
}

#' Impulse responses of the configured identifications.
#'
#' @param inputs Output of engine_inputs().
#' @param settings Output of engine_settings().
#' @return List: draws (named list of arrays draws x (H + 1) x n, from
#'   irf_draws()), quantiles (named list of tables, irf_quantiles()),
#'   horizon, seconds.
engine_irfs <- function(inputs, settings) {
  t0 <- proc.time()[["elapsed"]]
  s <- settings$irf
  draws <- irf_draws(inputs$beta, inputs$sigma, inputs$lags, s$horizon, inputs$variables,
                     s$girf_shocks, s$recursive)
  list(draws = draws, quantiles = lapply(draws, irf_quantiles, settings$quantiles),
       horizon = s$horizon, seconds = proc.time()[["elapsed"]] - t0)
}


# ---- Wrappers for the pipeline: read a run, compute, write ----

#' Folder of the engine outputs of a run (created if needed).
engine_dir <- function(run_dir) {
  d <- file.path(run_dir, "engine")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

#' Baseline of a run, written to <run_dir>/engine/baseline.rds.
#'
#' @param run_dir Folder of the estimation run.
#' @param series,lags,settings Series configuration, lags p, engine settings.
#' @return Path of the file written.
run_engine_baseline <- function(run_dir, series, lags, settings) {
  inputs <- engine_inputs(run_dir, series, lags, settings$horizon)
  out <- engine_baseline(inputs, settings)
  out$inputs <- inputs[c("last_observed", "quarters", "draws_md5", "covid_md5")]
  path <- file.path(engine_dir(run_dir), "baseline.rds")
  saveRDS(out, path)
  path
}

#' One scenario of a run, written to <run_dir>/engine/scenario_<id>.rds.
#'
#' @param id Scenario id.
#' @param scenario_file Path to scenarios/<id>.csv.
#' @param baseline_file Path written by run_engine_baseline() for this run.
#' @param run_dir,series,lags,settings As in run_engine_baseline().
#' @param soft The scenario's soft conditions (scenario_soft()); NULL or zero
#'   rows for a scenario with hard conditions only.
#' @param base_file NULL, or the file written by this function for the
#'   scenario's base (scenario_base()): a tilted scenario that fixes no cell,
#'   computed from the same draws with the same number of paths. Its tilting
#'   weights are reused (docs/methodology.md, section 9.4, "Scenarios built
#'   on a tilted forecast").
#' @return Path of the file written.
run_engine_scenario <- function(id, scenario_file, baseline_file, run_dir, series, lags, settings,
                                soft = NULL, base_file = NULL) {
  inputs <- engine_inputs(run_dir, series, lags, settings$horizon)
  baseline <- readRDS(baseline_file)
  if (!identical(baseline$inputs$draws_md5, inputs$draws_md5)) {
    stop("the baseline was computed from other draws; rerun the baseline first")
  }
  conditions <- read_scenario(scenario_file, inputs$variables, inputs$last_observed, settings$horizon, soft)
  targets <- soft_targets(scenario_file, inputs$variables, inputs$last_observed, settings$horizon, soft, series)
  base <- NULL
  if (!is.null(base_file)) {
    b <- readRDS(base_file)
    if (!identical(b$draws_md5, inputs$draws_md5)) stop("base ", b$id, " was computed from other draws; rerun it first")
    if (is.null(b$tilt)) stop("base ", b$id, " is not tilted")
    if (b$fixed_cells > 0) stop("base ", b$id, " fixes cells; a base may have soft conditions only")
    base <- list(id = b$id, weights = b$tilt$weights, ess = b$tilt$ess, ess_share = b$tilt$ess_share)
  }
  out <- engine_scenario(conditions, inputs, settings, baseline$mean, targets, base)
  out$id <- id
  out$draws_md5 <- inputs$draws_md5
  out$file <- scenario_file
  out$file_md5 <- unname(tools::md5sum(scenario_file))
  path <- file.path(engine_dir(run_dir), paste0("scenario_", id, ".rds"))
  saveRDS(out, path)
  path
}

#' Impulse responses of a run, written to <run_dir>/engine/irf.rds.
#'
#' @param run_dir,series,lags,settings As in run_engine_baseline().
#' @return Path of the file written.
run_engine_irfs <- function(run_dir, series, lags, settings) {
  inputs <- engine_inputs(run_dir, series, lags, settings$horizon)
  out <- engine_irfs(inputs, settings)
  path <- file.path(engine_dir(run_dir), "irf.rds")
  saveRDS(out, path)
  path
}

#' Generalized IRFs of the other variables, as quantile tables only.
#'
#' irf.rds keeps every draw of the configured shocks (the check reports read
#' them). The website also offers the generalized response to a shock in
#' every other variable; keeping their draws would cost about 10 MB per
#' shock, and only quantiles are shown. So each is computed draw by draw,
#' summarized at once with export_irf() (R/export.R), which takes the
#' 4-quarter changes and the spreads per draw before the quantiles, and its
#' draws are dropped. One shock at a time keeps memory small.
#'
#' @param inputs Output of engine_inputs() (beta, sigma, lags, variables).
#' @param horizon Last horizon H.
#' @param skip Shocks already in irf.rds (config irf$girf_shocks).
#' @param vars Variables of the bundle (export_variables()).
#' @param probs Probabilities (config: quantiles).
#' @return The table of export_irf() for the shocks not in `skip`.
irf_extra_quantiles <- function(inputs, horizon, skip, vars, probs) {
  shocks <- setdiff(inputs$variables, skip)
  tabs <- lapply(shocks, function(v) {
    d <- irf_draws(inputs$beta, inputs$sigma, inputs$lags, horizon, inputs$variables, v, NULL)
    export_irf(d, vars, probs)
  })
  out <- do.call(rbind, tabs)
  rownames(out) <- NULL
  out
}

#' The extra generalized IRFs of a run, written to <run_dir>/engine/irf_extra.rds.
#'
#' Run for the published run only (stage 5 of the pipeline, before the
#' export), since only the website uses them.
#'
#' @param run_dir,series,lags,settings As in run_engine_baseline().
#' @param spreads The spreads of the export (export_settings()), whose
#'   responses are part of the tables.
#' @param probs Probabilities (config: quantiles).
#' @return Path of the file written. The file holds the table, the MD5 of
#'   the draws it came from, the horizon and the timing.
run_engine_irf_extra <- function(run_dir, series, lags, settings, spreads, probs) {
  t0 <- proc.time()[["elapsed"]]
  inputs <- engine_inputs(run_dir, series, lags, settings$horizon)
  vars <- export_variables(series, spreads)
  s <- settings$irf
  out <- list(table = irf_extra_quantiles(inputs, s$horizon, s$girf_shocks, vars, probs),
              draws_md5 = inputs$draws_md5, horizon = s$horizon,
              seconds = proc.time()[["elapsed"]] - t0)
  path <- file.path(engine_dir(run_dir), "irf_extra.rds")
  saveRDS(out, path)
  path
}

#' Write <run_dir>/engine/meta.json: what the engine ran and on what.
#'
#' The engine settings are not in the run_id (they do not change the
#' estimation), so they are recorded here, with the MD5 of the draws and of
#' every scenario file used, the code version and the timings.
#'
#' @param run_dir Folder of the estimation run.
#' @param run_id Output of make_run_id() for the run.
#' @param settings Output of engine_settings().
#' @param files Paths written by the run_engine_*() functions of this run.
#' @return Path of meta.json.
write_engine_meta <- function(run_dir, run_id, settings, files) {
  baseline <- readRDS(files[basename(files) == "baseline.rds"])
  scen <- lapply(files[grepl("^scenario_.*\\.rds$", basename(files))], function(f) {
    s <- readRDS(f)
    tilt <- if (is.null(s$tilt)) NULL else
      list(targets = nrow(s$tilt$targets), paths = length(s$tilt$weights),
           ess = s$tilt$ess, ess_share = s$tilt$ess_share, ess_warning = s$tilt$ess_warning,
           iterations = s$tilt$iterations, kl = s$tilt$kl,
           largest_miss = max(abs(s$tilt$check$achieved - s$tilt$check$goal)))
    base <- if (is.null(s$base)) NULL else list(id = s$base$id, ess_share = s$base$ess_share)
    list(id = s$id, file = s$file, file_md5 = s$file_md5, fixed_cells = s$fixed_cells,
         fixed_max_error = s$fixed_max_error, tilting = tilt, base = base, seconds = s$seconds)
  })
  irf <- readRDS(files[basename(files) == "irf.rds"])
  meta <- list(
    run_id = run_id$run_id, created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    origin_last_quarter = baseline$inputs$last_observed,
    forecast_quarters = baseline$inputs$quarters[c(1, length(baseline$inputs$quarters))],
    draws_md5 = baseline$inputs$draws_md5, covid_md5 = baseline$inputs$covid_md5,
    n_draws = dim(baseline$mean)[1],
    settings = settings[c("horizon", "method", "seed", "seed_offset", "paths_per_draw", "quantiles", "tilting")],
    scenario_effect = paste("per draw: conditional mean given the hard cells minus unconditional mean",
                            "(Crump et al. 2021, Section 2.3); tilted scenarios: weighted mean of the",
                            "tilted paths minus the baseline mean averaged over draws"),
    scenarios = unname(scen),
    irf = list(horizon = irf$horizon, identifications = names(irf$draws),
               recursive = settings$irf$recursive),
    timings_seconds = list(baseline = baseline$seconds, irf = irf$seconds),
    files = as.list(stats::setNames(unname(tools::md5sum(files)), basename(files))),
    code = git_state(), r_version = R.version.string, rng_kind = RNGkind())
  path <- file.path(engine_dir(run_dir), "meta.json")
  jsonlite::write_json(meta, path, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
  path
}
