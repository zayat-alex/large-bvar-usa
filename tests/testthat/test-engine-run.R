# test-engine-run.R: the engine stage (R/engine_run.R): stored reference
# outputs, the definition of a scenario's effect, the baseline, settings, and
# the wrappers that read a run folder and write the engine's files.
#
# The known VAR(2) and its reference outputs are stored in
# fixtures/engine_reference.rds (written by scripts/make_engine_reference.R).

ref <- readRDS(test_path("fixtures", "engine_reference.rds"))
drop_time <- function(x) { x$seconds <- NULL; x }

test_that("the engine reproduces the stored reference outputs", {
  baseline <- drop_time(engine_baseline(ref$inputs, ref$settings))
  scenario <- drop_time(engine_scenario(ref$conditions, ref$inputs, ref$settings, baseline$mean))
  irf <- drop_time(engine_irfs(ref$inputs, ref$settings))
  expect_equal(baseline, ref$baseline, tolerance = 1e-8)
  expect_equal(scenario, ref$scenario, tolerance = 1e-8)
  expect_equal(irf, ref$irf, tolerance = 1e-8)
  # The Kalman method gives the same paths with the same seed (docs/methodology.md, section 9.2).
  kal <- ref$settings; kal$method <- "kalman"
  scenario_k <- engine_scenario(ref$conditions, ref$inputs, kal, baseline$mean)
  expect_equal(scenario_k$paths, ref$scenario$paths, tolerance = 1e-8)
  expect_equal(scenario_k$effect, ref$scenario$effect, tolerance = 1e-8)
})

test_that("a tilted scenario reproduces the reference and keeps its hard cells", {
  tilted <- drop_time(engine_scenario(ref$conditions, ref$inputs, ref$settings, ref$baseline$mean, ref$targets))
  expect_equal(tilted, ref$tilted, tolerance = 1e-8)
  # Tilting changes weights, not paths: the paths are the hard-only paths
  # with tilting's paths per draw, and the hard cells hold in every path.
  st <- ref$settings; st$paths_per_draw <- st$tilting$paths_per_draw
  hard <- engine_scenario(ref$conditions, ref$inputs, st, ref$baseline$mean)
  expect_identical(tilted$paths, hard$paths)
  expect_identical(tilted$quantiles_untilted, hard$quantiles)
  expect_lt(tilted$fixed_max_error, 1e-10)
  # Targets met; the tilted effect is the weighted mean minus the baseline mean.
  expect_lt(max(abs(tilted$tilt$check$achieved - tilted$tilt$check$goal)), 1e-10)
  w <- tilted$tilt$weights
  eff <- apply(tilted$paths, c(2, 3), function(v) sum(w * v)) - apply(ref$baseline$mean, c(2, 3), mean)
  expect_equal(tilted$tilt$effect, eff)
  expect_equal(tilted$tilt$ess_share, tilted$tilt$ess / length(w))
  # A low effective sample size is a warning, not an error.
  st <- ref$settings; st$tilting$ess_warn_share <- 0.99
  expect_warning(engine_scenario(ref$conditions, ref$inputs, st, ref$baseline$mean, ref$targets),
                 "effective sample size")
})

test_that("a scenario built on a tilted base reuses the base's weights on its own paths", {
  # The base: the unconditional paths tilted to the targets (no fixed cell).
  free <- ref$conditions; free[] <- NA
  base_run <- engine_scenario(free, ref$inputs, ref$settings, ref$baseline$mean, ref$targets)
  base <- list(id = "b", weights = base_run$tilt$weights, ess = base_run$tilt$ess,
               ess_share = base_run$tilt$ess_share)
  # No condition of its own: the scenario is the base exactly.
  same <- engine_scenario(free, ref$inputs, ref$settings, ref$baseline$mean, base = base)
  expect_identical(same$paths, base_run$paths)
  expect_identical(same$quantiles, base_run$quantiles)
  # With fixed cells: the paths are the hard scenario's (with the base's
  # number of paths, hence the same random shocks path by path), weighted
  # with the base's weights; the fixed cells hold in every path.
  based <- engine_scenario(ref$conditions, ref$inputs, ref$settings, ref$baseline$mean, base = base)
  st <- ref$settings; st$paths_per_draw <- st$tilting$paths_per_draw
  hard <- engine_scenario(ref$conditions, ref$inputs, st, ref$baseline$mean)
  expect_identical(based$paths, hard$paths)
  expect_identical(based$quantiles_unweighted, hard$quantiles)
  expect_identical(based$quantiles, engine_quantiles(hard$paths, st$quantiles, ref$inputs$quarters, base$weights))
  expect_lt(based$fixed_max_error, 1e-10)
  expect_equal(based$base$id, "b")
  # Not together with targets of its own; the weights must match the paths.
  expect_error(engine_scenario(ref$conditions, ref$inputs, ref$settings, ref$baseline$mean, ref$targets, base),
               "cannot have soft targets")
  short <- base; short$weights <- short$weights[-1]
  expect_error(engine_scenario(ref$conditions, ref$inputs, ref$settings, ref$baseline$mean, base = short),
               "another number of paths")
})

test_that("the baseline is the unconditional forecast", {
  st <- ref$settings; st$paths_per_draw <- 1L
  b <- engine_baseline(ref$inputs, st)
  # Same random numbers as forecast_unconditional() with the same seed.
  fu <- forecast_unconditional(ref$inputs$beta, ref$inputs$sigma, ref$inputs$history, 2, st$horizon, st$seed)
  expect_equal(b$paths, fu, tolerance = 1e-10)
  fu0 <- forecast_unconditional(ref$inputs$beta, ref$inputs$sigma, ref$inputs$history, 2, st$horizon,
                                st$seed, shocks = FALSE)
  expect_equal(b$mean, fu0, tolerance = 1e-10)
  # Paths per draw: rows are draws x paths, means one per draw.
  expect_equal(dim(ref$baseline$paths), c(3 * 2, 6, 3))
  expect_equal(dim(ref$baseline$mean), c(3, 6, 3))
})

test_that("a scenario's effect is the conditional mean minus the unconditional mean, per draw", {
  s <- ref$scenario
  for (i in 1:3) {
    exact <- conditional_moments_gaussian(ref$inputs$beta[i, , ], ref$inputs$sigma[i, , ], ref$inputs$history,
                                          2, 6, ref$conditions)$mean
    expect_equal(unname(s$mean[i, , ]), unname(exact), tolerance = 1e-10)
    expect_equal(s$effect[i, , ], s$mean[i, , ] - ref$baseline$mean[i, , ], tolerance = 1e-12)
  }
  # Fixed cells hold in every path and in every mean.
  expect_lt(s$fixed_max_error, 1e-10)
  expect_equal(s$fixed_cells, 3L)
  # No condition: the effect is exactly zero and the paths are the baseline's.
  free <- ref$conditions; free[] <- NA
  e0 <- engine_scenario(free, ref$inputs, ref$settings, ref$baseline$mean)
  expect_identical(max(abs(e0$effect)), 0)
  expect_equal(e0$paths, ref$baseline$paths, tolerance = 1e-12)
  # Null scenario: the same cells fixed at the baseline mean (one draw) give no effect.
  one <- ref$inputs
  one$beta <- one$beta[1, , , drop = FALSE]; one$sigma <- one$sigma[1, , , drop = FALSE]
  null <- ref$conditions; null[!is.na(null)] <- ref$baseline$mean[1, , ][!is.na(null)]
  en <- engine_scenario(null, one, ref$settings, ref$baseline$mean[1, , , drop = FALSE])
  expect_lt(max(abs(en$effect)), 1e-10)
  # Quantile tables have one row per quarter and variable.
  expect_equal(nrow(s$quantiles), 6 * 3)
  expect_equal(nrow(s$effect_quantiles), 6 * 3)
})

test_that("engine quantiles are the quantiles over draws", {
  x <- ref$scenario$paths
  q <- engine_quantiles(x, c(0.16, 0.5), ref$inputs$quarters)
  expect_equal(names(q), c("quarter", "variable", "q16", "q50"))
  row <- q[q$quarter == "2026Q4" & q$variable == "b", ]
  expect_equal(c(row$q16, row$q50), unname(quantile(x[, 2, "b"], c(0.16, 0.5))))
})

test_that("engine settings and scenario lists are checked", {
  root <- project_root()
  cfg <- yaml::read_yaml(file.path(root, "config", "run.yml"))
  series <- read_series_config(file.path(root, "config", "series.csv"))
  st <- engine_settings(cfg, series$short_name)
  expect_equal(st$seed, cfg$mcmc$seed + cfg$forecast$seed_offset)
  bad <- cfg; bad$forecast$conditional_method <- "other"
  expect_error(engine_settings(bad, series$short_name), "gaussian")
  bad <- cfg; bad$forecast$paths_per_draw <- 0
  expect_error(engine_settings(bad, series$short_name), "paths_per_draw")
  expect_equal(st$tilting$paths_per_draw, cfg$forecast$tilting$paths_per_draw)
  bad <- cfg; bad$forecast$tilting$ess_warn_share <- 1.5
  expect_error(engine_settings(bad, series$short_name), "ess_warn_share")
  bad <- cfg; bad$forecast$tilting <- NULL
  expect_error(engine_settings(bad, series$short_name), "tilting must be set")
  idx <- read_scenario_index(file.path(root, "scenarios", "scenarios.yml"))
  expect_true(all(file.exists(scenario_files(unlist(cfg$scenarios), idx, file.path(root, "scenarios")))))
  expect_error(scenario_files("no_such_scenario", idx, file.path(root, "scenarios")), "not in scenarios.yml")
})

test_that("the run-folder wrappers read a run and write the engine files", {
  run_dir <- file.path(tempfile("run_"))
  dir.create(run_dir)
  on.exit(unlink(run_dir, recursive = TRUE), add = TRUE)
  series <- data.frame(short_name = ref$inputs$variables, stringsAsFactors = FALSE)
  quarters <- as.Date(c("2025-04-01", "2025-07-01", "2025-10-01", "2026-01-01", "2026-04-01"))
  origin <- data.frame(quarter = quarters, ref$inputs$history, row.names = NULL)
  saveRDS(list(beta = ref$inputs$beta, sigma = ref$inputs$sigma, variables = ref$inputs$variables, lags = 2L),
          file.path(run_dir, "draws.rds"))
  saveRDS(list(est_data = origin, origin_data = origin), file.path(run_dir, "covid.rds"))
  scen <- tempfile(fileext = ".csv")
  on.exit(unlink(scen), add = TRUE)
  writeLines(c("quarter,c,a,b", "2026Q3,,2,", "2026Q4,1,,0"), scen)

  st <- ref$settings
  expect_error(run_engine_baseline(run_dir, series, 3, st), "lags")
  f_base <- run_engine_baseline(run_dir, series, 2, st)
  f_scen <- run_engine_scenario("test", scen, f_base, run_dir, series, 2, st)
  f_irf  <- run_engine_irfs(run_dir, series, 2, st)
  f_meta <- write_engine_meta(run_dir, list(run_id = "test_run"), st, c(f_base, f_irf, f_scen))
  expect_true(all(file.exists(c(f_base, f_scen, f_irf, f_meta))))
  expect_equal(basename(f_scen), "scenario_test.rds")

  base <- readRDS(f_base)
  expect_equal(base$paths, ref$baseline$paths, tolerance = 1e-10)
  expect_equal(base$inputs$last_observed, "2026Q2")
  sc <- readRDS(f_scen)
  # The CSV's columns come in another order; the conditions are the model's order.
  expect_equal(unname(sc$conditions[2, ]), c(NA, 0, 1))
  meta <- jsonlite::read_json(f_meta)
  expect_equal(meta$settings$seed, st$seed)
  expect_equal(meta$scenarios[[1]]$fixed_cells, 3L)
  expect_equal(meta$draws_md5, unname(tools::md5sum(file.path(run_dir, "draws.rds"))))
  expect_null(meta$scenarios[[1]]$tilting)

  # The same file with column a soft (its cell is then a target for the mean
  # level, not a fixed value): the scenario is tilted and meta.json says so.
  soft <- parse_soft(list(a = list(measure = "level", statistic = "mean")))
  f_soft <- run_engine_scenario("soft", scen, f_base, run_dir, series, 2, st, soft)
  sc_soft <- readRDS(f_soft)
  expect_true(all(is.na(sc_soft$conditions[, "a"])))
  expect_equal(sc_soft$tilt$targets$target, 2)
  expect_lt(max(abs(sc_soft$tilt$check$achieved - sc_soft$tilt$check$goal)), 1e-10)
  meta <- jsonlite::read_json(write_engine_meta(run_dir, list(run_id = "test_run"), st,
                                                c(f_base, f_irf, f_scen, f_soft)))
  expect_equal(meta$scenarios[[2]]$tilting$targets, 1L)

  # A base must fix no cell: "soft" above also fixes b and c, so it is refused.
  expect_error(run_engine_scenario("built", scen, f_base, run_dir, series, 2, st, base_file = f_soft),
               "fixes cells")
  # A base with soft conditions only, and a scenario built on it.
  soft_only <- tempfile(fileext = ".csv")
  on.exit(unlink(soft_only), add = TRUE)
  writeLines(c("quarter,a", "2026Q3,2"), soft_only)
  f_tilt <- run_engine_scenario("tilt", soft_only, f_base, run_dir, series, 2, st, soft)
  f_built <- run_engine_scenario("built", scen, f_base, run_dir, series, 2, st, base_file = f_tilt)
  built <- readRDS(f_built)
  expect_equal(built$base$id, "tilt")
  expect_identical(built$base$weights, readRDS(f_tilt)$tilt$weights)
  meta <- jsonlite::read_json(write_engine_meta(run_dir, list(run_id = "test_run"), st, c(f_base, f_irf, f_built)))
  expect_equal(meta$scenarios[[1]]$base$id, "tilt")

  # A baseline from other draws is refused.
  d <- readRDS(file.path(run_dir, "draws.rds")); d$sigma <- d$sigma * 1.01
  saveRDS(d, file.path(run_dir, "draws.rds"))
  expect_error(run_engine_scenario("test", scen, f_base, run_dir, series, 2, st), "other draws")
})
