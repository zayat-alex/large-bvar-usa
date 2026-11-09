# check_tilting.R: checks of entropic tilting on the model's own draws
# (R/tilting_checks.R), for both COVID runs.
#
# Run from the project root (so .Rprofile activates renv):
#   Rscript scripts/check_tilting.R run          (1,000 draws per run)
#   Rscript scripts/check_tilting.R run quick    (100 draws; tests the script)
#   Rscript scripts/check_tilting.R report
#
# Cases, each from the latest data (2026Q2 for vintage 2026m09):
#   1. The inflation anchor of Crump et al. (2021), Section 4.2, Figure 16b
#      (scenarios/anchor_pce_2pct): the median of 4-quarter headline PCE
#      inflation 12 quarters ahead tilted to 2.0%, with forecast$tilting
#      paths per draw.
#   2. The same anchor with 1, 10 and 20 paths per posterior draw: how the
#      effective sample size and the result depend on it.
#   3. Hard and soft together: scenarios/gs_inflation_oil_2026 with its core
#      PCE and core CPI columns as soft targets (median, yoy) and its oil
#      prices hard. The shipped scenario keeps every column hard.
# The paper's Figure 16a (inflation fixed exactly from 12 to 20 quarters
# ahead) is not reproduced: fixing a growth rate is a condition on a linear
# combination of cells, which scenario files do not support.
#
# Seeds: the engine's (mcmc$seed + forecast$seed_offset), as in the pipeline.
# Writes outputs/checks_tilting/ (tilting_checks.rds, PNG figures,
# tilting_checks.html).

for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)

args  <- commandArgs(trailingOnly = TRUE)
stage <- args[1]
quick <- isTRUE(args[2] == "quick")
if (is.na(stage)) stop("give a stage: run or report")

cfg      <- yaml::read_yaml("config/run.yml")
series   <- read_series_config("config/series.csv")
settings <- engine_settings(cfg, series$short_name)
out_dir  <- "outputs/checks_tilting"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
methods  <- unlist(cfg$covid$methods)
anchor_vars <- c("pce", "corepce", "rgdp", "urate", "gs2", "gs10")
mixed_vars  <- c("corepce", "corecpi", "wti", "rgdp", "urate", "gs2")
n_hist <- 12                                     # observed quarters shown before the forecast

#' Observed values of a fan's measure for the last n_hist quarters.
hist_measure <- function(history, v, measure) {
  if (measure == "yoy") yoy_history(history, v, n_hist) else utils::tail(history[, v], n_hist)
}

#' Figures from saved results (also redrawn by the report stage).
write_figures <- function(res) {
  for (m in names(res$methods)) {
    r <- res$methods[[m]]
    for (case in c("anchor", "mixed")) {
      cs <- r[[case]]
      vars <- names(cs$fans)
      png(file.path(out_dir, paste0(case, "_", m, ".png")), 1500, 950, res = 130)
      par(mfrow = c(2, 3), mar = c(4, 3, 2.5, 1), oma = c(2.5, 0, 0, 0))
      for (v in vars) {
        f <- cs$fans[[v]]
        tg <- cs$check[cs$check$variable == v, ]
        hard_h <- hard <- NULL
        # Hard cells are levels, so they are marked only on level panels.
        if (!is.null(cs$hard) && f$measure == "level" && any(!is.na(cs$hard[, v]))) {
          hard_h <- which(!is.na(cs$hard[, v])); hard <- cs$hard[hard_h, v]
        }
        unit <- if (f$measure == "yoy") "4-quarter % change" else "level"
        hv <- r$hist[[case]][[v]]
        if (v == "wti") {                         # oil: show dollars, not 100 ln
          f$untilted <- exp(f$untilted / 100); f$tilted <- exp(f$tilted / 100)
          hv <- exp(hv / 100)
          if (length(hard)) hard <- exp(hard / 100)
          unit <- "$ per barrel"
        }
        plot_tilt_fan(f, hv, paste0(v, " (", unit, ")"),
                      target_h = tg$h, target = tg$target, hard_h = hard_h, hard = hard)
      }
      legend_tilt_fan(hard = case == "mixed")
      dev.off()
    }
  }
}

# ---- Stage: run ----
if (stage == "run") {
  t_start <- proc.time()[["elapsed"]]
  idx <- read_scenario_index("scenarios/scenarios.yml")
  res <- list(run = list(time = Sys.time(), quick = quick, seed = settings$seed,
                         paths_per_draw = settings$tilting$paths_per_draw, tilting = settings$tilting),
              methods = list())
  for (m in methods) {
    run_id <- make_run_id(cfg, m, mcmc_settings(cfg), series)$run_id
    inputs <- engine_inputs(file.path(cfg$output_root, run_id), series, cfg$model$lags, settings$horizon)
    if (quick) {
      inputs$beta <- inputs$beta[1:100, , , drop = FALSE]
      inputs$sigma <- inputs$sigma[1:100, , , drop = FALSE]
    }
    free <- matrix(NA_real_, settings$horizon, length(inputs$variables),
                   dimnames = list(inputs$quarters, inputs$variables))
    base_mean <- forecast_conditional(inputs$beta, inputs$sigma, inputs$history, inputs$lags, settings$horizon,
                                      free, settings$seed, settings$method, shocks = FALSE)
    cat(m, ":", run_id, "|", dim(inputs$beta)[1], "draws\n")

    # Case 1: the anchor.
    a_file <- "scenarios/anchor_pce_2pct.csv"
    a_soft <- scenario_soft(idx, "anchor_pce_2pct")
    a_cond <- read_scenario(a_file, inputs$variables, inputs$last_observed, settings$horizon, a_soft)
    a_tg   <- soft_targets(a_file, inputs$variables, inputs$last_observed, settings$horizon, a_soft, series)
    anchor <- tilting_case(a_cond, a_tg, inputs, settings, base_mean, anchor_vars, series)
    cat("  anchor: ESS", round(anchor$concentration$ess), "of", anchor$concentration$paths,
        "| untilted median", round(anchor$check$untilted_median, 3), "\n")

    # Case 2: paths per draw.
    sens <- do.call(rbind, lapply(c(1L, 10L, 20L), function(k) {
      st <- settings; st$tilting$paths_per_draw <- k
      cs <- tilting_case(a_cond, a_tg, inputs, st, base_mean, "pce", series)
      f <- cs$fans$pce
      h <- a_tg$h[1]
      data.frame(paths_per_draw = k, cs$concentration, iterations = cs$iterations,
                 tilted_median = f$tilted[h, 4], tilted_q16 = f$tilted[h, 2], tilted_q84 = f$tilted[h, 6],
                 untilted_median = f$untilted[h, 4], seconds = cs$seconds)
    }))
    print(sens, digits = 3, row.names = FALSE)

    # Case 3: hard and soft together. The shipped file holds core inflation
    # as levels (100 ln of the price index, built from published 4-quarter
    # rates); a soft yoy column must hold the rates themselves, so they are
    # recovered from the levels (base quarter from the history or from the
    # file) and written to a temporary copy with oil unchanged.
    g_file <- tempfile("gs_inflation_oil_2026_soft_", fileext = ".csv")
    g_df <- utils::read.csv("scenarios/gs_inflation_oil_2026.csv", colClasses = "character", check.names = FALSE)
    g_lev <- read_scenario("scenarios/gs_inflation_oil_2026.csv", inputs$variables, inputs$last_observed,
                           settings$horizon)
    for (v in c("corepce", "corecpi")) {
      full <- c(inputs$history[, v], g_lev[, v])       # observed, then the scenario's levels
      h <- quarter_index(g_df$quarter) - quarter_index(inputs$last_observed)
      rows <- which(nzchar(trimws(g_df[[v]])))
      pos <- nrow(inputs$history) + h[rows]
      g_df[[v]][rows] <- format(round(100 * (exp((full[pos] - full[pos - 4]) / 100) - 1), 6))
    }
    utils::write.csv(g_df, g_file, row.names = FALSE, quote = FALSE)
    cat("  mixed: soft yoy targets recovered from the file's levels:",
        paste(g_df$quarter, g_df$corepce, g_df$corecpi, collapse = "; "), "\n")
    g_soft <- parse_soft(list(corepce = "yoy", corecpi = "yoy"))
    g_cond <- read_scenario(g_file, inputs$variables, inputs$last_observed, settings$horizon, g_soft)
    g_tg   <- soft_targets(g_file, inputs$variables, inputs$last_observed, settings$horizon, g_soft, series)
    mixed  <- tilting_case(g_cond, g_tg, inputs, settings, base_mean, mixed_vars, series, level_vars = "wti")
    mixed$hard <- g_cond
    cat("  mixed: ESS", round(mixed$concentration$ess), "of", mixed$concentration$paths,
        "| hard cells met to", signif(mixed$fixed_max_error, 3), "\n")

    hist <- list(
      anchor = lapply(stats::setNames(anchor_vars, anchor_vars), function(v)
        hist_measure(inputs$history, v, anchor$fans[[v]]$measure)),
      mixed = lapply(stats::setNames(mixed_vars, mixed_vars), function(v)
        hist_measure(inputs$history, v, mixed$fans[[v]]$measure)))
    res$methods[[m]] <- list(run_id = run_id, n_draws = dim(inputs$beta)[1], quarters = inputs$quarters,
                             last_observed = inputs$last_observed, anchor = anchor, sens = sens,
                             mixed = mixed, hist = hist)
  }
  res$run$seconds <- proc.time()[["elapsed"]] - t_start
  saveRDS(res, file.path(out_dir, "tilting_checks.rds"))
  write_figures(res)
  cat("done in", round(res$run$seconds / 60, 1), "minutes\n")
}

# ---- Stage: report ----
if (stage == "report") {
  write_figures(readRDS(file.path(out_dir, "tilting_checks.rds")))
  out_abs <- normalizePath(out_dir, winslash = "/")
  quarto::quarto_render("reports/tilting_checks.qmd", execute_params = list(output_dir = out_abs),
                        quiet = TRUE)
  target <- file.path(out_dir, "tilting_checks.html")
  if (file.exists(target)) file.remove(target)
  if (!file.rename("reports/tilting_checks.html", target)) stop("could not move the report to ", target)
  cat("Wrote", target, "\n")
}
