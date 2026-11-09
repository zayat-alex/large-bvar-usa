# check_irf_review.R: impulse-response review (spec 9.5; R/irf_review.R).
#
# Run from the project root (so .Rprofile activates renv):
#   Rscript scripts/check_irf_review.R figures   (response figures only)
#   Rscript scripts/check_irf_review.R run       (figures, spread widening, tables)
#   Rscript scripts/check_irf_review.R report
#
# Reads each production run's engine output (outputs/<run_id>/engine/irf.rds)
# and draws every identification's responses, both COVID methods overlaid,
# up to 20 quarters. Unit shocks, as computed by the engine; the tables
# scale them to the sizes of the expectations (oil +10%, S&P 500 -10%) and
# of the paper's figures (gs2 +25 bp). `run` also computes the Baa spread
# widening (Baa +1 pp with the 10-year yield at its baseline on impact, a
# conditional forecast) from each run's draws.
#
# Writes outputs/checks_irf_review/ (irf_review.rds, PNG figures,
# irf_review.html).

for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)

args  <- commandArgs(trailingOnly = TRUE)
stage <- args[1]
if (is.na(stage)) stop("give a stage: figures, run or report")

cfg     <- yaml::read_yaml("config/run.yml")
out_dir <- "outputs/checks_irf_review"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
max_h <- 20

# Responses shown per shock: the variables of the spec 9.5 expectations and
# the paper's Figures 6, 8 and 9.
shocks <- list(
  recursive_gs2 = list(title = "Recursive 2-year yield shock (+1 pp; slow variables fixed on impact)",
                       responses = c("rgdp", "rpce", "rnresinv", "houst", "indpro", "payems", "urate",
                                     "corepce", "gs5", "gs10", "aaa", "baa", "sp500", "usd", "wti", "rvol")),
  girf_gs2  = list(title = "2-year yield, generalized IRF (+1 pp; nothing else fixed)",
                   responses = c("rgdp", "rpce", "rnresinv", "houst", "indpro", "payems", "urate",
                                 "corepce", "gs5", "gs10", "aaa", "baa", "sp500", "usd", "wti", "rvol")),
  girf_rgdp = list(title = "Real GDP, generalized IRF (+1 log point)",
                   responses = c("rgdp", "rpce", "rnresinv", "payems", "urate", "indpro", "pce", "corepce",
                                 "gs2", "gs10", "baa", "sp500")),
  girf_wti  = list(title = "Oil price (WTI), generalized IRF (+1 log point)",
                   responses = c("wti", "pce", "corepce", "cpi", "corecpi", "rgdp", "rpce", "urate",
                                 "gs2", "gs10", "sp500", "gsci_nonenergy")),
  girf_sp500 = list(title = "S&P 500, generalized IRF (+1 log point)",
                    responses = c("sp500", "rgdp", "rpce", "rnresinv", "payems", "urate", "corepce",
                                  "gs2", "gs10", "baa", "rvol", "umich")),
  girf_baa  = list(title = "Baa yield, generalized IRF (+1 pp)",
                   responses = c("baa", "aaa", "gs2", "gs10", "rgdp", "rnresinv", "indpro", "payems",
                                 "urate", "corepce", "sp500", "rvol")))

if (stage %in% c("figures", "run")) {
  runs <- list.files(cfg$output_root, pattern = "^[0-9]{4}m[0-9]{2}_.+_full_[0-9a-f]{8}$", full.names = TRUE)
  methods <- vapply(runs, function(d) jsonlite::read_json(file.path(d, "meta.json"))$covid$method, character(1))
  irf <- setNames(lapply(runs, function(d) readRDS(file.path(d, "engine", "irf.rds"))), methods)
  method_order <- c("none", "precovid_params", "exclude_window")
  irf <- irf[method_order[method_order %in% names(irf)]]
  for (s in names(shocks)) {
    q <- lapply(irf, function(x) x$quantiles[[s]])
    n_row <- ceiling(length(shocks[[s]]$responses) / 4)
    png(file.path(out_dir, paste0("irf_", s, ".png")), 1700, 330 * n_row + 120, res = 130)
    plot_irf_methods(q, shocks[[s]]$responses, max_h, shocks[[s]]$title)
    dev.off()
    cat("wrote", file.path(out_dir, paste0("irf_", s, ".png")), "\n")
  }
}

# ---- Stage: run ----
if (stage == "run") {
  series <- read_series_config("config/series.csv")
  vars <- series$short_name
  is_log <- series$transform == "log"
  lags <- cfg$model$lags

  # Spread widening from each run's draws, and its cross-check on three draws.
  spread <- list(); spread_check <- list()
  for (k in seq_along(runs)) {
    m <- methods[[k]]
    dr <- readRDS(file.path(runs[k], "draws.rds"))
    hist <- est_data_to_matrix(readRDS(file.path(runs[k], "covid.rds"))$origin_data, series)
    sd_draws <- spread_irf_draws(dr$beta, dr$sigma, lags, cfg$irf$horizon, vars, "baa", "gs10", 1)
    spread[[m]] <- list(draws = sd_draws, quantiles = irf_quantiles(sd_draws, unlist(cfg$quantiles)))
    spread_check[[m]] <- max(sapply(c(1, 500, 1000), function(i)
      spread_vs_conditional(dr$beta[i, , ], dr$sigma[i, , ], hist, lags, 21, "baa", "gs10", 1)))
    cat(m, ": spread widening vs exact conditional means, largest difference", signif(spread_check[[m]], 3), "\n")
  }
  spread <- spread[names(irf)]
  q_all <- lapply(names(irf), function(m) c(irf[[m]]$quantiles, list(spread_baa = spread[[m]]$quantiles)))
  names(q_all) <- names(irf)
  d_all <- lapply(names(irf), function(m) c(irf[[m]]$draws, list(spread_baa = spread[[m]]$draws)))
  names(d_all) <- names(irf)
  png(file.path(out_dir, "irf_spread_baa.png"), 1700, 330 * 3 + 120, res = 130)
  plot_irf_methods(lapply(spread, `[[`, "quantiles"),
                   c("baa", "gs10", "aaa", "gs2", "rgdp", "rnresinv", "indpro", "payems", "urate", "corepce",
                     "sp500", "rvol"), max_h,
                   "Baa spread widening (Baa +1 pp, 10-year yield at baseline on impact; conditional forecast)")
  dev.off()

  # Expectations of spec 9.5: sign within 0-12 quarters, scaled shock sizes.
  ex <- function(shock, scale, v, sg) data.frame(shock = shock, variable = v, sign = sg, h_from = 0, h_to = 12,
                                                 scale = scale)
  rows <- rbind(
    ex("recursive_gs2", 1, c("houst", "indpro", "payems", "urate", "gs5", "gs10", "sp500", "usd", "rgdp"),
       c(-1, -1, -1, 1, 1, 1, -1, 1, -1)),
    ex("girf_wti", 10, c("pce", "corepce", "cpi", "corecpi", "rgdp", "urate"), c(1, 1, 1, 1, -1, 1)),
    ex("girf_sp500", -10, c("rgdp", "rpce", "payems", "urate", "gs2", "gs10", "rvol"), c(-1, -1, -1, 1, -1, -1, 1)),
    ex("girf_rgdp", 1, c("payems", "urate", "indpro", "pce", "corepce"), c(1, -1, 1, 1, 1)),
    ex("spread_baa", 1, c("rgdp", "indpro", "payems", "urate", "rnresinv"), c(-1, -1, -1, 1, -1)))
  expect <- do.call(rbind, lapply(names(q_all), function(m) do.call(rbind, lapply(unique(rows$shock), function(s)
    cbind(method = m, expectation_results(q_all[[m]][[s]], rows[rows$shock == s, ]))))))

  # Paper comparison (Crump et al. 2021, Figures 6, 8, 9): medians read from
  # the charts by eye, approximate. Paper shock sizes: GDP +1%, gs2 +25 bp.
  pr <- function(fig, shock, v, stat, h_from, h_to, scale, paper) data.frame(figure = fig, shock = shock,
    variable = v, stat = stat, h_from = h_from, h_to = h_to, scale = scale, paper = paper)
  paper <- rbind(
    pr("6", "girf_rgdp", "rgdp", "max", 0, 12, 1, 1.1), pr("6", "girf_rgdp", "rgdp", "at", 12, 12, 1, 0.5),
    pr("6", "girf_rgdp", "payems", "max", 0, 12, 1, 0.65), pr("6", "girf_rgdp", "urate", "min", 0, 12, 1, -0.35),
    pr("6", "girf_rgdp", "rnresinv", "max", 0, 12, 1, 2.9), pr("6", "girf_rgdp", "pce", "at", 12, 12, 1, 0.37),
    pr("6", "girf_rgdp", "corepce", "at", 12, 12, 1, 0.26), pr("6", "girf_rgdp", "gs2", "max", 0, 12, 1, 0.55),
    pr("6", "girf_rgdp", "gs10", "max", 0, 12, 1, 0.33), pr("6", "girf_rgdp", "baa", "max", 0, 12, 1, 0.33),
    pr("6", "girf_rgdp", "sp500", "max", 0, 4, 1, 2.8),
    pr("8", "girf_gs2", "rgdp", "max", 0, 4, 0.25, 0.09), pr("8", "girf_gs2", "rgdp", "at", 12, 12, 0.25, -0.08),
    pr("8", "girf_gs2", "urate", "min", 0, 4, 0.25, -0.04), pr("8", "girf_gs2", "urate", "at", 12, 12, 0.25, 0.06),
    pr("8", "girf_gs2", "baa", "max", 0, 12, 0.25, 0.14), pr("8", "girf_gs2", "sp500", "min", 0, 12, 0.25, -1.0),
    pr("9", "recursive_gs2", "rgdp", "at", 1, 1, 0.25, 0), pr("9", "recursive_gs2", "rgdp", "min", 0, 12, 0.25, -0.09),
    pr("9", "recursive_gs2", "urate", "at", 12, 12, 0.25, 0.07),
    pr("9", "recursive_gs2", "payems", "min", 0, 12, 0.25, -0.12),
    pr("9", "recursive_gs2", "baa", "max", 0, 12, 0.25, 0.14), pr("9", "recursive_gs2", "sp500", "min", 0, 12, 0.25, -0.8),
    pr("9", "recursive_gs2", "usd", "at", 12, 12, 0.25, 0.6))
  for (m in names(q_all)) {
    st <- lapply(seq_len(nrow(paper)), function(k)
      with(paper[k, ], paper_stat(q_all[[m]][[shock]], variable, stat, h_from, h_to, scale)))
    paper[[paste0(m, "_value")]] <- sapply(st, `[[`, "value")
    paper[[paste0(m, "_horizon")]] <- sapply(st, `[[`, "horizon")
  }

  # Bands one quarter after impact, recursive gs2 shock (the h1 blips).
  h1 <- do.call(rbind, lapply(names(q_all), function(m) {
    q <- q_all[[m]]$recursive_gs2
    d <- q[q$horizon == 1 & q$variable %in% c("rgdp", "urate", "payems", "indpro"), ]
    data.frame(method = m, variable = d$variable, q05 = d$q05, q16 = d$q16, q50 = d$q50, q84 = d$q84, q95 = d$q95)
  }))

  # Growth-rate reading: effect on YoY growth (log series) or level (rates).
  gspec <- list(recursive_gs2 = list(scale = 1, v = c("rgdp", "rpce", "payems", "urate", "corepce", "houst")),
                girf_wti = list(scale = 10, v = c("pce", "corepce", "cpi", "rgdp", "urate")),
                girf_sp500 = list(scale = -10, v = c("rgdp", "rpce", "payems", "urate", "corepce")),
                girf_rgdp = list(scale = 1, v = c("rgdp", "payems", "urate", "pce", "corepce")),
                spread_baa = list(scale = 1, v = c("rgdp", "rnresinv", "payems", "urate", "corepce")))
  growth <- do.call(rbind, lapply(names(d_all), function(m) do.call(rbind, lapply(names(gspec), function(s)
    cbind(method = m, shock = s, irf_growth_table(d_all[[m]][[s]], gspec[[s]]$v, is_log, c(4, 8, 12),
                                                  gspec[[s]]$scale))))))

  res <- list(time = Sys.time(), runs = setNames(basename(runs), methods), spread_check = spread_check,
              expect = expect, paper = paper, h1 = h1, growth = growth, horizon_shown = max_h)
  saveRDS(res, file.path(out_dir, "irf_review.rds"))
  print(expect[, c("method", "shock", "variable", "horizon", "median", "q16", "q84", "verdict")], digits = 2)
  print(paper, digits = 2)
  print(h1, digits = 2)
}

# ---- Stage: report ----
if (stage == "report") {
  out_abs <- normalizePath(out_dir, winslash = "/")
  quarto::quarto_render("reports/irf_review.qmd", execute_params = list(output_dir = out_abs), quiet = TRUE)
  target <- file.path(out_dir, "irf_review.html")
  if (file.exists(target)) file.remove(target)
  if (!file.rename("reports/irf_review.html", target)) stop("could not move the report to ", target)
  cat("Wrote", target, "\n")
}
