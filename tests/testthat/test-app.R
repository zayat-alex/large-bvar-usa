# test-app.R: the website (app/): it reads only the bundle, the numbers of
# each chart come from the bundle unchanged, differences are medians minus
# the baseline's median, the chart function draws what the settings ask,
# and every chart of the published bundle renders.
#
# The app's functions are loaded into their own environment, as Shiny does,
# not from the project's R/ folder. A small bundle is built from the stored
# engine reference outputs (fixtures/engine_reference.rds), as in
# test-export.R.

app_dir <- file.path(project_root(), "app")
app <- new.env()
for (f in list.files(file.path(app_dir, "R"), pattern = "[.]R$", full.names = TRUE)) sys.source(f, app)

ref <- readRDS(test_path("fixtures", "engine_reference.rds"))
ap_series <- data.frame(short_name = c("a", "b", "c"), variable = c("Series A", "Rate B", "Rate C"),
                        transform = c("log", "raw", "raw"))
ap_settings <- list(run = "test", horizon = 4L, objects = c("hard", "soft"),
                    spreads = data.frame(id = "c_b", a = "c", b = "b", label = "C minus B"), shown_in_levels = character(0))
ap_baseline <- ref$baseline
ap_baseline$inputs <- ref$inputs[c("last_observed", "quarters")]
ap_bundle <- build_bundle(ap_baseline, list(hard = ref$scenario, soft = ref$tilted), ref$irf, ref$inputs$history,
                          ap_series, data.frame(id = c("hard", "soft"), label = c("H", "S"), description = c("", "")),
                          ap_settings, ref$settings$quantiles)
ap_bundle$meta <- list(run_id = "test", data = list(last_observed = "2026Q2"), forecast = list())
ap_scen <- data.frame(id = c("baseline", "hard", "soft"), label = c("Base", "Hard", "Soft"),
                      colour = c("#0055A4", "#8E2C48", "#B5562B"))

test_that("the app reads only its own files and the bundle", {
  code <- unlist(lapply(list.files(app_dir, pattern = "[.]R$", recursive = TRUE, full.names = TRUE), readLines))
  expect_false(any(grepl("\\.\\./", code)))                        # nothing outside app/
  expect_false(any(grepl("source\\(", code)))                      # Shiny loads app/R/ itself
  expect_false(any(grepl("library\\(targets\\)|tar_read", code)))  # no pipeline at run time
})

test_that("the bundle is checked when loaded", {
  expect_true(app$check_bundle(ap_bundle))
  b <- ap_bundle; b$forecasts$q50 <- NULL
  expect_error(app$check_bundle(b), "lacks: q50")
  b <- ap_bundle; b$irf <- NULL
  expect_error(app$check_bundle(b), "no table 'irf'")
  expect_error(app$load_bundle(tempfile(fileext = ".rds")), "bundle not found")
})

test_that("each variable is shown in its display transform", {
  v <- ap_bundle$variables
  expect_equal(app$display_transform(v, "a", "yoy"), "yoy")
  expect_equal(app$display_transform(v, "a", "natural"), "natural")
  expect_equal(app$display_transform(v, "b", "yoy"), "level")
  expect_equal(app$display_transform(v, "c_b", "qoq_ann"), "level")
  v$display[v$variable == "a"] <- "natural"
  expect_equal(app$display_transform(v, "a", "yoy"), "natural")
})

test_that("chart series are the bundle's numbers, joined at the last observation", {
  s <- app$chart_series(ap_bundle, "a", "level", c("baseline", "hard"), history_quarters = 3)
  h <- ap_bundle$history[ap_bundle$history$variable == "a" & ap_bundle$history$transform == "level", ]
  expect_equal(s$history$value, utils::tail(h$value, 3))
  expect_equal(s$last$quarter, "2026Q2")
  for (obj in c("baseline", "hard")) {
    f <- s$forecasts[[obj]]
    expect_equal(unname(unlist(f[1, c("q05", "q50", "q95")])), rep(s$last$value, 3))
    want <- ap_bundle$forecasts[ap_bundle$forecasts$object == obj & ap_bundle$forecasts$variable == "a" &
                                ap_bundle$forecasts$transform == "level", ]
    expect_equal(f$q50[-1], want$q50)
    expect_equal(f$quarter[-1], want$quarter)
  }
  expect_equal(nrow(app$chart_series(ap_bundle, "a", "level", "baseline")$history), nrow(h))
})

test_that("differences are each median minus the baseline's median", {
  s <- app$chart_series(ap_bundle, "a", "yoy", c("baseline", "hard", "soft"))
  d <- app$median_differences(s, "hard")
  expect_named(d, c("baseline", "soft"))
  expect_equal(d$soft$value, s$forecasts$soft$q50 - s$forecasts$hard$q50)
  expect_equal(d$soft$value[1], 0)                                 # zero at the last observation
  expect_error(app$median_differences(s, "missing"), "not among")
})

test_that("the chart draws the lines and bands the settings ask for", {
  n_traces <- function(p) length(plotly::plotly_build(p)$x$data)
  s <- app$chart_series(ap_bundle, "a", "level", ap_scen$id)
  # Medians only: history + three medians.
  expect_equal(n_traces(app$fan_chart(s, ap_scen, character(0), "baseline")), 4)
  # Baseline bands as two shaded areas, another scenario's as four lines.
  expect_equal(n_traces(app$fan_chart(s, ap_scen, c("baseline", "hard"), "baseline")), 4 + 2 + 4)
  # Difference view: one line per scenario other than the baseline.
  expect_equal(n_traces(app$fan_chart(s, ap_scen, ap_scen$id, "soft", difference = TRUE)), 2)
  expect_equal(app$unit_label("rgdp", "yoy"), "% y/y")
  expect_equal(app$unit_label("rgdp", "qoq_ann", difference = TRUE), "pp (% q/q ann.)")
  expect_equal(app$unit_label("gs2", "level", difference = TRUE), "pp")
  expect_equal(app$unit_label("wti", "natural"), "$ per barrel")
})

test_that("impulse responses are the bundle's quantiles, scaled by the shock size", {
  irf <- ap_bundle$irf
  base <- irf[irf$identification == "recursive_a" & irf$variable == "a" & irf$transform == "level", ]
  x <- app$irf_series(ap_bundle, "recursive_a", "a", "level", horizon = 4)
  expect_equal(x$horizon, 0:4)
  expect_equal(x$q50, base$q50[base$horizon <= 4])
  # Linear: a shock of 0.25 scales every quantile; a negative shock swaps the bands.
  y <- app$irf_series(ap_bundle, "recursive_a", "a", "level", horizon = 4, size = 0.25)
  expect_equal(y$q84, 0.25 * x$q84)
  z <- app$irf_series(ap_bundle, "recursive_a", "a", "level", horizon = 4, size = -1)
  expect_equal(z$q05, -x$q95)
  expect_equal(z$q50, -x$q50)
  # A rate has no 4-quarter change: its level response is used.
  r <- app$irf_series(ap_bundle, "girf_a", "b", "yoy", horizon = 4)
  expect_equal(attr(r, "transform"), "level")
  expect_equal(app$irf_unit(ap_bundle$variables, "a", "level"), "%")
  expect_equal(app$irf_unit(ap_bundle$variables, "a", "yoy"), "pp, 4-quarter change")
  expect_equal(app$irf_unit(ap_bundle$variables, "b", "level"), "pp")
  d <- app$irf_data(ap_bundle, "girf_a", c("a", "b"), "yoy", horizon = 3, size = 2)
  expect_equal(nrow(d), 2 * 4)
  expect_equal(d$Median[d$Unit == "pp, 4-quarter change"],
               2 * irf$q50[irf$identification == "girf_a" & irf$variable == "a" & irf$transform == "yoy" &
                             irf$horizon <= 3])
})

test_that("shocks and their default responses come from the bundle", {
  s <- app$irf_shocks(ap_bundle)
  # No default shock in this bundle: recursive ones first, then generalized.
  expect_equal(s$id, c("recursive_a", "girf_a"))
  expect_equal(s$label, c("a (recursive)", "a (generalized)"))      # unlabelled variable: its name
  # A shock without a chosen set: its own variable, then the general set (those in the bundle).
  expect_equal(app$irf_default_responses(ap_bundle, "girf_a"), "a")
  expect_equal(app$irf_default_responses(ap_bundle, "girf_c_b"), "c_b")
})

test_that("the side panel groups the forecast tabs", {
  html <- as.character(app$side_panel())
  expect_match(html, "nav-group\">Forecasts")
  expect_match(html, "value=\"overview\" checked")
  # Overview, the topics, then Impulse responses, Settings and About.
  expect_equal(lengths(regmatches(html, gregexpr("name=\"page\"", html))), 1 + length(app$APP_TOPICS) + 3)
})

test_that("fixed cells are marked on level charts, at their fixed values", {
  m <- app$chart_markers(ap_bundle, "a", "level", c("baseline", "hard"))
  cn <- ap_bundle$conditions[ap_bundle$conditions$object == "hard" & ap_bundle$conditions$variable == "a" &
                               ap_bundle$conditions$transform == "level", ]
  expect_equal(m$value, cn$value)
  expect_true(all(m$kind == "fixed"))
  expect_equal(nrow(app$chart_markers(ap_bundle, "a", "yoy", c("baseline", "hard"))), 0)   # not on growth charts
  expect_equal(nrow(app$chart_markers(ap_bundle, "a", "level", "baseline")), 0)             # only shown scenarios
  # Median targets are marked where the chart shows their measure; mean targets are not.
  tg <- ap_bundle$targets
  b <- ap_bundle
  b$targets <- data.frame(object = "soft", variable = "a", quarter = "2027Q2", h = 4, measure = "yoy",
                          statistic = "median", target = 2.5)
  m <- app$chart_markers(b, "a", "yoy", "soft")
  expect_equal(m$value, 2.5); expect_equal(m$kind, "target")
  b$targets$statistic <- "mean"
  expect_equal(nrow(app$chart_markers(b, "a", "yoy", "soft")), 0)
})

test_that("the data table holds the charts' numbers and leaves out restricted history", {
  s <- list(show = c("baseline", "hard"), baseline = "baseline", mode = "forecast", growth = "yoy", history = "5y")
  d <- app$tab_data(ap_bundle, c("a", "b"), s, ap_scen)
  expect_setequal(unique(d$Series), c("Observed", "Base (median)", "Hard (median)"))
  f <- ap_bundle$forecasts
  want <- f[f$object == "hard" & f$variable == "a" & f$transform == "yoy", ]
  got <- d[which(d$Series == "Hard (median)" & d$Unit == "% y/y"), ]          # variable a (the yoy one)
  expect_equal(got$Value, want$q50)
  expect_equal(got$`90% low`, want$q05)
  # Restricted history is left out; the forecasts stay.
  old <- app$APP_NO_HISTORY_DOWNLOAD
  assign("APP_NO_HISTORY_DOWNLOAD", "b", envir = app)
  d2 <- app$tab_data(ap_bundle, "b", s, ap_scen)
  assign("APP_NO_HISTORY_DOWNLOAD", old, envir = app)
  expect_false("Observed" %in% d2$Series)
  expect_true("Base (median)" %in% d2$Series)
  # Difference view: medians minus the baseline's, no history, no bands.
  s$mode <- "difference"
  d3 <- app$tab_data(ap_bundle, "a", s, ap_scen)
  expect_equal(unique(d3$Series), "Hard minus Base")
  base <- f[f$object == "baseline" & f$variable == "a" & f$transform == "yoy", ]
  expect_equal(d3$Value, want$q50 - base$q50)
  expect_false("90% low" %in% names(d3))
})

test_that("every scenario of the bundle gets a fixed default colour", {
  cat <- app$scenario_catalog(ap_bundle)
  expect_equal(cat$id, c("baseline", "hard", "soft"))              # the default ("baseline") first
  expect_equal(cat$colour, app$APP_PALETTE[1:3])
  expect_equal(cat$default, c(TRUE, FALSE, FALSE))
  expect_equal(cat$label, c("Unconditional", "H", "S"))            # short labels from the bundle
  # More scenarios than palette slots: grey, not a recycled colour.
  many <- ap_bundle
  many$objects <- do.call(rbind, c(list(ap_bundle$objects), lapply(1:8, function(i) {
    o <- ap_bundle$objects[2, ]; o$id <- paste0("x", i); o
  })))
  cat <- app$scenario_catalog(many)
  expect_equal(cat$colour[9:11], rep(app$APP_PALETTE_EXTRA, 3))
  expect_false(anyDuplicated(cat$colour[1:8]) > 0)
})

# The published bundle (skipped if not published): labels and sources for
# every variable, every chart in every view, and the server's behaviour.
pub_path <- file.path(app_dir, "bundle", "bundle.rds")
pub <- if (file.exists(pub_path)) app$load_bundle(pub_path) else NULL

test_that("every variable of the published bundle has a label, a source and renders in every view", {
  skip_if(is.null(pub), "app/bundle/bundle.rds not published (scripts/publish_bundle.R)")
  expect_setequal(app$APP_VARIABLES$variable, pub$variables$variable)
  expect_true(all(app$APP_DEFAULT_SCENARIOS %in% pub$objects$id))
  expect_true(all(unlist(lapply(app$APP_TOPICS, `[[`, "variables")) %in% pub$variables$variable))
  expect_setequal(app$APP_SOURCES$variable, pub$variables$variable)
  expect_true(all(nzchar(app$APP_SOURCES$source)))
  page <- as.character(app$app_ui(pub))
  expect_match(page, "WTISPLC")
  expect_match(page, "Gürkaynak", fixed = TRUE)
  # Charts are built directly (the test harness would redraw every chart of
  # every tab on each input change): all variables in two views, a few in all.
  cat <- app$scenario_catalog(pub)
  s <- list(show = cat$id, bands = cat$id, baseline = "baseline", history = "5y")
  build <- function(v, mode, growth) plotly::plotly_build(app$make_chart(pub, v, c(s, mode = mode, growth = growth), cat))
  for (v in app$APP_VARIABLES$variable) for (mode in c("forecast", "difference")) expect_no_error(build(v, mode, "yoy"))
  for (v in c("rgdp", "wti", "gs2")) for (g in c("qoq_ann", "natural")) for (mode in c("forecast", "difference")) {
    expect_no_error(build(v, mode, g))
  }
  # Every shock in the bundle is listed (the default first), and its charts build.
  shocks <- app$irf_shocks(pub)
  expect_setequal(shocks$id, unique(pub$irf$identification))
  expect_equal(shocks$id[1], "recursive_gs2")
  expect_equal(shocks$label[1], "2y Treasury yield (recursive)")
  for (sh in shocks$id) for (v in app$irf_default_responses(pub, sh)) {
    expect_no_error(plotly::plotly_build(app$irf_chart(app$irf_series(pub, sh, v, "yoy", 20))))
  }
  # Text built from the bundle, not fixed: series count, lags, window, horizon.
  expect_match(page, "of 31 quarterly US series with 5 lags", fixed = TRUE)
  expect_match(page, "2020Q1–2021Q4 are left out", fixed = TRUE)
  expect_match(page, "the site shows 12", fixed = TRUE)
  # The oil scenarios' fixed cells appear on the oil chart, in dollars.
  m <- app$chart_markers(pub, "wti", "natural", "oil_rises")
  expect_equal(m$value[m$quarter == "2027Q1"], 110, tolerance = 1e-6)   # the CSV holds 4 decimals of 100 ln(price)
})

test_that("the server: settings, the Settings tab, tables and downloads", {
  skip_if(is.null(pub), "app/bundle/bundle.rds not published (scripts/publish_bundle.R)")
  ids <- app$APP_DEFAULT_SCENARIOS
  extra <- setdiff(pub$objects$id, ids)[1]
  shiny::testServer(app$app_server(pub), {
    session$setInputs(page = "overview", show = ids, bands = ids, baseline = "baseline", mode = "forecast",
                      growth = "yoy", history = "5y", overview_vars = c("rgdp", "sp500", "wti"),
                      irf_shock = "recursive_gs2", irf_units = "level", irf_horizon = "20", irf_size = 0.25,
                      irf_vars = c("rgdp", "gs2"))
    expect_equal(listed()$id, ids)
    # Impulse responses: a 25 bp recursive shock to the 2-year yield is 0.25 on impact.
    expect_match(as.character(output$grid_irf$html), "c_irf_gs2")
    expect_no_error(output$c_irf_rgdp)
    irf_csv <- utils::read.csv(output$dl_irf, check.names = FALSE)
    expect_equal(irf_csv$Median[irf_csv$Variable == "2y Treasury yield" & irf_csv$`Quarters after shock` == 0], 0.25)
    expect_equal(max(irf_csv$`Quarters after shock`), 20)
    expect_match(as.character(output$irf_size_unit$html), "pp")
    expect_no_error(output$c_overview_rgdp)
    expect_match(as.character(output$grid_overview$html), "c_overview_wti")
    expect_match(as.character(output$scenario_controls$html), "Bands")
    # Data table and CSV: the S&P 500's observed history is left out, its forecasts are not.
    expect_match(as.character(output$tbl_overview), "WTI crude oil")
    csv <- utils::read.csv(output$dl_overview, check.names = FALSE)
    expect_true(any(csv$Variable == "Real GDP" & csv$Series == "Observed"))
    expect_false(any(csv$Variable == "S&P 500" & csv$Series == "Observed"))
    expect_true(any(csv$Variable == "S&P 500" & grepl("median", csv$Series)))
    expect_true(all(csv$Unit[csv$Variable == "WTI crude oil"] == "$ per barrel"))
    # Settings tab, in one change: unlist a scenario, list another, recolour, rename;
    # and the difference view, which hides the Bands column.
    args <- list(list_oil_rises = FALSE, colour_baseline = "#000000", name_baseline = "Model forecast",
                 mode = "difference", data_overview_wti = 1)   # also opens the oil chart's Data pop-up
    args[[paste0("list_", extra)]] <- TRUE
    do.call(session$setInputs, args)
    session$elapse(600)
    expect_false("oil_rises" %in% listed()$id)
    expect_true(extra %in% listed()$id)
    expect_equal(listed()$colour[listed()$id == "baseline"], "#000000")
    expect_equal(listed()$label[listed()$id == "baseline"], "Model forecast")
    html <- as.character(output$scenario_controls$html)
    expect_match(html, "Model forecast")
    expect_false(grepl(">Bands<", html))
    # The Data pop-up: one variable, in the current (difference) view.
    expect_equal(modal_variable(), "wti")
    expect_match(as.character(output$modal_tbl), "WTI crude oil")
    pop <- utils::read.csv(output$modal_dl, check.names = FALSE)
    expect_true(all(pop$Variable == "WTI crude oil"))
    expect_true(all(grepl("minus Model forecast$", pop$Series)))
    csv <- utils::read.csv(output$dl_overview, check.names = FALSE)
    expect_true(all(grepl("minus Model forecast$", csv$Series)))
    # A blank name falls back; reset restores the defaults.
    session$setInputs(name_baseline = "  ", reset_scenarios = 1)
    session$elapse(600)
    expect_equal(listed()$id, ids)
    expect_equal(listed()$label[1], "Unconditional")
    expect_equal(listed()$colour[1], app$APP_PALETTE[1])
  })
})
