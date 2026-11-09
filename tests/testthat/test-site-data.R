# test-site-data.R: the static website's data (site/R/site_data.R). The
# JSON files must hold exactly what the Shiny app plots, at full precision
# (17 significant digits, which read back as the same doubles): the same transforms, history, forecast quantiles,
# markers, impulse responses and settings. The app's own functions
# (chart_series(), chart_markers(), irf_series()) are the reference.
#
# Uses the published bundle (app/bundle/bundle.rds); skipped if it has not
# been published (scripts/publish_bundle.R).

site_app <- new.env()
for (f in list.files(file.path(project_root(), "app", "R"), pattern = "[.]R$", full.names = TRUE)) sys.source(f, site_app)
site_fns <- new.env()
sys.source(file.path(project_root(), "site", "R", "site_data.R"), site_fns)

site_pub_path <- file.path(project_root(), "app", "bundle", "bundle.rds")
site_pub <- if (file.exists(site_pub_path)) site_app$load_bundle(site_pub_path) else NULL
# Largest difference relative to the value (at least 1). 17 significant
# digits read back exactly; the tolerance allows for the app's own
# arithmetic on the bundle (for example a growth rate computed again).
site_rel <- function(a, b) max(abs(a - b) / pmax(1, abs(b)), na.rm = TRUE)
site_tol <- 1e-13

# The written files, read back as the browser would get them.
site_json <- if (!is.null(site_pub)) {
  dir <- tempfile("site-data-")
  paths <- site_fns$write_site_data(site_pub, site_app, dir)
  list(core = jsonlite::read_json(paths[1], simplifyVector = TRUE),
       irf = jsonlite::read_json(paths[2], simplifyVector = TRUE))
}

test_that("quarter sequences are consecutive and cross years", {
  expect_equal(site_fns$site_quarter_seq("2026Q3", "2027Q2"), c("2026Q3", "2026Q4", "2027Q1", "2027Q2"))
  expect_equal(site_fns$site_quarter_seq("1973Q1", "1973Q1"), "1973Q1")
})

test_that("the site keeps exactly the transforms the app can display", {
  skip_if(is.null(site_pub), "app/bundle/bundle.rds not published (scripts/publish_bundle.R)")
  v <- site_pub$variables
  for (k in seq_len(nrow(v))) {
    shown <- unique(vapply(c("yoy", "qoq_ann", "natural"),
                           function(g) site_app$display_transform(v, v$variable[k], g), ""))
    expect_setequal(site_fns$site_transforms(v$display[k]), shown)
  }
  expect_error(site_fns$site_transforms("other"), "unknown display mark")
})

test_that("history and forecasts equal what the app plots", {
  skip_if(is.null(site_pub), "app/bundle/bundle.rds not published (scripts/publish_bundle.R)")
  core <- site_json$core
  ids <- site_pub$objects$id
  quarters <- core$meta$forecast_quarters
  n_checked <- 0
  for (k in seq_len(nrow(core$variables))) {
    v <- core$variables$variable[k]
    for (tr in site_fns$site_transforms(core$variables$display[k])) {
      ref <- site_app$chart_series(site_pub, v, tr, ids)
      h <- core$history[[v]][[tr]]
      q <- site_fns$site_quarter_seq(h$start, utils::tail(ref$history$quarter, 1))
      expect_equal(q, ref$history$quarter)
      expect_lte(site_rel(h$values, ref$history$value), site_tol)
      expect_identical(is.na(h$values), is.na(ref$history$value))
      for (obj in ids) {
        f <- ref$forecasts[[obj]][-1, ]                               # without the anchor (last observation)
        expect_equal(f$quarter, quarters)
        for (qq in site_fns$SITE_QUANTILES) {
          expect_lte(site_rel(core$forecasts[[obj]][[v]][[tr]][[qq]], f[[qq]]), site_tol)
        }
        n_checked <- n_checked + 1
      }
    }
  }
  # 34 variables in 1 or 3 transforms, for every scenario.
  expect_gt(n_checked, 100)
})

test_that("fixed cells and targets equal the app's chart markers", {
  skip_if(is.null(site_pub), "app/bundle/bundle.rds not published (scripts/publish_bundle.R)")
  core <- site_json$core
  ids <- site_pub$objects$id
  site_m <- rbind(data.frame(core$conditions, kind = "fixed"), data.frame(core$targets, kind = "target"))
  for (k in seq_len(nrow(core$variables))) {
    v <- core$variables$variable[k]
    for (tr in site_fns$site_transforms(core$variables$display[k])) {
      ref <- site_app$chart_markers(site_pub, v, tr, ids)
      got <- site_m[site_m$variable == v & site_m$transform == tr, ]
      expect_equal(nrow(got), nrow(ref))
      if (nrow(ref)) {
        ref <- ref[order(ref$object, ref$quarter, ref$kind), ]
        got <- got[order(got$object, got$quarter, got$kind), ]
        expect_equal(got$object, ref$object)
        expect_equal(got$quarter, ref$quarter)
        expect_equal(got$kind, ref$kind)
        expect_lte(site_rel(got$value, ref$value), site_tol)
      }
    }
  }
  # The published scenarios fix cells and set targets the charts show.
  expect_gt(nrow(core$conditions), 0)
  expect_gt(nrow(core$targets), 0)
})

test_that("impulse responses equal what the app plots, to the longest horizon offered", {
  skip_if(is.null(site_pub), "app/bundle/bundle.rds not published (scripts/publish_bundle.R)")
  irf <- site_json$irf
  hmax <- max(site_app$APP_IRF_HORIZONS)
  expect_equal(irf$horizon_max, hmax)
  expect_equal(irf$shocks, site_app$irf_shocks(site_pub))
  for (id in irf$shocks$id) {
    for (v in names(irf$responses[[id]])) {
      for (tr in c("level", "yoy")) {
        ref <- site_app$irf_series(site_pub, id, v, tr, hmax)
        got <- irf$responses[[id]][[v]][[attr(ref, "transform")]]
        # A variable has a 4-quarter response in the site exactly when the app shows one.
        expect_identical(tr %in% names(irf$responses[[id]][[v]]), attr(ref, "transform") == tr)
        expect_equal(ref$horizon, 0:hmax)
        for (qq in site_fns$SITE_QUANTILES) expect_lte(site_rel(got[[qq]], ref[[qq]]), site_tol)
      }
    }
    expect_setequal(names(irf$responses[[id]]), unique(site_pub$irf$variable[site_pub$irf$identification == id]))
  }
  expect_equal(length(irf$shocks$id), length(unique(site_pub$irf$identification)))
})

test_that("settings, variables and scenarios are the app's", {
  skip_if(is.null(site_pub), "app/bundle/bundle.rds not published (scripts/publish_bundle.R)")
  core <- site_json$core
  s <- core$settings
  expect_equal(s$palette, site_app$APP_PALETTE)
  expect_equal(s$overview_default, site_app$APP_OVERVIEW_DEFAULT)
  expect_equal(s$no_history_download, site_app$APP_NO_HISTORY_DOWNLOAD)
  expect_equal(s$history$name, names(site_app$APP_HISTORY))
  expect_equal(s$history$quarters, c(20, 40, 80, NA))
  expect_equal(s$topics$id, vapply(site_app$APP_TOPICS, `[[`, "", "id"))
  expect_equal(s$irf_horizons, site_app$APP_IRF_HORIZONS)
  expect_equal(core$scenarios, site_app$scenario_catalog(site_pub))
  expect_equal(core$variables$label, site_app$var_label(core$variables$variable))
  expect_equal(core$meta$last_observed, site_pub$meta$data$last_observed)
  expect_equal(core$meta$run_id, site_pub$meta$run_id)
})

test_that("a history with a missing quarter is an error", {
  skip_if(is.null(site_pub), "app/bundle/bundle.rds not published (scripts/publish_bundle.R)")
  b <- site_pub
  drop <- which(b$history$variable == "rgdp" & b$history$transform == "yoy")[10]
  b$history <- b$history[-drop, ]
  expect_error(site_fns$site_history(b, site_fns$site_variables(b, site_app)), "skips quarters")
})
