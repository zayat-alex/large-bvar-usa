# test-export.R: the export stage (R/export.R): display transforms by hand,
# growth rates as the inverse of path_from_growth(), the bundle built from
# the stored engine reference outputs (fixtures/engine_reference.rds) against
# the engine's own quantiles, spreads, the history, impulse responses, the
# bundle's schema and the export settings.

ref <- readRDS(test_path("fixtures", "engine_reference.rds"))
# The reference VAR has three variables: a log series and two rates.
ex_series <- data.frame(short_name = c("a", "b", "c"), variable = c("Series A", "Rate B", "Rate C"),
                        transform = c("log", "raw", "raw"))
ex_settings <- list(run = "test", horizon = 4L, objects = c("hard", "soft"),
                    spreads = data.frame(id = "c_b", a = "c", b = "b", label = "C minus B"),
                    shown_in_levels = "a")
ex_index <- data.frame(id = c("hard", "soft"), label = c("Hard", "Soft"), description = c("d1", "d2"))
ex_baseline <- ref$baseline
ex_baseline$inputs <- ref$inputs[c("last_observed", "quarters")]
ex_bundle <- build_bundle(ex_baseline, list(hard = ref$scenario, soft = ref$tilted), ref$irf,
                           ref$inputs$history, ex_series, ex_index, ex_settings, ref$settings$quantiles)
qcols <- sprintf("q%02d", round(100 * ref$settings$quantiles))
qm <- function(tab) unname(as.matrix(tab[qcols]))   # quantile columns as a plain matrix
qv <- function(x) unname(unlist(x))
pick <- function(tab, ...) {
  f <- list(...)
  keep <- Reduce(`&`, lapply(names(f), function(n) tab[[n]] == f[[n]]))
  tab[keep, , drop = FALSE]
}

test_that("display transforms match a hand calculation", {
  base <- 100 * log(c(97, 98, 99, 100))          # quarters T-3..T
  y <- matrix(100 * log(c(102, 104)), 1, dimnames = list(NULL, c("2026Q3", "2026Q4")))
  expect_equal(transform_paths(y, base, "level"), y)
  expect_equal(as.vector(transform_paths(y, base, "natural")), c(102, 104))
  expect_equal(as.vector(transform_paths(y, base, "qoq_ann")), 100 * (c(102 / 100, 104 / 102)^4 - 1))
  expect_equal(as.vector(transform_paths(y, base, "yoy")), 100 * (c(102 / 97, 104 / 98) - 1))
  expect_identical(colnames(transform_paths(y, base, "yoy")), c("2026Q3", "2026Q4"))
  # Missing base quarters make the growth rates that need them NA.
  expect_equal(unname(is.na(transform_paths(y, c(NA, NA, NA, base[4]), "yoy")[1, ])), c(TRUE, TRUE))
  expect_false(anyNA(transform_paths(y, c(NA, NA, NA, base[4]), "qoq_ann")))
  expect_error(transform_paths(y, base, "log_diff"), "unknown transform")
})

test_that("exact growth rates invert path_from_growth()", {
  h <- ref$inputs$history
  g <- c(2, 3.5, -1, 4, 0.5)
  for (type in c("qoq_ann", "yoy")) {
    y <- path_from_growth("a", g, type, h, ex_series)
    back <- transform_paths(matrix(y, 1), utils::tail(h[, "a"], 4), type)
    expect_equal(as.vector(back), g, tolerance = 1e-12)
  }
})

test_that("forecast quantiles in levels equal the engine's, weighted for a tilted object", {
  quarters <- ref$inputs$quarters[1:4]
  for (obj in c("baseline", "hard", "soft")) {
    eng <- switch(obj, baseline = ref$baseline, hard = ref$scenario, soft = ref$tilted)$quantiles
    for (v in c("a", "b", "c")) {
      got <- pick(ex_bundle$forecasts, object = obj, variable = v, transform = "level")
      want <- eng[eng$variable == v & eng$quarter %in% quarters, ]
      expect_equal(got$quarter, quarters)
      expect_equal(qm(got), qm(want))
    }
  }
})

test_that("growth quantiles are quantiles of each path's growth rate", {
  paths <- ref$baseline$paths
  base <- utils::tail(ref$inputs$history[, "a"], 4)
  h <- 3
  yoy <- 100 * (exp((paths[, h, "a"] - base[h]) / 100) - 1)       # base quarter T + 3 - 4 = T - 1
  got <- pick(ex_bundle$forecasts, object = "baseline", variable = "a", transform = "yoy")
  expect_equal(qv(got[h, qcols]), stats::quantile(yoy, ref$settings$quantiles, names = FALSE))
  # Weighted for the tilted object.
  w <- ref$tilted$tilt$weights
  qoq <- 100 * (exp(4 * (ref$tilted$paths[, h, "a"] - ref$tilted$paths[, h - 1, "a"]) / 100) - 1)
  got <- pick(ex_bundle$forecasts, object = "soft", variable = "a", transform = "qoq_ann")
  expect_equal(qv(got[h, qcols]), weighted_quantiles(qoq, w, ref$settings$quantiles))
})

test_that("spreads are formed path by path", {
  d <- ref$baseline$paths[, , "c"] - ref$baseline$paths[, , "b"]
  got <- pick(ex_bundle$forecasts, object = "baseline", variable = "c_b")
  expect_identical(unique(got$transform), "level")
  for (h in 1:4) {
    expect_equal(qv(got[h, qcols]), stats::quantile(d[, h], ref$settings$quantiles, names = FALSE))
  }
})

test_that("every object has the same rows, so the app can difference any two medians", {
  # The app shows differences as one object's median minus another's, row
  # by row; that needs the same variables, transforms and quarters in the
  # same order for every object.
  key <- function(obj) {
    f <- pick(ex_bundle$forecasts, object = obj)
    paste(f$variable, f$transform, f$quarter)
  }
  expect_identical(key("hard"), key("baseline"))
  expect_identical(key("soft"), key("baseline"))
  d <- pick(ex_bundle$forecasts, object = "hard")$q50 - pick(ex_bundle$forecasts, object = "baseline")$q50
  expect_false(anyNA(d))
})

test_that("the history is the observed data in each transform", {
  h <- ref$inputs$history
  lev <- pick(ex_bundle$history, variable = "a", transform = "level")
  expect_equal(lev$value, unname(h[, "a"]))
  expect_equal(lev$quarter, rownames(h))
  # Growth rates start once their base quarter is in the sample.
  yoy <- pick(ex_bundle$history, variable = "a", transform = "yoy")
  expect_equal(yoy$quarter, rownames(h)[5])
  expect_equal(yoy$value, 100 * (exp((h[5, "a"] - h[1, "a"]) / 100) - 1))
  expect_equal(pick(ex_bundle$history, variable = "c_b")$value, unname(h[, "c"] - h[, "b"]))
  # The forecast starts the quarter after the history ends.
  expect_identical(ex_bundle$forecasts$quarter[1], ref$inputs$quarters[1])
})

test_that("impulse responses: levels equal the engine's, yoy is the four-quarter change", {
  for (id in names(ref$irf$draws)) {
    got <- pick(ex_bundle$irf, identification = id, variable = "a", transform = "level")
    want <- ref$irf$quantiles[[id]]
    want <- want[want$variable == "a", ]
    expect_equal(qm(got), qm(want))
    r <- ref$irf$draws[[id]][, , "a"]
    d <- r[, 6] - r[, 2]                         # horizon 5 minus horizon 1
    yoy <- pick(ex_bundle$irf, identification = id, variable = "a", transform = "yoy")
    expect_equal(qv(yoy[yoy$horizon == 5, qcols]), stats::quantile(d, ref$settings$quantiles, names = FALSE))
    # Before four quarters have passed, the change is the response itself.
    expect_equal(qm(yoy[yoy$horizon <= 3, ]), qm(got[got$horizon <= 3, ]))
  }
  expect_identical(unique(pick(ex_bundle$irf, variable = "b")$transform), "level")
  expect_setequal(unique(ex_bundle$irf$scheme), c("girf", "recursive"))
})

test_that("generalized IRFs of the other variables are summarized from their draws", {
  vars <- export_variables(ex_series, ex_settings$spreads)
  probs <- ref$settings$quantiles
  x <- irf_extra_quantiles(ref$inputs, 8, skip = "a", vars, probs)
  expect_setequal(unique(x$identification), c("girf_b", "girf_c"))
  # The same quantiles as summarizing the full draws (4-quarter changes and
  # spreads taken per draw).
  d <- irf_draws(ref$inputs$beta, ref$inputs$sigma, ref$inputs$lags, 8, ref$inputs$variables, "c", NULL)
  want <- export_irf(d, vars, probs)
  got <- x[x$identification == "girf_c", ]
  rownames(got) <- NULL
  expect_equal(got, want)
  # Appended to the bundle's responses.
  b <- build_bundle(ex_baseline, list(hard = ref$scenario, soft = ref$tilted), ref$irf, ref$inputs$history,
                    ex_series, ex_index, ex_settings, probs, irf_extra = x)
  expect_setequal(unique(b$irf$identification), c(names(ref$irf$draws), "girf_b", "girf_c"))
})

test_that("the bundle holds only the tables the app reads", {
  expect_named(ex_bundle, c("variables", "objects", "history", "forecasts", "conditions", "targets", "irf"))
  id_cols <- c("object", "variable", "transform", "quarter", "h", "date")
  expect_named(ex_bundle$forecasts, c(id_cols, qcols))
  expect_named(ex_bundle$history, c("variable", "transform", "quarter", "date", "value"))
  expect_named(ex_bundle$irf, c("identification", "scheme", "shock", "variable", "transform", "horizon", qcols))
  expect_equal(ex_bundle$objects$id, c("baseline", "hard", "soft"))
  expect_equal(ex_bundle$objects$type, c("baseline", "hard", "tilted"))
  expect_equal(ex_bundle$objects$short_label, c("Unconditional", "Hard", "Soft"))   # the label when no short label
  expect_equal(ex_bundle$variables$variable, c("a", "b", "c", "c_b"))
  # Display in the app: rates and spreads as levels; a log series listed in
  # shown_in_levels in natural units, other log series by the growth setting.
  expect_equal(ex_bundle$variables$display, c("natural", "level", "level", "level"))
  expect_equal(export_variables(ex_series, ex_settings$spreads)$display, c("growth", "level", "level", "level"))
  # Rates have the level transform only; log series have all four.
  expect_setequal(unique(pick(ex_bundle$forecasts, variable = "a")$transform),
                  c("level", "natural", "qoq_ann", "yoy"))
  expect_identical(unique(pick(ex_bundle$forecasts, variable = "b")$transform), "level")
  expect_false(anyNA(ex_bundle$forecasts[qcols]))
  expect_equal(max(ex_bundle$forecasts$h), 4L)
  expect_s3_class(ex_bundle$forecasts$date, "Date")
  # Fixed cells within the exported quarters, in levels and natural units for log series.
  fixed <- which(!is.na(ref$conditions[1:4, ]), arr.ind = TRUE)
  expect_equal(nrow(pick(ex_bundle$conditions, object = "hard", transform = "level")), nrow(fixed))
  nat <- pick(ex_bundle$conditions, object = "hard", transform = "natural")
  expect_true(all(nat$variable == "a"))
  # Soft targets with their check.
  expect_equal(ex_bundle$targets$achieved, ref$tilted$tilt$check$achieved)
})

test_that("export settings are checked", {
  cfg <- list(covid = list(methods = list("none", "exclude_window")), forecast = list(horizon = 20),
              scenarios = list("s1", "s2"),
              export = list(run = "exclude_window", horizon = 12,
                            spreads = list(list(id = "c_b", of = list("c", "b"), label = "C minus B")),
                            shown_in_levels = list("a")))
  s <- export_settings(cfg, ex_series)
  expect_identical(s$horizon, 12L)
  expect_equal(s$spreads$a, "c")
  bad <- function(f) { x <- cfg; x$export <- f(x$export); export_settings(x, ex_series) }
  expect_error(bad(function(e) { e$run <- "precovid"; e }), "covid\\$methods")
  expect_error(bad(function(e) { e$horizon <- 24; e }), "between 1 and forecast")
  expect_equal(s$objects, c("s1", "s2"))                           # every configured scenario
  expect_error(bad(function(e) { e$objects <- list("s1"); e }), "no longer used")
  expect_error(bad(function(e) { e$spreads[[1]]$of <- list("a", "b"); e }), "rate")
  expect_error(bad(function(e) { e$spreads[[1]]$id <- "b"; e }), "differ from the variable names")
  expect_identical(s$shown_in_levels, "a")
  expect_error(bad(function(e) { e$shown_in_levels <- list("b"); e }), "must list log series")
})
