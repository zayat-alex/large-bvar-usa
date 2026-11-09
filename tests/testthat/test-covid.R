# test-covid.R: the COVID treatment stage (R/covid.R).
#
# A small hand-made table with known answers checks each rule; the real
# model_data, when the local data build exists, checks the sample sizes the
# estimation will use.

q <- function(...) parse_quarter(c(...))

# Toy model data: 16 quarters (2018Q1 to 2021Q4), two variables.
toy_quarters <- seq(q("2018Q1"), by = "3 months", length.out = 16)
toy_data <- data.frame(quarter = toy_quarters, x = 1:16, y = 101:116)

toy_cfg <- function(method, cutoff = "2019Q4", lags = 2) {
  list(model = list(lags = lags),
       covid = list(method = method, window = c("2020Q1", "2021Q4"), cutoff = cutoff))
}

test_that("method none returns the full sample for estimation and forecasting", {
  out <- covid_treatment(toy_data, toy_cfg("none"))
  expect_named(out, c("est_data", "origin_data", "est_weights", "notes"))
  expect_identical(out$est_data, toy_data)
  expect_identical(out$origin_data, toy_data)
  expect_null(out$est_weights)
  expect_equal(nrow(out$notes), 0)
  expect_named(out$notes, c("quarter", "quarter_label", "method", "action", "reason"))
})

test_that("method precovid_params estimates up to the cutoff and forecasts from the latest data", {
  out <- covid_treatment(toy_data, toy_cfg("precovid_params"))
  # Estimation sample: 2018Q1 to 2019Q4, 8 quarters.
  expect_equal(nrow(out$est_data), 8)
  expect_equal(max(out$est_data$quarter), q("2019Q4"))
  expect_equal(out$est_data$x, 1:8)
  # Forecast history: the full actual sample.
  expect_identical(out$origin_data, toy_data)
  # Notes: the 8 dropped quarters, 2020Q1 to 2021Q4.
  expect_equal(out$notes$quarter, toy_quarters[9:16])
  expect_equal(out$notes$quarter_label[c(1, 8)], c("2020Q1", "2021Q4"))
  expect_true(all(out$notes$method == "precovid_params"))
  expect_true(all(out$notes$action == "dropped from est_data"))
  expect_true(all(out$notes$reason == "after the estimation cutoff 2019Q4"))
})

test_that("origin_data is the same under both methods", {
  a <- covid_treatment(toy_data, toy_cfg("none"))
  b <- covid_treatment(toy_data, toy_cfg("precovid_params"))
  expect_identical(a$origin_data, b$origin_data)
})

test_that("invalid settings are rejected", {
  expect_error(covid_treatment(toy_data, toy_cfg("drop_2020")), "unknown covid\\$method")
  expect_error(covid_treatment(toy_data, toy_cfg("scaling")), "not implemented")
  expect_error(covid_treatment(toy_data, toy_cfg("interpolate_window")), "not implemented")
  expect_error(covid_treatment(toy_data, toy_cfg("precovid_params", cutoff = "2019-12")),
               "not a quarter label")
  expect_error(covid_treatment(toy_data, toy_cfg("precovid_params", cutoff = "2017Q4")),
               "before the sample start")
  expect_error(covid_treatment(toy_data, toy_cfg("precovid_params", cutoff = "2021Q4")),
               "not before the last quarter")
  # Cutoff 2018Q2 keeps 2 quarters, not more than 2 lags.
  expect_error(covid_treatment(toy_data, toy_cfg("precovid_params", cutoff = "2018Q2")),
               "not more than the 2 lags")
  expect_error(covid_treatment(toy_data, toy_cfg("none", lags = 0)), "positive whole number")
})

test_that("check_covid_output catches broken results", {
  good <- covid_none(toy_data)
  expect_true(check_covid_output(good, toy_data, lags = 2))

  bad <- good
  bad$est_data$x[3] <- NA
  expect_error(check_covid_output(bad, toy_data, 2), "missing values")

  bad <- good
  bad$origin_data <- bad$origin_data[, c("quarter", "y", "x")]
  expect_error(check_covid_output(bad, toy_data, 2), "same order")

  bad <- good
  bad$est_data <- bad$est_data[-5, ]            # a gap in the quarters
  expect_error(check_covid_output(bad, toy_data, 2), "consecutive quarters")

  bad <- good
  bad$origin_data <- bad$origin_data[1:10, ]    # does not reach the last quarter
  expect_error(check_covid_output(bad, toy_data, 2), "must end at the last quarter")

  bad <- good
  bad$est_weights <- rep(2, 16)
  expect_error(check_covid_output(bad, toy_data, 2), "est_weights must be 0 or 1")

  bad <- good
  bad$est_weights <- rep(1, 15)
  expect_error(check_covid_output(bad, toy_data, 2), "est_weights must be 0 or 1")

  bad <- good[c("est_data", "origin_data", "notes")]
  expect_error(check_covid_output(bad, toy_data, 2), "must have elements")
})

test_that("on the real data the sample sizes are as expected", {
  cfg  <- yaml::read_yaml(file.path(project_root(), "config", "run.yml"))
  path <- file.path(project_root(), cfg$data$output_dir, "data.rds")
  skip_if_not(file.exists(path), "data build not present (outputs are not in git)")
  model_data <- readRDS(path)$model_data
  last <- max(model_data$quarter)
  n_all <- nrow(model_data)

  cfg$covid$method <- "none"
  out_none <- covid_treatment(model_data, cfg)
  expect_equal(nrow(out_none$est_data), n_all)

  cfg$covid$method <- "precovid_params"
  cfg$covid$cutoff <- "2019Q4"
  out_pre <- covid_treatment(model_data, cfg)
  # 1973Q1 to 2019Q4 is 47 years of 4 quarters.
  expect_equal(nrow(out_pre$est_data), 188)
  expect_equal(max(out_pre$est_data$quarter), q("2019Q4"))
  # Dropped: every quarter from 2020Q1 to the last quarter.
  expect_equal(out_pre$notes$quarter, seq(q("2020Q1"), last, by = "3 months"))
  expect_identical(out_pre$origin_data, model_data)

  cfg$covid$method <- "exclude_window"
  cfg$covid$window <- c("2020Q1", "2021Q4")
  out_ex <- covid_treatment(model_data, cfg)
  expect_identical(out_ex$est_data, model_data)
  # Eight quarters' equations left out, the rest used.
  expect_equal(sum(out_ex$est_weights == 0), 8)
  expect_equal(out_ex$notes$quarter, seq(q("2020Q1"), q("2021Q4"), by = "3 months"))
})

test_that("method exclude_window keeps every quarter but marks the window's equations", {
  cfg <- toy_cfg("exclude_window")
  cfg$covid$window <- c("2019Q1", "2019Q4")
  out <- covid_treatment(toy_data, cfg)
  expect_identical(out$est_data, toy_data)
  expect_identical(out$origin_data, toy_data)
  expect_equal(out$est_weights, as.numeric(!(toy_quarters >= q("2019Q1") & toy_quarters <= q("2019Q4"))))
  expect_equal(out$notes$quarter_label, c("2019Q1", "2019Q2", "2019Q3", "2019Q4"))
  expect_equal(covid_keep(out), out$est_weights == 1)
  expect_null(covid_keep(covid_treatment(toy_data, toy_cfg("none"))))
})

test_that("exclude_window rejects a window outside the usable sample", {
  cfg <- toy_cfg("exclude_window")
  cfg$covid$window <- c("2018Q2", "2019Q1")                # inside the initial lags
  expect_error(covid_treatment(toy_data, cfg), "initial lags")
  cfg$covid$window <- c("2020Q1", "2021Q4")                # reaches the last quarter
  expect_error(covid_treatment(toy_data, cfg), "must end before")
  cfg$covid$window <- c("2020Q1", "2019Q1")
  expect_error(covid_treatment(toy_data, cfg), "start before it ends")
  cfg$covid$window <- "2020Q1"
  expect_error(covid_treatment(toy_data, cfg), "two quarters")
})
