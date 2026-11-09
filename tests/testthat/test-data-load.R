# test-data-load.R: loaders and quarter handling (R/data_load.R).
#
# Small hand-made inputs with known answers, so each rule of the docs/data.md
# (section 1) is checked on its own: quarter conventions, QD layout, the
# whole-quarter rule, and reading ALFRED vintages as of a date.

test_that("quarter helpers follow the first-day-of-quarter convention", {
  d <- as.Date(c("1973-01-01", "1973-02-15", "1973-03-31", "1973-04-01", "2025-12-31"))
  expect_equal(quarter_start(d),
               as.Date(c("1973-01-01", "1973-01-01", "1973-01-01", "1973-04-01", "2025-10-01")))
  expect_equal(quarter_end(as.Date(c("1973-01-01", "1973-05-10", "2024-07-01", "2024-10-01"))),
               as.Date(c("1973-03-31", "1973-06-30", "2024-09-30", "2024-12-31")))
  expect_equal(format_quarter(as.Date(c("1973-01-01", "2026-04-01"))), c("1973Q1", "2026Q2"))
  # Edge cases: no quarters gives no labels; a missing date gives NA.
  expect_identical(format_quarter(as.Date(character(0))), character(0))
  expect_identical(format_quarter(as.Date(c("1973-01-01", NA))), c("1973Q1", NA))
  expect_equal(parse_quarter(c("1973Q1", "2019:Q4")), as.Date(c("1973-01-01", "2019-10-01")))
  expect_error(parse_quarter("1973Q5"), "not a quarter label")
})

test_that("QD date stamps (first day of last month) map to quarters", {
  expect_equal(qd_date_to_quarter(c("3/1/1973", "6/1/1973", "12/1/2025")),
               as.Date(c("1973-01-01", "1973-04-01", "2025-10-01")))
  expect_error(qd_date_to_quarter("2/1/1973"), "last month of a quarter")
  expect_error(qd_date_to_quarter("3/15/1973"), "last month of a quarter")
})

test_that("read_qd_vintage drops the factors and transform rows", {
  path <- tempfile(fileext = ".csv")
  writeLines(c(
    "sasdate,GDPC1,S&P 500,UNRATE",
    "factors,0,1,1",
    "transform,5,5,2",
    "3/1/1973,100.5,110.2,4.9",
    "6/1/1973,101.0,,4.8",
    ","
  ), path)
  qd <- read_qd_vintage(path, c("GDPC1", "S&P 500"))
  expect_equal(qd$quarter, as.Date(c("1973-01-01", "1973-04-01")))
  expect_equal(qd$GDPC1, c(100.5, 101.0))
  expect_equal(qd[["S&P 500"]], c(110.2, NA))         # empty cell -> NA
  expect_false("UNRATE" %in% names(qd))
  expect_error(read_qd_vintage(path, "PAYEMS"), "columns not in")
})

test_that("monthly_to_quarterly averages whole quarters only", {
  m <- data.frame(
    date  = seq(as.Date("2025-01-01"), by = "month", length.out = 8),
    value = c(1, 2, 3,  4, NA, 6,  7, 8)
  )
  q <- monthly_to_quarterly(m)
  # 2025Q1 = mean(1, 2, 3); 2025Q2 has a missing month; 2025Q3 has two months.
  expect_equal(q$quarter, as.Date("2025-01-01"))
  expect_equal(q$value, 2)
  expect_error(monthly_to_quarterly(rbind(m, m[1, ])), "duplicated months")
})

test_that("alfred_as_of returns the values published on the date", {
  v <- data.frame(
    series_id      = "X",
    date           = as.Date(c("2019-01-01", "2019-01-01", "2019-01-01", "2019-04-01")),
    value          = c(10, 11, 12, 20),
    realtime_start = as.Date(c("2018-05-01", "2020-07-30", "2021-07-29", "2019-07-26")),
    realtime_end   = as.Date(c("2020-07-29", "2021-07-28", "2026-10-17", "2026-10-17"))
  )
  # Before the 2019Q2 value was published, only 2019Q1 exists.
  expect_equal(alfred_as_of(v, "2019-06-01")$value, 10)
  # Inclusive ends: on 2020-07-29 the first value still holds; a day later the second.
  expect_equal(alfred_as_of(v, "2020-07-29")$value, c(10, 20))
  expect_equal(alfred_as_of(v, "2020-07-30")$value, c(11, 20))
  expect_equal(alfred_as_of(v, "2026-09-30")$value, c(12, 20))
  expect_error(alfred_as_of(v, "2027-01-01"), "outside the downloaded real-time window")
  overlap <- v
  overlap$realtime_end[1] <- as.Date("2020-12-31")      # two values valid on one date
  expect_error(alfred_as_of(overlap, "2020-09-01"), "more than one value")
})

test_that("read_worldbank_index finds the non-energy column in the stored release", {
  cfg  <- yaml::read_yaml(file.path(project_root(), "config", "run.yml"))
  path <- file.path(project_root(), cfg$data$raw_dir, cfg$data$worldbank_file)
  skip_if_not(file.exists(path), "World Bank raw file not present (raw data are not in git)")
  wb <- read_worldbank_index(path)
  expect_equal(min(wb$date), as.Date("1960-01-01"))
  expect_true(all(!is.na(wb$value)))
  expect_true(all(diff(wb$date) <= 31))                 # consecutive months, no gaps
})

test_that("read_gsw_daily skips the notes and daily_to_quarterly keeps whole quarters", {
  path <- tempfile(fileext = ".csv")
  writeLines(c(
    "\"Note: staff research product\"",
    "",
    "Series,Compounding Convention,Mnemonic(s)",
    "",
    "Date,SVENPY02,SVENY02",
    "2024-12-31,4.0,3.9",
    "2025-01-02,4.2,4.1",
    "2025-02-03,NA,4.0",
    "2025-03-31,4.4,4.3",
    "2025-04-01,4.6,4.5"
  ), path)
  g <- read_gsw_daily(path, "SVENPY02")
  expect_equal(g$value, c(4.0, 4.2, 4.4, 4.6))          # NA day dropped
  expect_error(read_gsw_daily(path, "SVENPY10"), "not in")
  # 2024Q4: data start inside the quarter -> dropped. 2025Q1 complete -> mean(4.2, 4.4).
  # 2025Q2: not ended by complete_through -> dropped.
  qtr <- daily_to_quarterly(g, as.Date("2025-04-01"))
  expect_equal(qtr$quarter, as.Date("2025-01-01"))
  expect_equal(qtr$value, 4.3)
})
