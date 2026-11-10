# test-scenarios.R: scenario files and helpers (R/scenarios.R): the rules of
# a scenario file, the conditions matrix it becomes, independence of the
# column order (the forecast must not depend on how the file is laid out),
# unit conversion and paths built from growth rates; and the scenario index
# and example files in scenarios/.

sc_vars <- c("a", "b", "c")

#' Write lines to a temporary CSV and return its path.
sc_csv <- function(...) {
  path <- tempfile(fileext = ".csv")
  writeLines(c(...), path)
  path
}

#' Read a scenario with the test settings: data end in 2026Q2, H = 6.
sc_read <- function(path, variables = sc_vars) read_scenario(path, variables, "2026Q2", 6)

#' Drop the quarter row names (the engine's conditions have column names only).
unname_rows <- function(m) {
  rownames(m) <- NULL
  m
}

test_that("quarters are counted and listed", {
  expect_equal(quarter_index("2027Q1") - quarter_index("2026Q2"), 3L)
  expect_true(all(is.na(quarter_index(c("2027Q5", "2027-Q1", "2027:Q1", "27Q1", "", NA)))))
  expect_equal(forecast_quarters("2026Q3", 3), c("2026Q4", "2027Q1", "2027Q2"))
  expect_error(forecast_quarters("2026-06", 3), "not a quarter")
})

test_that("a valid file becomes the conditions matrix, blank cells free", {
  path <- sc_csv("quarter,b,a",
                 "2026Q3,1.5,",
                 "2027Q1,,2",
                 "2027Q4, 0.25 ,NA")
  m <- sc_read(path)
  expect_equal(dim(m), c(6L, 3L))
  expect_equal(rownames(m), forecast_quarters("2026Q2", 6))
  expect_equal(colnames(m), sc_vars)
  expected <- matrix(NA_real_, 6, 3, dimnames = dimnames(m))
  expected["2026Q3", "b"] <- 1.5
  expected["2027Q1", "a"] <- 2
  expected["2027Q4", "b"] <- 0.25
  expect_identical(m, expected)
})

test_that("each rule of a scenario file is enforced", {
  expect_error(sc_read(sc_csv("date,a", "2026Q3,1")), "no `quarter` column")
  expect_error(sc_read(sc_csv("quarter,a,zz", "2026Q3,1,2")), "unknown column.*zz")
  expect_error(sc_read(sc_csv("quarter,a,a", "2026Q3,1,2")), "repeated column.*a")
  expect_error(sc_read(sc_csv("quarter,a", "2026Q3,1", "2026-Q4,1")), "not a quarter label.*2026-Q4")
  expect_error(sc_read(sc_csv("quarter,a", "2026Q3,1", "2026Q3,2")), "repeated quarter.*2026Q3")
  # 2026Q2 is observed; 2028Q1 is quarter 7 with H = 6.
  expect_error(sc_read(sc_csv("quarter,a", "2026Q2,1")), "outside the forecast range 2026Q3 to 2027Q4.*2026Q2")
  expect_error(sc_read(sc_csv("quarter,a", "2028Q1,1")), "outside the forecast range.*2028Q1")
  expect_error(sc_read(sc_csv("quarter,a", "2026Q3,1.2.3")), "non-numeric.*a: 1.2.3")
  expect_error(sc_read(sc_csv("quarter,a", "2026Q3,Inf")), "non-numeric")
  expect_error(sc_read(sc_csv("quarter,a,b,c", "2026Q3,1,2,3")), "every variable is fixed.*2026Q3")
  # Several problems are reported together, and the message names the file.
  path <- sc_csv("quarter,a,zz", "2030Q1,x,1")
  err <- tryCatch(sc_read(path), error = conditionMessage)
  expect_match(err, basename(path), fixed = TRUE)
  expect_match(err, "unknown column")
  expect_match(err, "outside the forecast range")
  expect_match(err, "non-numeric")
})

test_that("results do not depend on the column order of the file", {
  # The same conditions written with the columns in two orders give the
  # same matrix and, on a known VAR(2) (as in test-forecast-conditional.R),
  # the same paths.
  path1 <- sc_csv("quarter,a,b,c", "2026Q3,2,,", "2026Q4,,0,1", "2027Q2,1,,")
  path2 <- sc_csv("c,quarter,b,a", "1,2026Q4,0,", ",2027Q2,,1", ",2026Q3,,2")
  m1 <- sc_read(path1)
  m2 <- sc_read(path2)
  expect_identical(m1, m2)
  B1 <- matrix(c(0.6, 0.1, 0.0, 0.2, 0.5, 0.1, 0.0, 0.1, 0.7), 3)
  B2 <- matrix(c(0.2, 0.0, 0.05, -0.1, 0.2, 0.0, 0.0, 0.0, 0.1), 3)
  B <- rbind(c(0.5, -0.2, 0.1), t(B1), t(B2))
  Sigma <- matrix(c(1, 0.4, 0.2, 0.4, 0.5, 0.1, 0.2, 0.1, 0.8), 3)
  hist <- matrix(c(1.0, 0.2, -0.5, 0.8, 1.1, 0.3, 0.0, -0.4, 0.6, 1.2, 0.9, 0.1, 0.4, 0.7, -0.2),
                 5, 3, dimnames = list(NULL, sc_vars))
  beta <- array(B, c(1, dim(B)))
  sigma <- array(Sigma, c(1, dim(Sigma)))
  for (method in c("gaussian", "kalman")) {
    f1 <- forecast_conditional(beta, sigma, hist, 2, 6, unname_rows(m1), seed = 7, method = method)
    f2 <- forecast_conditional(beta, sigma, hist, 2, 6, unname_rows(m2), seed = 7, method = method)
    expect_identical(f1, f2)
    expect_equal(f1[1, 1, 1], 2)
  }
})

test_that("values convert to model units", {
  series <- read_series_config(file.path(project_root(), "config", "series.csv"))
  expect_equal(to_model_units("rgdp", 24500, series), 100 * log(24500))
  expect_equal(to_model_units("urate", 4.3, series), 4.3)
  expect_error(to_model_units("wti", c(80, 0), series), "positive")
  expect_error(to_model_units("nope", 1, series), "unknown variable")
})

test_that("paths from growth rates reproduce the growth rates exactly", {
  series <- read_series_config(file.path(project_root(), "config", "series.csv"))
  hist <- matrix(100 * log(c(100, 101, 102, 103, 104)), 5, 1,
                 dimnames = list(c("2025Q2", "2025Q3", "2025Q4", "2026Q1", "2026Q2"), "rgdp"))
  # Quarter on quarter at an annual rate: (x_t / x_{t-1})^4 - 1 = g.
  g <- c(2, 3, -1, 0.5)
  y <- path_from_growth("rgdp", g, "qoq_ann", hist, series)
  expect_equal(names(y), c("2026Q3", "2026Q4", "2027Q1", "2027Q2"))
  x <- exp(c(hist[5, 1], y) / 100)
  expect_equal(unname(100 * ((x[-1] / x[-5])^4 - 1)), g)
  # Year on year: x_t / x_{t-4} - 1 = g, using the observed quarters as base.
  y <- path_from_growth("rgdp", c(NA, 3.4, NA, NA, NA, 2.2), "yoy", hist, series)
  expect_true(all(is.na(y[c(1, 3, 4, 5)])))
  expect_equal(100 * (exp((y[2] - hist[3, 1]) / 100) - 1), c("2026Q4" = 3.4))   # 2026Q4 over 2025Q4
  expect_equal(100 * (exp((y[6] - y[2]) / 100) - 1), c("2027Q4" = 2.2))         # 2027Q4 over 2026Q4
  # A growth rate whose base quarter is free cannot be turned into a level.
  expect_error(path_from_growth("rgdp", c(NA, 1), "qoq_ann", hist, series), "free")
  expect_error(path_from_growth("urate", 1, "yoy", hist, series), "not a log series")
  expect_error(path_from_growth("rgdp", 1, "annual", hist, series), "type must be")
})

test_that("the scenario index and the example files are well formed", {
  root <- project_root()
  series <- read_series_config(file.path(root, "config", "series.csv"))
  horizon <- yaml::read_yaml(file.path(root, "config", "run.yml"))$forecast$horizon
  idx <- read_scenario_index(file.path(root, "scenarios", "scenarios.yml"))
  for (id in idx$id) {
    path <- file.path(root, "scenarios", paste0(id, ".csv"))
    expect_true(file.exists(path), info = id)
    # Check the file against the quarter before its first row, so the test
    # does not depend on the data in outputs/.
    q <- utils::read.csv(path, colClasses = "character")$quarter
    k <- min(quarter_index(q)) - 1L
    before <- paste0(k %/% 4L, "Q", k %% 4L + 1L)
    soft <- scenario_soft(idx, id)
    m <- read_scenario(path, series$short_name, before, horizon, soft)
    tg <- soft_targets(path, series$short_name, before, horizon, soft, series)
    expect_gt(sum(!is.na(m)) + nrow(tg), 0)
  }
  # Index errors.
  bad <- tempfile(fileext = ".yml")
  writeLines(c("scenarios:", "  - {id: a, label: A, description: one}", "  - {id: a, label: B, description: two}"), bad)
  expect_error(read_scenario_index(bad), "repeated")
  writeLines(c("scenarios:", "  - {id: Bad-Id, label: A, description: x}"), bad)
  expect_error(read_scenario_index(bad), "lower-case")
  writeLines(c("scenarios:", "  - {id: a, label: A}"), bad)
  expect_error(read_scenario_index(bad), "lacks: description")
  writeLines(c("scenarios:", "  - {id: a, label: A, description: x, colour: red}"), bad)
  expect_error(read_scenario_index(bad), "unknown field\\(s\\): colour")
  writeLines(c("scenarios:", "  - {id: a, label: A, description: x, soft: {corepce: growth}}"), bad)
  expect_error(read_scenario_index(bad), "entry 1.*measure must be")
  writeLines("scenarios: []", bad)
  expect_equal(nrow(read_scenario_index(bad)), 0L)
  # Optional short label, the label when absent.
  writeLines(c("scenarios:", "  - {id: a, label: Long A, short_label: A, description: one}",
               "  - {id: b, label: Long B, description: two}"), bad)
  expect_equal(read_scenario_index(bad)$short_label, c("A", "Long B"))
})

test_that("a scenario's base is a tilted scenario without a base, and the scenario has no soft targets", {
  bad <- tempfile(fileext = ".yml")
  on.exit(unlink(bad))
  tilted <- "  - {id: t, label: T, description: x, soft: {a: yoy}}"
  writeLines(c("scenarios:", tilted, "  - {id: h, label: H, description: x, base: t}",
               "  - {id: u, label: U, description: x}"), bad)
  idx <- read_scenario_index(bad)
  expect_equal(idx$base, c(NA, "t", NA))
  expect_equal(scenario_base(idx, "h"), "t")
  expect_true(is.na(scenario_base(idx, "u")))
  # Each rule.
  writeLines(c("scenarios:", "  - {id: h, label: H, description: x, base: h}"), bad)
  expect_error(read_scenario_index(bad), "its own base")
  writeLines(c("scenarios:", "  - {id: h, label: H, description: x, base: nope}"), bad)
  expect_error(read_scenario_index(bad), "base nope is not in the index")
  writeLines(c("scenarios:", "  - {id: u, label: U, description: x}", "  - {id: h, label: H, description: x, base: u}"), bad)
  expect_error(read_scenario_index(bad), "base u is not tilted")
  writeLines(c("scenarios:", tilted, "  - {id: t2, label: T2, description: x, soft: {a: yoy}, base: t}",
               "  - {id: h, label: H, description: x, base: t2}"), bad)
  expect_error(read_scenario_index(bad), "base t2 has a base of its own")
  expect_error(read_scenario_index(bad), "t2: a scenario with a base cannot have soft conditions")
  writeLines(c("scenarios:", "  - {id: h, label: H, description: x, base: [t, u]}"), bad)
  expect_error(read_scenario_index(bad), "base must be one scenario id")
})

test_that("soft conditions are parsed, checked and read as targets", {
  # Short and long forms; the statistic defaults to the median.
  s <- parse_soft(list(a = "yoy", c = list(measure = "level", statistic = "mean")))
  expect_equal(s$variable, c("a", "c"))
  expect_equal(s$measure, c("yoy", "level"))
  expect_equal(s$statistic, c("median", "mean"))
  expect_equal(nrow(parse_soft(NULL)), 0L)
  expect_equal(parse_soft(list(a = "annual_avg"))$measure, "annual_avg")
  expect_error(parse_soft(list(a = "growth")), "measure must be")
  expect_error(parse_soft(list(a = list(measure = "yoy", statistic = "mode"))), "statistic must be")
  expect_error(parse_soft(list(a = list(measure = "yoy", band = 1))), "unknown field")
  expect_error(parse_soft("a"), "must map")

  # The index keeps each scenario's soft table.
  idx_path <- tempfile(fileext = ".yml")
  writeLines(c("scenarios:", "  - {id: s1, label: S, description: x, soft: {a: yoy}}",
               "  - {id: s2, label: T, description: two}"), idx_path)
  idx <- read_scenario_index(idx_path)
  expect_equal(scenario_soft(idx, "s1")$measure, "yoy")
  expect_equal(nrow(scenario_soft(idx, "s2")), 0L)
  expect_error(scenario_soft(idx, "s3"), "not in the index")

  # A file with a hard column (b) and a soft column (a): the conditions keep
  # only the hard cells; the targets come from the soft column.
  series <- data.frame(short_name = sc_vars, transform = c("log", "log", "raw"), stringsAsFactors = FALSE)
  path <- sc_csv("quarter,a,b", "2026Q3,,1.5", "2027Q2,2.1,", "2027Q4,2.0,0.5")
  soft <- parse_soft(list(a = "yoy"))
  m <- read_scenario(path, sc_vars, "2026Q2", 6, soft)
  expect_true(all(is.na(m[, "a"])))
  expect_equal(unname(m[c(1, 6), "b"]), c(1.5, 0.5))
  tg <- soft_targets(path, sc_vars, "2026Q2", 6, soft, series)
  expect_equal(tg$quarter, c("2027Q2", "2027Q4"))
  expect_equal(tg$h, c(4L, 6L))
  expect_equal(tg$target, c(2.1, 2.0))
  expect_equal(unique(tg$statistic), "median")
  expect_equal(nrow(soft_targets(path, sc_vars, "2026Q2", 6, NULL, series)), 0L)
  # Without the soft argument the same column would be read as hard levels.
  expect_equal(unname(sc_read(path)[4, "a"]), 2.1)

  # Errors: growth of a rate; a soft column missing or empty.
  expect_error(soft_targets(sc_csv("quarter,c", "2026Q3,1"), sc_vars, "2026Q2", 6,
                            parse_soft(list(c = "yoy")), series), "only to log series")
  expect_error(read_scenario(sc_csv("quarter,b", "2026Q3,1"), sc_vars, "2026Q2", 6, soft),
               "without a column: a")
  expect_error(soft_targets(sc_csv("quarter,a,b", "2026Q3,,1"), sc_vars, "2026Q2", 6, soft, series),
               "has no target")
  # A quarter with every variable filled is allowed when one of them is soft.
  full <- sc_csv("quarter,a,b,c", "2026Q3,2,1,1")
  expect_error(sc_read(full), "every variable is fixed")
  expect_equal(sum(!is.na(read_scenario(full, sc_vars, "2026Q2", 6, soft))), 2L)
})
