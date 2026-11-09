# test-forecast-review.R: the measures of the forecast reasonableness review
# (R/forecast_review.R). The review of the production runs is
# scripts/check_forecast_review.R; read_spf_survey() is exercised there,
# since the SPF files are raw data outside the repository.

test_that("qoq_annualized_exact gives the annualized percent change", {
  # 100 x log levels rising by 0.5 per quarter: exact annualized growth is
  # 100 (exp(0.02) - 1), slightly above the log approximation of 2.
  history <- cbind(g = c(100, 100.5), r = c(4, 4))
  paths <- array(NA_real_, c(2, 3, 2), dimnames = list(NULL, NULL, c("g", "r")))
  for (i in 1:2) paths[i, , "g"] <- 100 + 0.5 * (2:4)
  paths[, , "r"] <- 5
  q <- qoq_annualized_exact(paths, history, c(TRUE, FALSE))
  expect_equal(q[, , "g"], matrix(100 * (exp(0.02) - 1), 2, 3))
  expect_true(all(is.na(q[, , "r"])))
})

test_that("compare_spf places the SPF median in the model's distribution", {
  m <- array(rep(1:100, 2), c(100, 2, 1), dimnames = list(NULL, NULL, "x"))
  spf <- data.frame(variable = "x", horizon = 1:2, spf = c(25, 101))
  out <- compare_spf(m, spf)
  expect_equal(out$spf_percentile, c(25, 100))
  expect_equal(out$q50, rep(50.5, 2))
  expect_equal(out$spf_minus_median, c(25 - 50.5, 101 - 50.5))
})

test_that("historical_ranges and range_check use YoY growth for log series and levels for rates", {
  y <- cbind(g = 100 + 0.5 * (1:40), r = rep(c(3, 5), 20))
  rownames(y) <- format_quarter(seq(as.Date("2010-01-01"), by = "quarter", length.out = 40))
  rg <- historical_ranges(y, c(TRUE, FALSE), end = rownames(y)[40])
  expect_equal(rg$mean, c(2, 4))
  expect_equal(c(rg$low[1], rg$high[1]), c(2, 2))
  annual <- array(NA_real_, c(4, 1, 2), dimnames = list(NULL, NULL, c("g", "r")))
  annual[, 1, "g"] <- c(1, 2, 2, 3)
  annual[, 1, "r"] <- c(2, 4, 4, 6)
  rc <- range_check(annual, rg, 1)
  expect_equal(rc$share_below, c(0.25, 0.25))
  expect_equal(rc$share_above, c(0.25, 0.25))
  expect_equal(rc$median_outside, c(FALSE, FALSE))
})

test_that("negative_share, measure_quantiles and seed_difference", {
  p <- array(c(-1, 1, 2, 3, -2, -1, 1, 2), c(4, 2, 1), dimnames = list(NULL, NULL, "y"))
  expect_equal(negative_share(p, "y")$share_negative, c(0.25, 0.5))
  q <- measure_quantiles(p, "y", 1)
  expect_equal(q$q50, median(c(-1, 1, 2, 3)))
  # Weighted: all weight on the last path gives its value at every quantile.
  qw <- measure_quantiles(p, "y", 2, weights = c(0, 0, 0, 1))
  expect_equal(unlist(qw[, c("q05", "q50", "q95")]), c(q05 = 2, q50 = 2, q95 = 2))
  d <- seed_difference(p, p + 0.5, "y")
  expect_equal(d$diff, c(0.5, 0.5))
})
