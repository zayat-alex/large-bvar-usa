# test-data-assemble.R: assembly and transformation (R/data_assemble.R) and
# the generic validation helpers (R/data_validate.R).

q <- function(...) parse_quarter(c(...))

toy_list <- list(
  a = data.frame(quarter = q("2000Q1", "2000Q2", "2000Q3", "2000Q4"), value = c(100, 101, 102, 103)),
  b = data.frame(quarter = q("2000Q2", "2000Q3", "2000Q4", "2001Q1"), value = c(5, 5.5, 6, 6.5))
)

test_that("assemble_levels ends at the last quarter shared by all series", {
  lv <- assemble_levels(toy_list, "2000Q1")
  expect_equal(lv$quarter, q("2000Q1", "2000Q2", "2000Q3", "2000Q4"))
  expect_equal(names(lv), c("quarter", "a", "b"))
  # b starts later: its first quarter stays NA, nothing is filled in.
  expect_equal(lv$b, c(NA, 5, 5.5, 6))
})

test_that("assemble_levels accepts an explicit end and rejects a missing last quarter", {
  expect_equal(nrow(assemble_levels(toy_list, "2000Q1", "2000Q3")), 3)
  expect_error(assemble_levels(toy_list, "2000Q1", "2001Q1"), "no value in the last quarter")
})

test_that("transform_to_model_units applies 100*ln to log series only", {
  lv  <- assemble_levels(toy_list, "2000Q1")
  cfg <- data.frame(short_name = c("a", "b"), transform = c("log", "raw"))
  md  <- transform_to_model_units(lv, cfg)
  expect_equal(md$a, 100 * log(c(100, 101, 102, 103)))
  expect_equal(md$b, lv$b)
  lv$a[2] <- 0
  expect_error(transform_to_model_units(lv, cfg), "non-positive values in log series a at 2000Q2")
})

test_that("check_breaks flags a planted jump and nothing else", {
  set.seed(1)
  md  <- data.frame(quarter = seq(as.Date("1990-01-01"), by = "3 months", length.out = 120),
                    x = cumsum(rnorm(120, 0.5, 0.2)))
  md$x[80:120] <- md$x[80:120] + 10          # a level shift of 50 standard deviations
  res <- check_breaks(md, data.frame(short_name = "x"), k = 4)
  expect_equal(res$data$quarter, format_quarter(md$quarter[80]))
})

test_that("backfill_before_start fills only quarters before the first observation", {
  x      <- data.frame(quarter = q("2000Q3", "2000Q4"), value = c(5, 6))
  filler <- data.frame(quarter = q("2000Q1", "2000Q2", "2000Q3", "2000Q4"), value = c(1, 2, 99, 99))
  b <- backfill_before_start(x, filler)
  expect_equal(b$series$quarter, q("2000Q1", "2000Q2", "2000Q3", "2000Q4"))
  expect_equal(b$series$value, c(1, 2, 5, 6))          # own values kept where they exist
  expect_equal(b$added$quarter, q("2000Q1", "2000Q2"))
})
