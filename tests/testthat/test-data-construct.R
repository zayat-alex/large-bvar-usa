# test-data-construct.R: constructed series (R/data_assemble.R).
#
# The realized volatility formula of the docs/data.md (3.4) is checked on inputs
# small enough to compute by hand.

test_that("realized volatility matches a hand calculation", {
  # Closes: 100 on the last day of 2024Q4, then 101 and 100 in 2025Q1.
  # Returns in 2025Q1: r1 = ln(1.01), r2 = ln(100/101) = -r1, so mean(r^2) = r1^2
  # and rvol = 100 * sqrt(252) * ln(1.01).
  daily <- data.frame(date  = as.Date(c("2024-12-31", "2025-01-02", "2025-01-03")),
                      close = c(100, 101, 100))
  rv <- construct_realized_volatility(daily, complete_through = as.Date("2025-04-15"))
  expect_equal(rv$quarter, as.Date("2025-01-01"))
  expect_equal(rv$value, 100 * sqrt(252) * log(1.01))
  expect_equal(rv$n_days, 2L)
})

test_that("realized volatility uses the mean, not the sum, of squared returns", {
  # Same return size every day: rvol does not depend on the number of days.
  make <- function(n) {
    data.frame(date  = as.Date("2024-12-31") + 0:n,
               close = 100 * 1.01^(0:n))
  }
  a <- construct_realized_volatility(make(10), as.Date("2025-12-31"))
  b <- construct_realized_volatility(make(40), as.Date("2025-12-31"))
  expect_equal(a$value[1], b$value[1])
})

test_that("realized volatility keeps whole quarters only", {
  daily <- data.frame(
    date  = as.Date(c("2025-02-10", "2025-02-11", "2025-05-02", "2025-05-05", "2025-08-01")),
    close = c(100, 101, 102, 101, 103)
  )
  # 2025Q1: data start inside the quarter -> dropped.
  # 2025Q2: complete -> kept. 2025Q3: not ended by complete_through -> dropped.
  rv <- construct_realized_volatility(daily, complete_through = as.Date("2025-08-05"))
  expect_equal(rv$quarter, as.Date("2025-04-01"))
  # Once the quarter has ended, 2025Q3 is kept.
  rv2 <- construct_realized_volatility(daily, complete_through = as.Date("2025-09-30"))
  expect_equal(rv2$quarter, as.Date(c("2025-04-01", "2025-07-01")))
  expect_error(construct_realized_volatility(transform(daily, close = -close), Sys.Date()),
               "positive")
})
