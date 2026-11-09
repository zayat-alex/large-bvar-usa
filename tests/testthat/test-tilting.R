# test-tilting.R: entropic tilting (R/tilting.R), the properties of
# docs/methodology.md section 9.4 that tilting must satisfy: targets met,
# weights valid, equal weights when nothing has to change, the closed form
# for a Gaussian sample, minimal Kullback-Leibler divergence, infeasible
# targets rejected; growth quantities by hand; weighted quantiles.

#' Paths of one variable "a" (log series in model units) over H = 6 quarters,
#' a random walk with drift, plus a second variable "r" (a rate).
tl_paths <- function(J = 4000, seed = 1) {
  set.seed(seed)
  a <- 1000 + t(apply(matrix(rnorm(J * 6, 0.5, 1), J, 6), 1, cumsum))
  r <- matrix(rnorm(J * 6, 3, 1), J, 6)
  arr <- array(NA_real_, c(J, 6, 2), dimnames = list(NULL, paste0("h", 1:6), c("a", "r")))
  arr[, , "a"] <- a
  arr[, , "r"] <- r
  arr
}
tl_hist <- matrix(c(997, 998, 998.5, 999.2, 1000, 3, 3.1, 2.9, 3, 3), 5, 2,
                  dimnames = list(c("2025Q2", "2025Q3", "2025Q4", "2026Q1", "2026Q2"), c("a", "r")))
tl_target <- function(variable, h, measure, statistic, target) {
  data.frame(variable = variable, h = h, measure = measure, statistic = statistic, target = target,
             stringsAsFactors = FALSE)
}

test_that("growth quantities match a hand calculation, with base quarters from the history", {
  p <- tl_paths(3)
  # yoy at h = 2: base is the history row T + 2 - 4 = T - 2 (2025Q4).
  expect_equal(tilt_quantity(p, "a", 2, "yoy", tl_hist), 100 * (exp((p[, 2, "a"] - 998.5) / 100) - 1))
  # yoy at h = 6: base is forecast quarter 2.
  expect_equal(tilt_quantity(p, "a", 6, "yoy", tl_hist), 100 * (exp((p[, 6, "a"] - p[, 2, "a"]) / 100) - 1))
  # qoq_ann at h = 1: base is the last observed quarter; exact annualized rate.
  x1 <- exp(p[, 1, "a"] / 100); x0 <- exp(1000 / 100)
  expect_equal(tilt_quantity(p, "a", 1, "qoq_ann", tl_hist), 100 * ((x1 / x0)^4 - 1))
  expect_equal(tilt_quantity(p, "r", 3, "level", tl_hist), p[, 3, "r"])
})

test_that("annual-average growth matches a hand calculation, with base quarters from the history", {
  p <- tl_paths(3)
  # h = 6: current year = forecast quarters 3..6; previous year = forecast
  # quarters 1, 2 and the history rows T (2026Q2) and T - 1 (2026Q1).
  x <- function(y) exp(y / 100)
  cur  <- x(p[, 3, "a"]) + x(p[, 4, "a"]) + x(p[, 5, "a"]) + x(p[, 6, "a"])
  prev <- x(999.2) + x(1000) + x(p[, 1, "a"]) + x(p[, 2, "a"])
  expect_equal(tilt_quantity(p, "a", 6, "annual_avg", tl_hist), 100 * (cur / prev - 1))
  # h = 3 reaches back to T - 4 (2025Q2), the first history row.
  expect_length(tilt_quantity(p, "a", 3, "annual_avg", tl_hist), 3)
  # h = 2 would need T - 5, which this history does not have.
  expect_error(tilt_quantity(p, "a", 2, "annual_avg", tl_hist), "needs 6 observed quarters")
})

test_that("with constant growth, annual-average growth equals four-quarter growth", {
  # y grows by 0.6 log points a quarter, history and forecast alike, so every
  # four-quarter average is the previous one times exp(2.4 / 100).
  hist <- matrix(1000 + 0.6 * (1:8), 8, 1, dimnames = list(paste0("q", 1:8), "a"))
  # Two identical paths: with a single path, R drops the path dimension.
  p <- array(rep(1000 + 0.6 * (9:14), each = 2), c(2, 6, 1), dimnames = list(NULL, NULL, "a"))
  for (h in 1:6) {
    expect_equal(tilt_quantity(p, "a", h, "annual_avg", hist), tilt_quantity(p, "a", h, "yoy", hist))
  }
  expect_equal(tilt_quantity(p, "a", 4, "annual_avg", hist), rep(100 * (exp(0.024) - 1), 2))
})

test_that("a median target on annual-average growth is met by the tilted weights", {
  p <- tl_paths()
  s <- tilt_quantity(p, "a", 6, "annual_avg", tl_hist)
  target <- unname(quantile(s, 0.35))
  out <- tilt_paths(p, tl_target("a", 6, "annual_avg", "median", target), tl_hist)
  expect_lt(abs(sum(out$weights[s <= target]) - 0.5), 1e-10)
})

test_that("tilted weights meet the targets, are positive and sum to one", {
  p <- tl_paths()
  tg <- rbind(tl_target("a", 6, "yoy", "median", 2.5), tl_target("r", 4, "level", "mean", 3.4))
  out <- tilt_paths(p, tg, tl_hist)
  expect_true(all(out$weights > 0))
  expect_equal(sum(out$weights), 1)
  expect_lt(max(abs(out$check$achieved - out$check$goal)), 1e-10)
  # The median target: the weighted median of the quantity is the target.
  s <- tilt_quantity(p, "a", 6, "yoy", tl_hist)
  expect_lt(abs(sum(out$weights[s <= 2.5]) - 0.5), 1e-10)
  expect_gt(out$ess, 1)
  expect_lt(out$ess, nrow(p))
})

test_that("no target, or a target already met, gives equal weights", {
  p <- tl_paths(500)
  none <- tilt_weights(matrix(numeric(0), 500, 0), numeric(0))
  expect_equal(none$weights, rep(1 / 500, 500))
  # A mean target equal to the sample mean is met by equal weights.
  s <- p[, 3, "r"]
  met <- tilt_weights(matrix(s), mean(s))
  expect_equal(met$iterations, 0L)
  expect_equal(met$weights, rep(1 / 500, 500))
  expect_equal(tilt_ess(met$weights), 500)
})

test_that("tilting the mean of a Gaussian sample shifts the mean and keeps the variance", {
  # For x ~ N(mu, s2), weights proportional to exp(gamma x) give N(mu + gamma s2, s2),
  # so a mean target mu' needs gamma = (mu' - mu) / s2 and leaves the variance unchanged.
  set.seed(2)
  x <- rnorm(2e5, 1, 2)
  out <- tilt_weights(matrix(x), 1.5)
  expect_equal(sum(out$weights * x), 1.5, tolerance = 1e-10)
  expect_equal(out$gamma, 0.5 / 4, tolerance = 0.02)
  w_var <- sum(out$weights * (x - 1.5)^2)
  expect_equal(w_var, var(x), tolerance = 0.01)
})

test_that("the tilted weights have the smallest divergence among feasible reweightings", {
  p <- tl_paths(300)
  tg <- tl_target("r", 2, "level", "mean", 3.3)
  mom <- tilt_moments(p, tg, tl_hist)
  sol <- tilt_weights(mom$g, mom$gbar)
  kl <- function(w) sum(w * log(length(w) * w))
  # Perturb along directions that keep both constraints (sum w = 1, sum w g = gbar).
  A <- cbind(1, mom$g[, 1])
  set.seed(3)
  for (k in 1:5) {
    d <- rnorm(300)
    d <- d - A %*% solve(crossprod(A), crossprod(A, d))      # project out the constraints
    w2 <- sol$weights + 0.2 * min(sol$weights) / max(abs(d)) * as.vector(d)
    expect_equal(sum(w2 * mom$g[, 1]), mom$gbar, tolerance = 1e-10)
    expect_gt(kl(w2), kl(sol$weights))
  }
  expect_equal(sol$kl, kl(sol$weights))
})

test_that("targets that cannot be met are errors", {
  p <- tl_paths(500)
  # A median above every draw, and a mean outside the range of the draws.
  expect_error(tilt_paths(p, tl_target("r", 1, "level", "median", 100), tl_hist), "one side")
  expect_error(tilt_paths(p, tl_target("r", 1, "level", "mean", -50), tl_hist), "one side")
  # Two different means for the same quantity: each feasible, jointly not.
  tg <- rbind(tl_target("r", 1, "level", "mean", 3.2), tl_target("r", 1, "level", "mean", 3.5))
  expect_error(tilt_paths(p, tg, tl_hist), "redundant or conflict|did not converge")
})

test_that("weighted quantiles are type 7 with equal weights and follow Akinshin (2023)", {
  set.seed(4)
  for (n in c(1, 2, 7, 100)) {
    x <- rnorm(n)
    probs <- c(0, 0.05, 0.16, 0.5, 0.84, 0.95, 1)
    expect_equal(weighted_quantiles(x, rep(1, n), probs), unname(quantile(x, probs)), tolerance = 1e-10)
  }
  # Akinshin (2023), Example 13: weighted median 11/3.
  expect_equal(weighted_quantiles(1:5, c(0.3, 0.1, 0, 0.1, 0.4), 0.5), 11 / 3)
  # Zero weights drop out; one dominant weight gives that value.
  x <- c(5, 1, 9, 3, 7)
  expect_equal(weighted_quantiles(x, c(1, 0, 1, 1, 0), c(0.25, 0.5)), unname(quantile(c(5, 9, 3), c(0.25, 0.5))))
  expect_equal(weighted_quantiles(x, c(0, 0, 1, 0, 0), c(0.1, 0.9)), c(9, 9))
  expect_equal(weighted_quantiles(x, NULL, 0.5), 5)
  expect_error(weighted_quantiles(x, c(-1, 1, 1, 1, 1), 0.5), "non-negative")
})
