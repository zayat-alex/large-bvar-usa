# test-hyper-stability.R: the pieces of the hyperparameter stability check
# (R/hyper_stability_checks.R).
#
# Small simulated data, as in test-estimate-bvar.R. The check on the real
# data is scripts/check_hyper_stability.R.

set.seed(321)
hs_y <- cbind(
  a = 100 + cumsum(rnorm(80, 0.5, 1)),
  b = 50 + cumsum(rnorm(80, 0.2, 0.5)),
  c = 20 + rnorm(80, 0, 3)
)
hs_series <- data.frame(short_name = c("a", "b", "c"), prior = c("RW", "RW", "WN"),
                        stringsAsFactors = FALSE)
hs_prior_cfg <- list(
  lambda = list(mode = 0.2, sd = 0.4, min = 0.0001, max = 5),
  alpha = 2,
  soc = list(mode = 1, sd = 1, min = 0.0001, max = 50),
  soc_ybar = "presample",
  constant_var = 1e7,
  psi = list(shape = 0.0004, scale = 0.0004, start = "ar_variance", bound_factor = 100)
)
hs_setup <- estimation_setup(hs_y, hs_series, hs_prior_cfg, lags = 2)

test_that("sample_to_origin keeps the quarters up to the origin", {
  d <- data.frame(quarter = parse_quarter(c("2015Q4", "2016Q1", "2016Q2", "2016Q3")), x = 1:4)
  out <- sample_to_origin(d, "2016Q2")
  expect_equal(format_quarter(out$quarter), c("2015Q4", "2016Q1", "2016Q2"))
  expect_error(sample_to_origin(d, "2017Q1"), "not a quarter of the data")
})

test_that("the parallel mode search gives the sequential result", {
  seq_mode <- find_hyper_mode(hs_setup, start_factors = c(1, 0.5, 2))
  par_mode <- find_hyper_mode_parallel(hs_setup, c(1, 0.5, 2), 2000,
                                       r_dir = file.path(project_root(), "R"))
  expect_equal(par_mode$theta_mode, seq_mode$theta_mode)
  expect_equal(par_mode$table, seq_mode$table)
  expect_equal(par_mode$best, seq_mode$best)
})

test_that("hyper_shift_table measures shifts in posterior s.d.", {
  theta <- log(hs_setup$start)
  mode_now  <- list(theta_mode = theta)
  mode_past <- list(theta_mode = theta + c(0.2, -0.1, 0, 0, 0))
  prop <- list(approx_sd = setNames(rep(0.1, 5), names(theta)))
  # Current draws with s.d. 0.1 on the log scale: two points at +-0.1 around the mode.
  hyper_now <- exp(rbind(theta + 0.1, theta - 0.1))
  colnames(hyper_now) <- names(theta)
  s <- hyper_shift_table(mode_past, prop, mode_now, prop, hyper_now, hs_setup, c("a", "b", "c"))
  sd_now <- sd(c(0.1, -0.1))
  expect_equal(s$log_ratio, c(0.2, -0.1, 0, 0, 0))
  expect_equal(s$shift_in_sd_now, c(0.2, -0.1, 0, 0, 0) / sd_now)
  expect_equal(s$shift_in_sd_past, c(2, -1, 0, 0, 0))
  expect_equal(s$hyperparameter, c("lambda", "soc", "psi_a", "psi_b", "psi_c"))
})

test_that("draws_at_mode is reproducible and matches niw_posterior's mean", {
  theta <- find_hyper_mode(hs_setup, start_factors = 1)$theta_mode
  d1 <- draws_at_mode(theta, hs_setup, 2000, seed = 7)
  d2 <- draws_at_mode(theta, hs_setup, 2000, seed = 7)
  expect_identical(d1, d2)
  post <- niw_posterior(complete_hyper(exp(theta), hs_setup), hs_setup)
  B_mean <- apply(d1$beta, c(2, 3), mean)
  # Monte Carlo error of a mean of 2,000 draws is sd / sqrt(2000); allow 5 of them.
  B_sd <- apply(d1$beta, c(2, 3), sd)
  expect_lt(max(abs(B_mean - post$B_hat) / B_sd), 5 / sqrt(2000))
})

test_that("compare_bands reads shifts and widths as shares of the band", {
  set.seed(1)
  ref <- array(rnorm(5000 * 2 * 1), c(5000, 2, 1), dimnames = list(NULL, NULL, "a"))
  alt <- ref + 0.5            # same spread, median shifted by 0.5
  alt[, 2, ] <- ref[, 2, ] * 2  # horizon 2: every quantile doubled, so twice the width
  out <- compare_bands(ref, alt, ref, "a", 1:2)
  width <- diff(quantile(ref[, 1, 1], c(0.05, 0.95), names = FALSE))
  expect_equal(out$median_diff_share[1], 0.5 / width, tolerance = 1e-12)
  expect_equal(out$width_ratio[1], 1, tolerance = 1e-12)
  expect_equal(out$width_ratio[2], 2, tolerance = 1e-12)
  expect_equal(out$median_diff_share_seed, c(0, 0))
  # A base value shifts every quantile but not the comparison.
  out_b <- compare_bands(ref, alt, ref, "a", 1:2, base = c(a = 3))
  expect_equal(out_b$median_diff_share, out$median_diff_share)
  expect_equal(out_b$q50_ref, out$q50_ref - 3)
})

test_that("growth_measure gives YoY and QoQ annualized growth from log levels", {
  # History and paths on a straight line rising by 0.5 per quarter (100 x log):
  # YoY growth 2 percent, QoQ annualized 2 percent; the raw series keeps its level.
  history <- cbind(g = 100 + 0.5 * (1:8), r = rep(4, 8))
  paths <- array(NA_real_, c(2, 3, 2), dimnames = list(NULL, NULL, c("g", "r")))
  for (i in 1:2) paths[i, , "g"] <- 100 + 0.5 * (9:11)
  paths[, , "r"] <- 5
  is_log <- c(TRUE, FALSE)
  a <- growth_measure(paths, history, is_log, "annual")
  q <- growth_measure(paths, history, is_log, "quarterly")
  expect_equal(a[, , "g"], matrix(2, 2, 3))
  expect_equal(q[, , "g"], matrix(2, 2, 3))
  expect_equal(a[, , "r"], matrix(5, 2, 3))
  expect_true(all(is.na(q[, , "r"])))
})

test_that("irf_annual_measure turns a level response into the effect on YoY growth", {
  # A permanent 1-point level shift from impact raises YoY growth by 1 for
  # four quarters, then not at all; a raw series keeps its response.
  irf <- array(1, c(1, 7, 2), dimnames = list(NULL, NULL, c("g", "r")))
  out <- irf_annual_measure(irf, c(TRUE, FALSE))
  expect_equal(out[1, , "g"], c(1, 1, 1, 1, 0, 0, 0))
  expect_equal(out[1, , "r"], rep(1, 7))
})

test_that("historical_scale and summarize_medians read changes against history", {
  set.seed(2)
  y <- cbind(g = 100 + cumsum(rnorm(60, 0.5, 1)), r = 5 + cumsum(rnorm(60, 0, 0.2)))
  rownames(y) <- format_quarter(seq(as.Date("2000-01-01"), by = "quarter", length.out = 60))
  s <- historical_scale(y, c(TRUE, FALSE), end = rownames(y)[40])
  expect_equal(s$sd[1], sd(diff(y[1:40, "g"], lag = 4)))
  expect_equal(s$qoq_sd[1], sd(4 * diff(y[1:40, "g"])))
  expect_true(is.na(s$qoq_sd[2]))
  cmp <- data.frame(variable = c("g", "g"), horizon = 1:2, median_ref = 0, median_alt = 0,
                    diff = c(0.5, -1), diff_seed = c(0.1, 0.2))
  out <- summarize_medians(cmp, s, list(`1-2` = 1:2))
  expect_equal(out$max_abs_diff, 1)
  expect_equal(out$diff_in_sd, 1 / s$sd[1])
  expect_equal(out$diff_in_range, 1 / s$range_10_90[1])
  expect_equal(out$noise_in_sd, 0.2 / s$sd[1])
})
