# test-estimation-diagnostics.R: residual diagnostics (R/estimation_diagnostics.R).
# Known VAR(1) with two variables; the check on the production runs is
# scripts/check_estimation_diagnostics.R.

set.seed(11)
ed_B <- rbind(c(0.5, 1), c(0.9, 0.1), c(0, 0.6))   # constant; lag-1 rows for a, b
ed_Sigma <- matrix(c(1, 0.5, 0.5, 2), 2)
ed_T <- 300
ed_eps <- matrix(rnorm(ed_T * 2), ed_T) %*% chol(ed_Sigma)
ed_y <- matrix(0, ed_T, 2, dimnames = list(NULL, c("a", "b")))
for (t in 2:ed_T) ed_y[t, ] <- c(1, ed_y[t - 1, ]) %*% ed_B + ed_eps[t, ]
rownames(ed_y) <- format_quarter(seq(as.Date("1950-01-01"), by = "quarter", length.out = ed_T))

test_that("var_residuals recovers the shocks of a known VAR, and averages the draws", {
  # Three draws around the true B whose mean is the true B but whose
  # element-wise median is not: the residuals must be the simulated shocks
  # exactly (posterior mean, not median).
  d <- c(0.1, -0.02, -0.08)
  draws <- list(lags = 1,
                beta = aperm(array(c(ed_B + d[1], ed_B + d[2], ed_B + d[3]), c(3, 2, 3)), c(3, 1, 2)),
                sigma = aperm(array(rep(ed_Sigma, 3), c(2, 2, 3)), c(3, 1, 2)))
  r <- var_residuals(ed_y, draws)
  expect_equal(unname(r$e), unname(ed_eps[-1, ]), tolerance = 1e-12)
  expect_equal(rownames(r$e), rownames(ed_y)[-1])
  expect_equal(unname(r$z[, 2]), unname(ed_eps[-1, 2]) / sqrt(2), tolerance = 1e-12)
})

test_that("residual_outliers lists the cells beyond the threshold, largest first", {
  z <- matrix(c(0, 6, -7, 1), 2, dimnames = list(c("2000Q1", "2000Q2"), c("a", "b")))
  out <- residual_outliers(z, 5)
  expect_equal(out$z, c(-7, 6))
  expect_equal(out$quarter, c("2000Q1", "2000Q2"))
  expect_equal(nrow(residual_outliers(z, 10)), 0)
})

test_that("residual_autocorrelation separates white noise from an AR(1)", {
  set.seed(5)
  wn <- rnorm(400)
  ar <- as.numeric(stats::filter(rnorm(400), 0.6, method = "recursive"))
  a <- residual_autocorrelation(cbind(wn = wn, ar = ar), 8)
  expect_lt(abs(a$acf_1[a$variable == "wn"]), 0.15)
  expect_gt(a$acf_1[a$variable == "ar"], 0.45)
  expect_lt(a$ljung_box_p[a$variable == "ar"], 1e-6)
})

test_that("rows_before keeps the quarters before a cut-off", {
  m <- matrix(1:4, 4, dimnames = list(c("2019Q3", "2019Q4", "2020Q1", "2020Q2"), "x"))
  expect_equal(rownames(rows_before(m, "2020Q1")), c("2019Q3", "2019Q4"))
})
