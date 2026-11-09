# test-estimate-priors.R: prior construction (R/estimate_priors.R) and the
# independent marginal likelihood (R/estimate_ml.R).
#
# Uses a small simulated data set, so the tests are fast and need no data
# build. The full checks on the real data are in scripts/check_priors.R.

# Toy data: 3 variables, 80 quarters. Two random walks around positive levels
# and one white-noise series (like realized volatility).
set.seed(123)
toy_y <- cbind(
  a = 100 + cumsum(rnorm(80, 0.5, 1)),
  b = 50 + cumsum(rnorm(80, 0.2, 0.5)),
  c = 20 + rnorm(80, 0, 3)
)
toy_series <- data.frame(short_name = c("a", "b", "c"), prior = c("RW", "RW", "WN"),
                         stringsAsFactors = FALSE)
toy_prior_cfg <- list(
  lambda = list(mode = 0.2, sd = 0.4, min = 0.0001, max = 5),
  alpha = 2,
  soc = list(mode = 1, sd = 1, min = 0.0001, max = 50),
  soc_ybar = "estimation_sample",
  constant_var = 1e7,
  psi = list(shape = 0.0004, scale = 0.0004, start = "ar_variance", bound_factor = 100)
)

test_that("prior_mean_vector gives 1 for RW and 0 for WN", {
  expect_equal(prior_mean_vector(toy_series), c(a = 1, b = 1, c = 0))
  bad <- toy_series
  bad$prior[1] <- "AR"
  expect_error(prior_mean_vector(bad), "RW or WN")
})

test_that("ar_residual_variance matches lm()", {
  x <- toy_y[, "a"]
  p <- 3
  n <- length(x)
  fit <- lm(x[(p + 1):n] ~ x[p:(n - 1)] + x[(p - 1):(n - 2)] + x[(p - 2):(n - 3)])
  ref <- sum(residuals(fit)^2) / fit$df.residual
  expect_equal(ar_residual_variance(x, p)$variance, ref)
})

test_that("var_data_matrices lays out the constant and lags as documented", {
  y  <- matrix(c(1, 2, 3, 4, 10, 20, 30, 40), ncol = 2)   # 4 quarters, 2 variables
  dm <- var_data_matrices(y, lags = 2)
  expect_equal(dm$Y, unname(y[3:4, ]))
  # Columns: constant, lag 1 of (v1, v2), lag 2 of (v1, v2).
  expect_equal(dm$X, rbind(c(1, 2, 20, 1, 10), c(1, 3, 30, 2, 20)))
})

test_that("prior_mean_matrix and minnesota_omega follow the formulas", {
  b <- prior_mean_matrix(c(1, 0), lags = 2)
  expect_equal(dim(b), c(5, 2))
  expect_equal(b[2, 1], 1)
  expect_equal(sum(b != 0), 1)
  omega <- minnesota_omega(lambda = 0.2, alpha = 2, psi = c(4, 1), lags = 2, constant_var = 1e7)
  expect_equal(omega, c(1e7, 0.04 / 4, 0.04 / 1, 0.04 / (4 * 4), 0.04 / (4 * 1)))
})

test_that("log_ml_niw (rescaled form) equals the formula computed directly", {
  # Well-conditioned case, so the direct inversion is reliable.
  dm    <- var_data_matrices(scale(toy_y), lags = 1)
  b     <- prior_mean_matrix(c(1, 1, 0), lags = 1)
  omega <- c(10, rep(0.5, 3))
  psi   <- c(1, 2, 0.5)
  d     <- 5
  A     <- crossprod(dm$X) + diag(1 / omega)
  B_hat <- solve(A, crossprod(dm$X, dm$Y) + diag(1 / omega) %*% b)
  S     <- crossprod(dm$Y - dm$X %*% B_hat) + t(B_hat - b) %*% diag(1 / omega) %*% (B_hat - b)
  n <- 3; N <- nrow(dm$Y)
  direct <- -n * N / 2 * log(pi) + log_mv_gamma(n, (N + d) / 2) - log_mv_gamma(n, d / 2) -
    n / 2 * sum(log(omega)) - n / 2 * log(det(A)) + d / 2 * sum(log(psi)) -
    (N + d) / 2 * log(det(diag(psi) + S))
  expect_equal(log_ml_niw(dm$Y, dm$X, b, omega, psi, d)$log_ml, direct, tolerance = 1e-10)
})

test_that("build_priors sets psi start and bounds and the estimated hyperparameters", {
  built <- build_priors(toy_prior_cfg, toy_series, toy_y, lags = 2)
  expect_equal(built$b, c(a = 1, b = 1, c = 0))
  expect_equal(built$psi_table$start, built$psi_table$ar_variance)
  expect_equal(built$psi_table$min, built$psi_table$start / 100)
  expect_equal(built$psi_table$max, built$psi_table$start * 100)
  expect_equal(built$priors$hyper, c("lambda", "soc", "psi"))
  expect_equal(built$priors$psi$coef, list(k = 0.0004, theta = 0.0004))
})

test_that("the presample soc dummy uses the first p rows", {
  cfg <- toy_prior_cfg
  cfg$soc_ybar <- "presample"
  built <- build_priors(cfg, toy_series, toy_y, lags = 2)
  expect_equal(built$soc_ybar$ybar0, colMeans(toy_y[1:2, ]))
  dum <- built$priors$soc$fun(Y = NULL, lags = 2, par = 2)
  expect_equal(dum$Y, diag(colMeans(toy_y[1:2, ])) / 2)
  expect_equal(dim(dum$X), c(3, 7))
})

test_that("BVAR's log posterior matches the direct implementation (BVAR 1.0.5)", {
  lags <- 2
  ybar <- colMeans(toy_y[lags + seq_len(lags), ])        # BVAR's centring of the soc dummies
  # Case 1: all variables random walk. The dummy observations are fitted
  # exactly by b, so BVAR and the textbook formula must agree.
  ser_rw <- toy_series
  ser_rw$prior <- "RW"
  built  <- build_priors(toy_prior_cfg, ser_rw, toy_y, lags)
  fit    <- fit_for_mode(toy_y, lags, built$priors, seed = 1)
  h      <- c(lambda = 0.3, soc = 0.5, psi1 = 1, psi2 = 0.3, psi3 = 9)
  bv     <- bvar_log_posterior(fit, h)
  own    <- log_posterior_direct(toy_y, lags, 0.3, 0.5, c(1, 0.3, 9), built$b, toy_prior_cfg,
                                 ybar, "posterior")$log_posterior
  expect_equal(bv, own, tolerance = 1e-8)

  # Case 2: one white-noise variable. BVAR evaluates the dummy-only marginal
  # likelihood at the prior mean b; this test records that behaviour, so a
  # change in the package is noticed.
  built  <- build_priors(toy_prior_cfg, toy_series, toy_y, lags)
  fit    <- fit_for_mode(toy_y, lags, built$priors, seed = 1)
  bv     <- bvar_log_posterior(fit, h)
  own_b  <- log_posterior_direct(toy_y, lags, 0.3, 0.5, c(1, 0.3, 9), built$b, toy_prior_cfg,
                                 ybar, "prior_mean")$log_posterior
  own_tb <- log_posterior_direct(toy_y, lags, 0.3, 0.5, c(1, 0.3, 9), built$b, toy_prior_cfg,
                                 ybar, "posterior")$log_posterior
  expect_equal(bv, own_b, tolerance = 1e-8)
  expect_gt(abs(bv - own_tb), 1e-3)
})
