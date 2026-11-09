# test-forecast.R: the unconditional forecast simulator (R/forecast.R) and
# the exact forecast moments used to check it (R/forecast_checks.R).
#
# A small VAR(2) with two variables and known coefficients, so the right
# answers are known exactly.

fc_B1 <- matrix(c(0.6, 0.1, 0.2, 0.5), 2)        # (B_1)_ij: effect of variable j on i
fc_B2 <- matrix(c(0.2, 0.0, -0.1, 0.3), 2)
fc_c  <- c(0.5, -0.2)
fc_B  <- rbind(fc_c, t(fc_B1), t(fc_B2))          # K x n: constant, lag 1 rows, lag 2 rows
fc_Sigma <- matrix(c(1, 0.4, 0.4, 0.5), 2)
fc_hist <- matrix(c(1, 2, 3, 0.5, 1.5, 1.0), 3, 2, dimnames = list(NULL, c("a", "b")))
one_draw <- function(B, S) list(beta = array(B, c(1, dim(B))), sigma = array(S, c(1, dim(S))))

test_that("the no-shock path follows the VAR recursion", {
  d <- one_draw(fc_B, fc_Sigma)
  fc <- forecast_unconditional(d$beta, d$sigma, fc_hist, lags = 2, horizon = 3, seed = 1, shocks = FALSE)
  y1 <- fc_c + fc_B1 %*% fc_hist[3, ] + fc_B2 %*% fc_hist[2, ]
  y2 <- fc_c + fc_B1 %*% y1 + fc_B2 %*% fc_hist[3, ]
  y3 <- fc_c + fc_B1 %*% y2 + fc_B2 %*% y1
  expect_equal(unname(fc[1, , ]), unname(rbind(t(y1), t(y2), t(y3))))
  expect_equal(dimnames(fc)[[3]], c("a", "b"))
})

test_that("only the last p rows of the history matter", {
  d <- one_draw(fc_B, fc_Sigma)
  h2 <- fc_hist; h2[1, ] <- c(100, -100)
  expect_equal(forecast_unconditional(d$beta, d$sigma, fc_hist, 2, 4, seed = 3),
               forecast_unconditional(d$beta, d$sigma, h2, 2, 4, seed = 3))
})

test_that("the exact forecast covariance has the textbook form", {
  # h = 1: Sigma; h = 2: Sigma + B_1 Sigma B_1'.
  ex <- forecast_moments_exact(fc_B, fc_Sigma, fc_hist, lags = 2, horizon = 2)
  expect_equal(ex$cov[[1]], fc_Sigma)
  expect_equal(ex$cov[[2]], fc_Sigma + fc_B1 %*% fc_Sigma %*% t(fc_B1))
})

test_that("simulated forecasts match the exact mean and covariance", {
  res <- check_known_var(fc_B, fc_Sigma, fc_hist, lags = 2, horizon = 8, n_sim = 20000, seed = 11)
  expect_lt(max(abs(res$table$z_mean)), 4)
  expect_true(all(abs(res$table$var_ratio - 1) < 4 * res$table$var_ratio_mc_se))
  # Correlation between the two variables at the last horizon.
  r_exact <- cov2cor(res$cov_last_exact)[1, 2]
  r_sim   <- cov2cor(res$cov_last_sim)[1, 2]
  expect_lt(abs(r_sim - r_exact), 0.03)
})

test_that("BVAR::predict() given the same draws has the same no-noise path", {
  # BVAR's mean path (compute_fcast) and the project's no-shock path must
  # coincide; its noise is added afterwards (see R/forecast_checks.R).
  d <- one_draw(fc_B, fc_Sigma)
  own <- forecast_unconditional(d$beta, d$sigma, fc_hist, 2, 5, seed = 1, shocks = FALSE)[1, , ]
  pkg <- BVAR:::compute_fcast(Y = unname(fc_hist), K = 5, M = 2, N = 3, lags = 2, horizon = 5,
                              beta_comp = BVAR:::get_beta_comp(fc_B, K = 5, M = 2, lags = 2),
                              beta_const = fc_c)
  expect_equal(unname(own), unname(pkg))
})
