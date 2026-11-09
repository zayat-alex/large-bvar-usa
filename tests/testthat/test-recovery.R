# test-recovery.R: helpers of the simulation recovery check
# (R/recovery_checks.R): the simulator, band coverage and coefficient groups.

rc_B1 <- matrix(c(0.5, 0.1, 0.0, 0.3), 2)
rc_B <- rbind(c(1, -0.5), t(rc_B1))                  # VAR(1), K x n = 3 x 2
rc_Sigma <- matrix(c(1, 0.3, 0.3, 0.5), 2)
rc_init <- matrix(c(2, 1), 1, 2, dimnames = list(NULL, c("a", "b")))

test_that("simulate_var follows the VAR recursion with Cholesky shocks", {
  sim <- simulate_var(rc_B, rc_Sigma, rc_init, 1, 5, seed = 3)
  expect_equal(dim(sim), c(5L, 2L))
  # Recover the shocks from the recursion: e_t = y_t - c - B_1 y_{t-1}.
  y <- rbind(rc_init, sim)
  e <- t(sapply(2:6, function(t) y[t, ] - rc_B[1, ] - drop(y[t - 1, ] %*% rc_B[-1, ])))
  # With the same seed, the shocks are L u_t with u_t from rnorm(n) per quarter.
  set.seed(3)
  u <- matrix(rnorm(10), 5, 2, byrow = TRUE)
  expect_equal(unname(e), u %*% t(t(chol(rc_Sigma))), tolerance = 1e-12)
})

test_that("simulated moments match the stationary distribution", {
  # Long simulation of a stable VAR(1): sample covariance of the shocks
  # matches Sigma and the sample mean matches (I - B_1)^-1 c (B stores B_1').
  sim <- simulate_var(rc_B, rc_Sigma, rc_init, 1, 40000, seed = 11)
  y <- rbind(rc_init, sim)
  e <- y[-1, ] - matrix(rc_B[1, ], nrow(y) - 1, 2, byrow = TRUE) - y[-nrow(y), ] %*% rc_B[-1, ]
  expect_equal(unname(cov(e)), rc_Sigma, tolerance = 0.03)
  mu <- solve(diag(2) - rc_B1, rc_B[1, ])
  expect_equal(unname(colMeans(sim)), mu, tolerance = 0.03)
})

test_that("band coverage flags the truth inside the 68% and 90% bands", {
  draws <- array(rep(1:100, 3), c(100, 3))            # 100 draws of three cells, values 1..100
  cov <- band_coverage(draws, c(50, 10, 99))
  expect_equal(as.vector(cov$in68), c(TRUE, FALSE, FALSE))
  expect_equal(as.vector(cov$in90), c(TRUE, TRUE, FALSE))
  # Shape is kept for matrices of cells.
  d <- array(rnorm(200 * 6), c(200, 2, 3))
  cov <- band_coverage(d, matrix(0, 2, 3))
  expect_equal(dim(cov$in90), c(2L, 3L))
  expect_true(all(cov$in90))
})

test_that("coefficient groups follow the layout of B", {
  g <- coefficient_groups(n = 2, lags = 2)            # K = 5
  expect_equal(g[1, ], c("constant", "constant"))
  expect_equal(diag(g[2:3, ]), c("own first lag", "own first lag"))
  expect_equal(diag(g[4:5, ]), c("own longer lags", "own longer lags"))
  expect_equal(c(g[3, 1], g[2, 2], g[5, 1], g[4, 2]), rep("cross lags", 4))
  expect_equal(sum(g == "cross lags"), 2 * 2 * (2 - 1))
})
