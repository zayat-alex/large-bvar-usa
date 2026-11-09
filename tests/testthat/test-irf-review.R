# test-irf-review.R: the spread widening and the expectation verdicts of the
# impulse-response review (R/irf_review.R). Known VAR(2) with three
# variables; the review of the production runs is scripts/check_irf_review.R.

ir_vars <- c("a", "b", "c")
ir_B <- rbind(c(0.1, 0.2, 0), diag(c(0.5, 0.4, 0.3)), matrix(c(0.1, 0, 0.05, 0, 0.1, 0, 0.02, 0, 0.1), 3))
ir_Sigma <- matrix(c(1, 0.6, 0.2, 0.6, 1.5, 0.3, 0.2, 0.3, 0.8), 3)
set.seed(3)
ir_hist <- matrix(rnorm(30), 10, 3, dimnames = list(NULL, ir_vars))

test_that("spread_impact moves one variable by the size and holds the other", {
  a <- spread_impact(ir_Sigma, ir_vars, "b", "a", size = 2)
  expect_equal(a[["b"]], 2)
  expect_equal(a[["a"]], 0, tolerance = 1e-14)
  # Conditional mean of a normal vector given e_b = 2, e_a = 0.
  S <- c(2, 1)
  expect_equal(unname(a), as.vector(ir_Sigma[, S] %*% solve(ir_Sigma[S, S], c(2, 0))))
})

test_that("the spread response equals the conditional minus the unconditional forecast mean", {
  expect_lt(spread_vs_conditional(ir_B, ir_Sigma, ir_hist, 2, 8, "b", "a", 1), 1e-12)
})

test_that("spread_irf_draws propagates each draw's own impact", {
  beta <- aperm(array(c(ir_B, ir_B * 0.9), c(7, 3, 2)), c(3, 1, 2))
  sigma <- aperm(array(c(ir_Sigma, ir_Sigma * 2), c(3, 3, 2)), c(3, 1, 2))
  d <- spread_irf_draws(beta, sigma, 2, 6, ir_vars, "b", "a")
  expect_equal(unname(d[2, , ]), unname(propagate_impact(ir_B * 0.9, spread_impact(ir_Sigma * 2, ir_vars, "b", "a"),
                                                         2, 6, ir_vars)))
})

test_that("expectation_results gives the verdicts and scales the shock", {
  q <- data.frame(horizon = rep(0:2, 3), variable = rep(c("up", "flat", "down"), each = 3),
                  q05 = 0, q16 = c(0.1, 0.2, 0.3, -0.2, -0.1, -0.3, -0.5, -0.6, -0.4),
                  q25 = 0, q50 = c(0.2, 0.5, 0.4, 0.1, 0.2, -0.1, -0.3, -0.4, -0.2),
                  q75 = 0, q84 = c(0.3, 0.8, 0.5, 0.4, 0.5, 0.1, -0.1, -0.2, -0.1), q95 = 0)
  rows <- data.frame(variable = c("up", "flat", "down", "down"), sign = c(1, 1, 1, -1),
                     h_from = 0, h_to = 2, scale = c(1, 1, 1, -10))
  out <- expectation_results(q, rows)
  expect_equal(out$verdict, c("as expected", "weak", "opposite", "opposite"))
  expect_equal(out$horizon[1], 1)
  expect_equal(out$median[1], 0.5)
  # Scaled by -10 (a negative shock), "down" rises by 4 at horizon 1 with band
  # [2, 6]: "opposite" when a fall is expected (row 4 above), "as expected"
  # when a rise is expected (below).
  rows2 <- data.frame(variable = "down", sign = 1, h_from = 0, h_to = 2, scale = -10)
  out2 <- expectation_results(q, rows2)
  expect_equal(c(out2$median, out2$q16, out2$q84), c(4, 2, 6))
  expect_equal(out2$verdict, "as expected")
})
