# test-irf.R: impulse responses (R/irf.R) and their checks (R/irf_checks.R).
#
# A VAR(2) with four variables and known coefficients: two "slow" variables
# (a, b), a policy variable (r) and a "fast" variable (f), so the recursive
# identification has blocks with more than one ordering.

ir_vars <- c("a", "b", "r", "f")
ir_B1 <- matrix(c(0.6, 0.1, 0.0, 0.1,
                  0.2, 0.5, 0.1, 0.0,
                  0.0, 0.1, 0.7, 0.2,
                  0.1, 0.0, 0.1, 0.4), 4)                     # (B_1)_ij, filled by column
ir_B2 <- matrix(c(0.2, 0.0, 0.05, 0.0,
                  -0.1, 0.2, 0.0, 0.05,
                  0.0, 0.0, 0.1, 0.0,
                  0.0, 0.1, 0.0, 0.2), 4)
ir_B <- rbind(c(0.5, -0.2, 0.1, 0.3), t(ir_B1), t(ir_B2))     # K x n
ir_Sigma <- matrix(c(1.0, 0.4, 0.3, 0.2,
                     0.4, 0.6, 0.2, 0.1,
                     0.3, 0.2, 0.8, 0.3,
                     0.2, 0.1, 0.3, 0.9), 4)
ir_hist <- matrix(c(1.0, 0.2, -0.5, 0.8, 1.1,
                    0.3, 0.0, -0.4, 0.6, 1.2,
                    0.9, 0.1, 0.4, 0.7, -0.2,
                    0.5, 0.6, 0.2, -0.1, 0.3), 5, 4, dimnames = list(NULL, ir_vars))
ir_H <- 8
ir_rec <- list(slow = c("a", "b"), policy = "r", fast = "f")

test_that("responses follow the VAR recursion from a unit impact", {
  a <- cholesky_impact(ir_Sigma, ir_vars, ir_vars, "b")
  expect_equal(unname(a[c("a", "b")]), c(0, 1))              # a ordered before b: no impact
  r <- cholesky_irf(ir_B, ir_Sigma, 2, ir_H, ir_vars, ir_vars, "b")
  # By hand: r_0 = a, r_1 = B_1 r_0, r_h = B_1 r_{h-1} + B_2 r_{h-2}.
  hand <- matrix(0, ir_H + 1, 4)
  hand[1, ] <- a
  hand[2, ] <- ir_B1 %*% hand[1, ]
  for (h in 2:ir_H) hand[h + 1, ] <- ir_B1 %*% hand[h, ] + ir_B2 %*% hand[h - 1, ]
  expect_equal(unname(r), hand, tolerance = 1e-12)
})

test_that("the GIRF equals the Cholesky IRF with the shocked variable first", {
  g <- girf(ir_B, ir_Sigma, 2, ir_H, ir_vars, ir_vars)
  for (j in ir_vars) {
    ch <- cholesky_irf(ir_B, ir_Sigma, 2, ir_H, ir_vars, c(j, setdiff(ir_vars, j)), j)
    expect_equal(g[, , j], ch, tolerance = 1e-12)
    expect_equal(unname(g[1, j, j]), 1)
  }
  # Impact: the regression of each forecast error on that of the shocked variable.
  expect_equal(unname(g[1, , "r"]), ir_Sigma[, 3] / ir_Sigma[3, 3])
})

test_that("the GIRF equals the conditional-forecast effect (Crump et al. 2021, Section 2.3)", {
  g <- girf(ir_B, ir_Sigma, 2, ir_H, ir_vars, ir_vars)
  for (j in ir_vars) {
    e <- irf_by_conditioning(ir_B, ir_Sigma, ir_hist, 2, ir_H, setNames(1, j))
    expect_equal(e, g[, , j], tolerance = 1e-10)
  }
  # By linearity, a shock of 0.25 gives a quarter of the response.
  e <- irf_by_conditioning(ir_B, ir_Sigma, ir_hist, 2, ir_H, c(r = 0.25))
  expect_equal(e, 0.25 * g[, , "r"], tolerance = 1e-10)
})

test_that("the recursive IRF equals conditioning slow variables at baseline", {
  r <- recursive_irf(ir_B, ir_Sigma, 2, ir_H, ir_vars, ir_rec$slow, ir_rec$policy, ir_rec$fast)
  expect_equal(unname(r[1, c("a", "b")]), c(0, 0))
  expect_equal(unname(r[1, "r"]), 1)
  e <- irf_by_conditioning(ir_B, ir_Sigma, ir_hist, 2, ir_H, c(a = 0, b = 0, r = 1))
  expect_equal(e, r, tolerance = 1e-10)
  # It differs from the GIRF of r, which lets the slow variables move.
  expect_gt(max(abs(r - girf(ir_B, ir_Sigma, 2, ir_H, ir_vars, "r")[, , 1])), 0.05)
  # The order within the slow block does not matter.
  r2 <- recursive_irf(ir_B, ir_Sigma, 2, ir_H, ir_vars, rev(ir_rec$slow), ir_rec$policy, ir_rec$fast)
  expect_equal(r2, r, tolerance = 1e-12)
  res <- irf_equivalence(ir_B, ir_Sigma, ir_hist, 2, ir_H, ir_vars, ir_rec)
  expect_true(all(res$rel_diff < 1e-10))
})

test_that("blocks of a recursive identification are checked", {
  expect_error(check_recursive_blocks(c("a", "b"), "r", c("f", "a"), ir_vars), "more than one block")
  expect_error(check_recursive_blocks(c("a"), "r", "f", ir_vars), "in no block: b")
  expect_error(check_recursive_blocks(c("a", "b", "x"), "r", "f", ir_vars), "unknown")
  expect_error(check_recursive_blocks(c("a", "b"), c("r", "f"), character(0), ir_vars), "exactly one")
  expect_error(girf(ir_B, ir_Sigma, 2, ir_H, ir_vars, "x"), "unknown")
})

test_that("irf_draws and irf_quantiles have the documented shapes", {
  beta <- array(rep(ir_B, each = 3), c(3, dim(ir_B)))
  sigma <- array(rep(ir_Sigma, each = 3), c(3, dim(ir_Sigma)))
  sigma[2, , ] <- sigma[2, , ] * 2                         # unit shocks: the scale of Sigma cancels
  d <- irf_draws(beta, sigma, 2, ir_H, ir_vars, c("a", "r"), ir_rec)
  expect_named(d, c("girf_a", "girf_r", "recursive_r"))
  expect_equal(dim(d$girf_a), c(3, ir_H + 1, 4))
  expect_equal(d$girf_r[2, , ], d$girf_r[1, , ], tolerance = 1e-12)
  q <- irf_quantiles(d$recursive_r, c(0.05, 0.5, 0.95))
  expect_equal(nrow(q), (ir_H + 1) * 4)
  expect_named(q, c("horizon", "variable", "q05", "q50", "q95"))
})

test_that("the project's IRFs equal BVAR::irf() after reordering and rescaling", {
  res <- irf_bvar_agreement(ir_B, ir_Sigma, ir_hist, 2, ir_H, ir_vars, ir_rec)
  expect_true(all(res$max_abs_diff < 1e-12))
})

test_that("IRF settings are read from the configuration", {
  cfg <- yaml::read_yaml(file.path(project_root(), "config", "run.yml"))
  series <- read_series_config(file.path(project_root(), "config", "series.csv"))
  s <- irf_settings(cfg, series$short_name)
  expect_equal(s$recursive$policy, "gs2")
  # Slow block = the variables listed before the policy variable in the series order.
  expect_equal(s$recursive$slow, series$short_name[seq_len(match("gs2", series$short_name) - 1)])
  bad <- cfg
  bad$irf$recursive$fast <- head(unlist(bad$irf$recursive$fast), -1)
  expect_error(irf_settings(bad, series$short_name), "in no block")
})
