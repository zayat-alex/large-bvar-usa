# test-forecast-conditional.R: the two conditional forecast methods
# (R/forecast_conditional.R: both methods, the dispatcher and the input
# check), and the textbook reference
# (R/conditional_checks.R).
#
# A VAR(2) with three variables and known coefficients, so the conditional
# moments can be computed exactly by the textbook formula and compared.

cf_B1 <- matrix(c(0.6, 0.1, 0.0, 0.2, 0.5, 0.1, 0.0, 0.1, 0.7), 3)   # (B_1)_ij
cf_B2 <- matrix(c(0.2, 0.0, 0.05, -0.1, 0.2, 0.0, 0.0, 0.0, 0.1), 3)
cf_B  <- rbind(c(0.5, -0.2, 0.1), t(cf_B1), t(cf_B2))                 # K x n
cf_Sigma <- matrix(c(1, 0.4, 0.2, 0.4, 0.5, 0.1, 0.2, 0.1, 0.8), 3)
cf_hist <- matrix(c(1.0, 0.2, -0.5, 0.8, 1.1, 0.3, 0.0, -0.4, 0.6, 1.2, 0.9, 0.1, 0.4, 0.7, -0.2),
                  5, 3, dimnames = list(NULL, c("a", "b", "c")))
cf_H <- 6
cf_methods <- c("gaussian", "kalman")

#' Conditions with a single cell, several cells in one quarter, a skipped
#' quarter and a quarter in which every variable is fixed.
cf_conditions <- function() {
  x <- matrix(NA_real_, cf_H, 3, dimnames = list(NULL, c("a", "b", "c")))
  x[1, "a"] <- 2
  x[3, c("b", "c")] <- c(0, 1)
  x[5, ] <- c(1, 1, 1)
  x
}
cf_free <- function() matrix(NA_real_, cf_H, 3, dimnames = list(NULL, c("a", "b", "c")))

#' Several posterior "draws": the known VAR with its coefficients scaled.
cf_draws <- function(n) {
  beta <- array(NA_real_, c(n, dim(cf_B)))
  sigma <- array(NA_real_, c(n, dim(cf_Sigma)))
  for (i in seq_len(n)) {
    Bi <- cf_B
    Bi[-1, ] <- Bi[-1, ] * (1 - 0.05 * (i - 1))
    beta[i, , ] <- Bi
    sigma[i, , ] <- cf_Sigma * (1 + 0.1 * (i - 1))
  }
  list(beta = beta, sigma = sigma)
}

test_that("conditions are checked and the method must be named", {
  d <- cf_draws(1)
  expect_error(check_conditions(cf_conditions()[-1, ], cf_hist, cf_H), "rows")
  bad <- cf_conditions(); colnames(bad) <- c("b", "a", "c")
  expect_error(check_conditions(bad, cf_hist, cf_H), "same order")
  bad <- cf_conditions(); bad[2, 2] <- Inf
  expect_error(check_conditions(bad, cf_hist, cf_H), "Inf")
  expect_error(forecast_conditional(d$beta, d$sigma, cf_hist, 2, cf_H, cf_conditions(), seed = 1),
               "give the method")
  expect_error(forecast_conditional(d$beta, d$sigma, cf_hist, 2, cf_H, cf_conditions(), seed = 1,
                                    method = "other"), "unknown")
})

test_that("the building blocks match the companion form", {
  # Psi_j from the recursion equals the leading block of F^j.
  Psi <- ma_coefficients(cf_B, 2, 4)
  F_mat <- companion_matrix(cf_B, 2)
  expect_equal(Psi[[4]], (F_mat %*% F_mat %*% F_mat)[1:3, 1:3])
  # The state-space transition reproduces one step of the VAR.
  ss <- var_state_space(cf_B, cf_Sigma, cf_hist, 2)
  y1 <- cf_B[1, ] + cf_B1 %*% cf_hist[5, ] + cf_B2 %*% cf_hist[4, ]
  expect_equal(as.vector((ss$G %*% ss$S_T)[1:3]), as.vector(y1))
  expect_equal(ss$Q[1:3, 1:3], cf_Sigma)
})

test_that("with no conditions both methods give the unconditional forecast", {
  d <- cf_draws(4)
  unc <- forecast_unconditional(d$beta, d$sigma, cf_hist, 2, cf_H, seed = 3)
  unc_mean <- forecast_unconditional(d$beta, d$sigma, cf_hist, 2, cf_H, seed = 3, shocks = FALSE)
  for (m in cf_methods) {
    expect_equal(forecast_conditional(d$beta, d$sigma, cf_hist, 2, cf_H, cf_free(), seed = 3, method = m),
                 unc, tolerance = 1e-10)
    expect_equal(forecast_conditional(d$beta, d$sigma, cf_hist, 2, cf_H, cf_free(), seed = 3, method = m,
                                      shocks = FALSE), unc_mean, tolerance = 1e-10)
  }
})

test_that("conditional moments equal the textbook conditional normal", {
  cond <- cf_conditions()
  tb <- condition_textbook(forecast_joint_moments_exact(cf_B, cf_Sigma, cf_hist, 2, cf_H), cond)
  g <- conditional_moments_gaussian(cf_B, cf_Sigma, cf_hist, 2, cf_H, cond)
  k <- conditional_moments_kalman(cf_B, cf_Sigma, cf_hist, 2, cf_H, cond)
  expect_equal(unname(g$mean), unname(tb$mean), tolerance = 1e-10)
  expect_equal(unname(k$mean), unname(tb$mean), tolerance = 1e-10)
  expect_equal(g$cov, tb$cov, tolerance = 1e-10)
  for (h in seq_len(cf_H)) {
    expect_equal(g$var[[h]], tb$var[[h]], tolerance = 1e-10)
    expect_equal(k$var[[h]], tb$var[[h]], tolerance = 1e-10)
  }
  # A condition only in a later quarter: the earlier quarters still respond.
  late <- cf_free(); late[4, "b"] <- 3
  tb_late <- condition_textbook(forecast_joint_moments_exact(cf_B, cf_Sigma, cf_hist, 2, cf_H), late)
  k_late <- conditional_moments_kalman(cf_B, cf_Sigma, cf_hist, 2, cf_H, late)
  expect_equal(unname(k_late$mean), unname(tb_late$mean), tolerance = 1e-10)
})

test_that("fixed cells are met in every path and have zero variance", {
  cond <- cf_conditions()
  d <- cf_draws(3)
  fixed <- which(!is.na(cond))
  for (m in cf_methods) {
    paths <- forecast_conditional(d$beta, d$sigma, cf_hist, 2, cf_H, cond, seed = 4, method = m, n_paths = 50)
    for (i in seq_len(dim(paths)[1])) expect_lt(max(abs(paths[i, , ][fixed] - cond[fixed])), 1e-10)
  }
  # Conditional variances: zero at fixed cells, never above the unconditional ones.
  unc <- forecast_joint_moments_exact(cf_B, cf_Sigma, cf_hist, 2, cf_H)
  k <- conditional_moments_kalman(cf_B, cf_Sigma, cf_hist, 2, cf_H, cond)
  var_cond <- t(sapply(k$var, diag))                      # H x n
  var_unc <- matrix(diag(unc$cov), cf_H, 3, byrow = TRUE)
  expect_lt(max(abs(var_cond[fixed])), 1e-10)
  expect_true(all(var_cond <= var_unc + 1e-10))
})

test_that("with the same seed the two methods give the same paths", {
  d <- cf_draws(3)
  g <- forecast_conditional(d$beta, d$sigma, cf_hist, 2, cf_H, cf_conditions(), seed = 8,
                            method = "gaussian", n_paths = 20)
  k <- forecast_conditional(d$beta, d$sigma, cf_hist, 2, cf_H, cf_conditions(), seed = 8,
                            method = "kalman", n_paths = 20)
  expect_equal(g, k, tolerance = 1e-10)
})

test_that("simulated conditional paths match the exact conditional moments", {
  cond <- cf_conditions()
  tb <- condition_textbook(forecast_joint_moments_exact(cf_B, cf_Sigma, cf_hist, 2, cf_H), cond)
  b1 <- array(cf_B, c(1, dim(cf_B))); s1 <- array(cf_Sigma, c(1, dim(cf_Sigma)))
  n_sim <- 20000
  free <- which(is.na(cond))
  sd_exact <- sqrt(matrix(diag(tb$cov), cf_H, 3, byrow = TRUE))
  for (m in cf_methods) {
    paths <- forecast_conditional(b1, s1, cf_hist, 2, cf_H, cond, seed = 12, method = m, n_paths = n_sim)
    sim_mean <- apply(paths, c(2, 3), mean)
    sim_var  <- apply(paths, c(2, 3), var)
    z <- (sim_mean - tb$mean)[free] / (sd_exact[free] / sqrt(n_sim))
    ratio <- sim_var[free] / sd_exact[free]^2
    expect_lt(max(abs(z)), 4)
    expect_true(all(abs(ratio - 1) < 4 * sqrt(2 / (n_sim - 1))))
  }
})

# ---- Engine properties (R/engine_checks.R) ----

test_that("null scenario, linearity, superposition and ordering hold exactly", {
  base <- forecast_unconditional(array(cf_B, c(1, dim(cf_B))), array(cf_Sigma, c(1, dim(cf_Sigma))),
                                 cf_hist, 2, cf_H, seed = 1, shocks = FALSE)[1, , ]
  template <- cf_conditions()
  delta <- template - base                              # NA at free cells
  partner <- cf_free(); partner[2, "a"] <- 0; partner[5, "a"] <- 0
  partner_delta <- partner; partner_delta[2, "a"] <- -1; partner_delta[5, "a"] <- 0.5
  for (m in cf_methods) {
    res <- engine_properties(cf_B, cf_Sigma, cf_hist, 2, cf_H, template, delta, partner, partner_delta,
                             alphas = c(-1, 0.5, 3), perms = list(reversed = 3:1, rotated = c(2, 3, 1)),
                             method = m)
    expect_setequal(unique(res$property), c("null scenario", "linearity", "superposition", "ordering"))
    expect_true(all(res$rel_diff < 1e-10), info = m)
  }
})

test_that("permute_variables reorders a VAR consistently", {
  perm <- c(3, 1, 2)
  pv <- permute_variables(cf_B, cf_Sigma, cf_hist, cf_conditions(), perm, lags = 2)
  d0 <- forecast_unconditional(array(cf_B, c(1, dim(cf_B))), array(cf_Sigma, c(1, dim(cf_Sigma))),
                               cf_hist, 2, 4, seed = 1, shocks = FALSE)[1, , ]
  d1 <- forecast_unconditional(array(pv$B, c(1, dim(pv$B))), array(pv$Sigma, c(1, dim(pv$Sigma))),
                               pv$history, 2, 4, seed = 1, shocks = FALSE)[1, , ]
  expect_equal(d1, d0[, perm])
})

# ---- BVAR's conditional forecasts (cross-check; BVAR 1.0.5) ----

test_that("BVAR's conditional draws match the exact moments when two or more cells are fixed", {
  # BVAR requires at least one free variable in every quarter.
  cond <- cf_conditions()
  cond[5, "c"] <- NA
  p <- bvar_conditional_one_draw(cf_B, cf_Sigma, cf_hist, 2, cf_H, cond, n_paths = 4000, seed = 21)
  fixed <- which(!is.na(cond))
  expect_lt(max(apply(p, 1, function(x) max(abs(x[fixed] - cond[fixed])))), 1e-10)
  tb <- condition_textbook(forecast_joint_moments_exact(cf_B, cf_Sigma, cf_hist, 2, cf_H), cond)
  st <- paths_vs_exact_stats(p, tb, cond)
  expect_gt(st$chi2_p_value, 0.001)
  expect_lt(abs(st$var_ratio_max - 1), 0.15)
  # BVAR's Cholesky responses are Psi_{h-1} L.
  expect_lt(bvar_irf_agreement(cf_B, cf_Sigma, cf_hist, 2, cf_H), 1e-12)
})

test_that("BVAR 1.0.5 mishandles a single fixed cell (recorded behaviour)", {
  # Variable b in quarter 1: d = sqrt(Sigma_bb) = sqrt(0.5), 1 / d = 1.41, so
  # diag(1 / d) is the 1 x 1 identity and the cell moves by d r, not r.
  res <- bvar_single_cell(cf_B, cf_Sigma, cf_hist, 2, cf_H,
                          data.frame(variable = "b", quarter = 1, deviation = 1))
  expect_equal(res$bvar_deviation, res$d_times_r, tolerance = 1e-10)
  expect_gt(abs(res$bvar_deviation - res$requested), 0.2)
})
