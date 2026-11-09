# test-estimate-exclusion.R: estimation with some quarters' equations left out
# (COVID treatment exclude_window; R/estimate_ml.R, R/estimate_bvar.R,
# R/estimate_priors.R). Small simulated data with known answers.

set.seed(321)
ex_y <- cbind(
  a = 100 + cumsum(rnorm(60, 0.5, 1)),
  b = 50 + cumsum(rnorm(60, 0.2, 0.5)),
  c = 20 + rnorm(60, 0, 3)
)
ex_series <- data.frame(short_name = c("a", "b", "c"), prior = c("RW", "RW", "WN"),
                        stringsAsFactors = FALSE)
ex_prior_cfg <- list(
  lambda = list(mode = 0.2, sd = 0.4, min = 0.0001, max = 5),
  alpha = 2,
  soc = list(mode = 1, sd = 1, min = 0.0001, max = 50),
  soc_ybar = "presample",
  constant_var = 1e7,
  psi = list(shape = 0.0004, scale = 0.0004, start = "ar_variance", bound_factor = 100)
)
ex_keep <- rep(TRUE, 60)
ex_keep[31:36] <- FALSE                          # six quarters' equations left out

test_that("var_data_matrices drops exactly the excluded equations and keeps their lags", {
  full <- var_data_matrices(ex_y, 2)
  part <- var_data_matrices(ex_y, 2, ex_keep)
  used <- ex_keep[3:60]                           # equation rows are quarters 3..60
  expect_equal(part$Y, full$Y[used, , drop = FALSE])
  expect_equal(part$X, full$X[used, , drop = FALSE])
  # The first equation after the gap (quarter 37) has quarters 36 and 35 as lags.
  k <- which(which(used) + 2 == 37)
  expect_equal(part$X[k, 2:4], unname(ex_y[36, ]))
  expect_equal(part$X[k, 5:7], unname(ex_y[35, ]))
  # With nothing excluded, the result is unchanged.
  expect_identical(var_data_matrices(ex_y, 2, rep(TRUE, 60)), full)
})

test_that("the AR residual variance uses only the kept equations", {
  x <- ex_y[, "a"]
  v <- ar_residual_variance(x, 2, ex_keep)
  # By hand: regression of x_t on (1, x_{t-1}, x_{t-2}) over the kept t.
  t_used <- (3:60)[ex_keep[3:60]]
  X <- cbind(1, x[t_used - 1], x[t_used - 2])
  r <- x[t_used] - X %*% solve(crossprod(X), crossprod(X, x[t_used]))
  expect_equal(v$variance, sum(r^2) / (length(t_used) - 3))
  expect_equal(v$n_obs, length(t_used))
  expect_equal(ar_residual_variance(x, 2), ar_residual_variance(x, 2, rep(TRUE, 60)))
})

test_that("setup, posterior and log posterior with excluded equations", {
  s_full <- estimation_setup(ex_y, ex_series, ex_prior_cfg, 2)
  s_part <- estimation_setup(ex_y, ex_series, ex_prior_cfg, 2, keep = ex_keep)
  expect_null(s_full$keep)                        # unchanged for runs without exclusions
  expect_identical(s_part$keep, ex_keep)
  expect_identical(s_part$ybar, s_full$ybar)      # centring on the presample, unaffected

  # The posterior on the kept equations equals the posterior computed from
  # the data matrices with those rows removed by hand.
  h <- complete_hyper(s_part$start, s_part)
  post <- niw_posterior(h, s_part)
  dm <- var_data_matrices(ex_y, 2)
  used <- ex_keep[3:60]
  dum <- soc_dummy_obs(s_part$ybar, h[["soc"]], 2)
  Y <- rbind(dum$Y, dm$Y[used, ]); X <- rbind(dum$X, dm$X[used, ])
  b <- prior_mean_matrix(s_part$b_own, 2)
  omega <- minnesota_omega(h[["lambda"]], 2, h[grep("^psi", names(h))], 2, 1e7)
  B_hand <- solve(crossprod(X) + diag(1 / omega), crossprod(X, Y) + diag(1 / omega) %*% b)
  expect_equal(post$B_hat, B_hand, tolerance = 1e-8)
  expect_equal(post$dof, nrow(Y) + 3 + 2)

  # Log posterior: the same as computing it on the kept rows directly.
  lp <- log_post_inside(h, s_part)
  lp_direct <- log_posterior_direct(ex_y, 2, h[["lambda"]], h[["soc"]], h[grep("^psi", names(h))],
                                    s_part$b_own, ex_prior_cfg, s_part$ybar, keep = ex_keep)$log_posterior
  expect_equal(lp, lp_direct)
  expect_false(isTRUE(all.equal(lp, log_post_inside(h, s_full))))
})
