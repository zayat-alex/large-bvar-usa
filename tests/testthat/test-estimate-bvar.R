# test-estimate-bvar.R: the estimator (R/estimate_bvar.R).
#
# Small simulated data and targets with known answers, so the tests are fast
# and need no data build.

set.seed(321)
est_y <- cbind(
  a = 100 + cumsum(rnorm(80, 0.5, 1)),
  b = 50 + cumsum(rnorm(80, 0.2, 0.5)),
  c = 20 + rnorm(80, 0, 3)
)
est_series <- data.frame(short_name = c("a", "b", "c"), prior = c("RW", "RW", "WN"),
                         stringsAsFactors = FALSE)
est_prior_cfg <- list(
  lambda = list(mode = 0.2, sd = 0.4, min = 0.0001, max = 5),
  alpha = 2,
  soc = list(mode = 1, sd = 1, min = 0.0001, max = 50),
  soc_ybar = "presample",
  constant_var = 1e7,
  psi = list(shape = 0.0004, scale = 0.0004, start = "ar_variance", bound_factor = 100)
)
est_setup <- estimation_setup(est_y, est_series, est_prior_cfg, lags = 2)

test_that("run_metropolis recovers a known normal target", {
  # Target: independent normals with means (1, -2) and sds (1, 2).
  target <- function(th) dnorm(th[[1]], 1, 1, log = TRUE) + dnorm(th[[2]], -2, 2, log = TRUE)
  W  <- diag(c(1, 4))
  mh <- run_metropolis(target, c(x = 0, y = 0), W, n_draw = 30000, n_burn = 5000, seed = 1)
  expect_equal(colMeans(mh$theta), c(x = 1, y = -2), tolerance = 0.1)
  expect_equal(unname(apply(mh$theta, 2, sd)), c(1, 2), tolerance = 0.1)
  expect_gt(mh$acceptance, 0.15)
  expect_lt(mh$acceptance, 0.5)
})

test_that("log_post_theta adds the Jacobian of theta = log(h)", {
  theta <- log(est_setup$start)
  expect_equal(log_post_theta(theta, est_setup) - log_post_hyper(exp(theta), est_setup),
               sum(theta), tolerance = 1e-8)
  # Outside the bounds the density is zero.
  theta_out <- theta
  theta_out[["lambda"]] <- log(10)
  expect_equal(log_post_theta(theta_out, est_setup), -Inf)
})

test_that("fixed hyperparameters are held and restored in order", {
  fixed <- c(psi1 = 2, psi3 = 5)
  s <- fix_hyperparameters(est_setup, fixed)
  expect_equal(names(s$start), c("lambda", "soc", "psi2"))
  full <- complete_hyper(c(lambda = 0.3, soc = 0.7, psi2 = 1), s)
  expect_equal(full, c(lambda = 0.3, soc = 0.7, psi1 = 2, psi2 = 1, psi3 = 5))
  expect_error(fix_hyperparameters(est_setup, c(nu = 1)), "unknown hyperparameter")
})

test_that("niw_posterior's B_hat equals the marginal likelihood code's", {
  h <- c(lambda = 0.3, soc = 0.5, psi1 = 1, psi2 = 0.3, psi3 = 9)
  post <- niw_posterior(h, est_setup)
  dm  <- var_data_matrices(est_y, 2)
  dum <- soc_dummy_obs(est_setup$ybar, 0.5, 2)
  ml  <- log_ml_niw(rbind(dum$Y, dm$Y), rbind(dum$X, dm$X), prior_mean_matrix(est_setup$b_own, 2),
                    minnesota_omega(0.3, 2, c(1, 0.3, 9), 2, 1e7), c(1, 0.3, 9), d = 5)
  expect_equal(post$B_hat, ml$B_hat, tolerance = 1e-8)
  # Degrees of freedom: 78 data rows + 3 dummy rows + n + 2.
  expect_equal(post$dof, 78 + 3 + 3 + 2)
})

test_that("draw_coefficients matches the exact posterior moments", {
  h <- c(lambda = 0.3, soc = 0.5, psi1 = 1, psi2 = 0.3, psi3 = 9)
  post <- niw_posterior(h, est_setup)
  set.seed(5)
  d <- draw_coefficients(post, 3000)
  # Mean of B within 4.5 Monte Carlo standard errors of B_hat, every element.
  z <- (apply(d$beta, c(2, 3), mean) - post$B_hat) / (apply(d$beta, c(2, 3), sd) / sqrt(3000))
  expect_lt(max(abs(z)), 4.5)
  # Mean of Sigma against (Psi + S) / (dof - n - 1).
  sigma_exact <- post$scale_matrix / (post$dof - 3 - 1)
  expect_equal(apply(d$sigma, c(2, 3), mean), sigma_exact, tolerance = 0.03)
})

test_that("find_hyper_mode reaches the same mode from every start", {
  m <- find_hyper_mode(est_setup, start_factors = c(1, 0.5, 2))
  expect_true(all(m$table$convergence == 0))
  hcols <- names(est_setup$start)
  spread <- apply(m$table[, hcols], 2, function(x) (max(x) - min(x)) / mean(x))
  expect_lt(max(spread), 1e-3)
})

test_that("mode_objective maps a point one rounding step past a bound onto the bound", {
  theta <- log(est_setup$start)
  theta[["soc"]] <- log(est_setup$upper[["soc"]])
  past <- theta
  past[["soc"]] <- past[["soc"]] + 4.4e-16 * 10      # just outside, as optim can evaluate
  expect_equal(log_post_theta(past, est_setup), -Inf)  # the sampler still rejects it
  expect_equal(mode_objective(past, est_setup), log_post_theta(theta, est_setup))
  # Inside the bounds the objective is log_post_theta itself.
  inside <- log(est_setup$start)
  expect_identical(mode_objective(inside, est_setup), log_post_theta(inside, est_setup))
})
