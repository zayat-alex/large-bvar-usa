# test-prior-sensitivity.R: the variants of the hyperparameter sensitivity
# check (R/prior_sensitivity_checks.R). The check on the real data is
# scripts/check_prior_sensitivity.R.

test_that("sensitivity_variants moves one hyperparameter at a time by the stated amount", {
  theta <- c(lambda = log(0.4), soc = log(0.3), psi1 = 0.5, psi2 = 1.5)
  sd_theta <- c(lambda = 0.06, soc = 0.12, psi1 = 0.2, psi2 = 0.2)
  v <- sensitivity_variants(theta, sd_theta)
  info <- attr(v, "info")
  expect_equal(names(v), info$variant)
  expect_equal(exp(v[["lambda x 2"]][["lambda"]]), 0.8)
  expect_equal(exp(v[["mu x 0.5"]][["soc"]]), 0.15)
  expect_equal(v[["psi x 2"]][c("psi1", "psi2")], theta[c("psi1", "psi2")] + log(2))
  expect_equal(v[["mu +2 sd"]][["soc"]], theta[["soc"]] + 0.24)
  expect_equal(v[["lambda -2 sd"]][["lambda"]], theta[["lambda"]] - 0.12)
  # Everything else stays at the mode, and the info table states the change.
  for (k in seq_along(v)) {
    changed <- names(theta)[v[[k]] != theta]
    expect_true(length(changed) >= 1)
    expect_equal(unique(round(v[[k]][changed] - theta[changed], 12)), round(info$theta_change[k], 12))
  }
})
