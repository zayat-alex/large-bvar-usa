# test-full-size-checks.R: the companion matrix (R/estimate_run.R).

test_that("companion_matrix reproduces a VAR(2) with known coefficients", {
  # y_t = B1 y_{t-1} + B2 y_{t-2}, n = 2. In B (K x n), row 1 is the constant,
  # rows 2-3 hold lag 1 of variables 1 and 2, rows 4-5 lag 2, so B_s = t(B[rows, ]).
  B1 <- matrix(c(0.5, 0.1, 0.0, 0.3), 2)
  B2 <- matrix(c(0.2, 0.0, 0.05, 0.1), 2)
  B  <- rbind(c(0, 0), t(B1), t(B2))
  F_mat <- companion_matrix(B, 2)
  z <- c(1, -1, 0.5, 0.2)                     # (y_{t-1}', y_{t-2}')'
  expect_equal(as.vector(F_mat %*% z), c(as.vector(B1 %*% z[1:2] + B2 %*% z[3:4]), z[1:2]))

  # The eigenvalues of F are the roots of det(I x^2 - B1 x - B2) = 0, a quartic
  # whose coefficients follow from the 2 x 2 determinant.
  det_poly <- function(x) det(diag(2) * x^2 - B1 * x - B2)
  xs <- c(-2, -1, 0, 1, 2)
  coefs <- solve(outer(xs, 0:4, `^`), sapply(xs, det_poly))
  roots <- polyroot(coefs)
  expect_equal(companion_max_modulus(B, 2), max(Mod(roots)), tolerance = 1e-10)
})

test_that("companion_max_modulus is 1 for a random walk", {
  B <- rbind(c(0, 0), diag(2), matrix(0, 2, 2))
  expect_equal(companion_max_modulus(B, 2), 1)
})
