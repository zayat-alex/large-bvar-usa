# estimate_ml.R: the marginal likelihood of the data and the posterior of
# the hyperparameters, written directly from the formulas.
#
# Why it exists: the hyperparameters (lambda, mu, psi) are chosen by the
# marginal likelihood of the data (Giannone, Lenza and Primiceri 2015, Section
# 2 and Appendix A), and every posterior draw depends on it. The estimator
# (R/estimate_bvar.R) uses these functions. They share no code with the BVAR
# package, whose value they were checked against (R/estimate_checks.R,
# bvar_log_posterior()).
#
# Notation (as in Giannone, Lenza and Primiceri 2015):
#   Y (N x n), X (N x K): left-hand-side data and regressors, X = [1, y_{t-1}, ..., y_{t-p}]
#   b (K x n): prior mean of the coefficients
#   Omega = diag(omega) (K x K): prior covariance factor of the coefficients
#   Psi = diag(psi) (n x n): scale of the inverse-Wishart prior on Sigma
#   d = n + 2: inverse-Wishart degrees of freedom
#
# Marginal likelihood of the normal-inverse-Wishart VAR:
#   p(Y) = pi^(-nN/2) * Gamma_n((N+d)/2) / Gamma_n(d/2)
#          * |Omega|^(-n/2) * |X'X + Omega^-1|^(-n/2)
#          * |Psi|^(d/2) * |Psi + S|^(-(N+d)/2)
#   B_hat = (X'X + Omega^-1)^-1 (X'Y + Omega^-1 b)
#   S     = (Y - X B_hat)'(Y - X B_hat) + (B_hat - b)' Omega^-1 (B_hat - b)
# With dummy-observation priors: p(Y) = p([Y+; Y]) / p(Y+).
#
# The formulas and settings are set out in docs/methodology.md (sections 4-5);
# code comments cite its equation labels (M1), (M2), ...
#
# Inputs:  numeric data matrices and hyperparameter values.
# Outputs: log densities, returned with their components for inspection.


# ---- Data matrices ----

#' Left-hand-side and regressor matrices of a VAR(p) with a constant.
#'
#' Implements docs/methodology.md (M2).
#'
#' Excluding quarters (`keep`): the equation of quarter t is the row
#' (y_t, x_t). Excluding quarter t drops that row only; the quarter's data
#' still enter later rows as lags. This is how the COVID treatment
#' `exclude_window` leaves the pandemic quarters out of the estimation
#' (docs/methodology.md section 11).
#'
#' @param y Numeric matrix T x n (one column per variable).
#' @param lags Number of lags p.
#' @param keep NULL (every equation), or a logical vector of length T: TRUE for
#'   the quarters whose equation is used. Its first p entries are ignored
#'   (those quarters only serve as initial lags).
#' @return List: Y (N x n, rows p+1..T that are kept), X (N x K, columns:
#'   constant, then lag 1 of all variables, lag 2 of all variables, ..., lag p).
var_data_matrices <- function(y, lags, keep = NULL) {
  n_total <- nrow(y)
  n_var   <- ncol(y)
  rows    <- (lags + 1):n_total
  if (!is.null(keep)) {
    if (length(keep) != n_total) stop("keep must have one entry per row of y")
    rows <- rows[keep[rows]]
  }

  # Step 1: left-hand side.
  Y <- y[rows, , drop = FALSE]

  # Step 2: regressors, one block of n columns per lag.
  X <- matrix(1, nrow = length(rows), ncol = 1 + n_var * lags)
  for (s in seq_len(lags)) {
    cols <- 1 + (s - 1) * n_var + seq_len(n_var)
    X[, cols] <- y[rows - s, , drop = FALSE]
  }
  list(Y = unname(Y), X = X)
}


# ---- Prior pieces ----

#' Diagonal of Omega, the prior covariance factor of the coefficients.
#'
#' Constant: constant_var. Variable j at lag s: lambda^2 / (s^alpha * psi_j).
#'
#' Implements docs/methodology.md (M8).
#'
#' @return Numeric vector of length K = 1 + n p.
minnesota_omega <- function(lambda, alpha, psi, lags, constant_var) {
  omega <- constant_var
  for (s in seq_len(lags)) {
    omega <- c(omega, lambda^2 / (s^alpha * psi))
  }
  omega
}

#' Prior mean matrix b: own first lag = b_own, everything else 0.
#'
#' Implements docs/methodology.md (M6).
#'
#' @param b_own Numeric vector (length n): 1 for random walk, 0 for white noise.
#' @param lags Number of lags p.
#' @return Numeric matrix K x n.
prior_mean_matrix <- function(b_own, lags) {
  n_var <- length(b_own)
  b <- matrix(0, nrow = 1 + n_var * lags, ncol = n_var)
  for (j in seq_len(n_var)) {
    b[1 + j, j] <- b_own[j]       # row 1 + j: lag 1 of variable j
  }
  b
}

#' Sum-of-coefficients dummy observations.
#'
#' Y+ = diag(ybar) / mu; X+ = [0, Y+, ..., Y+] (Y+ repeated p times).
#'
#' Implements docs/methodology.md (M9).
#'
#' @param ybar Numeric vector (length n): the centring means.
#' @param mu Tightness of the prior.
#' @param lags Number of lags p.
#' @return List: Y (n x n), X (n x K).
soc_dummy_obs <- function(ybar, mu, lags) {
  y_plus <- diag(ybar, nrow = length(ybar)) / mu
  x_plus <- cbind(0, do.call(cbind, rep(list(y_plus), lags)))
  list(Y = y_plus, X = x_plus)
}


# ---- Marginal likelihood ----

#' Log of the multivariate gamma function Gamma_n(a).
log_mv_gamma <- function(n, a) {
  n * (n - 1) / 4 * log(pi) + sum(lgamma(a + (1 - seq_len(n)) / 2))
}

#' Log determinant of a symmetric positive definite matrix.
log_det <- function(m) {
  as.numeric(determinant(m, logarithm = TRUE)$modulus)
}

#' Log marginal likelihood of a normal-inverse-Wishart VAR.
#'
#' Implements docs/methodology.md (M12)-(M14).
#'
#' @param Y,X Data matrices (N x n, N x K).
#' @param b Prior mean, K x n.
#' @param omega Diagonal of Omega (length K).
#' @param psi Diagonal of Psi (length n).
#' @param d Inverse-Wishart degrees of freedom.
#' @param beta_at "posterior" (the formula above) or "prior_mean" (B_hat
#'   replaced by b and the penalty term dropped, the shortcut BVAR uses for the
#'   dummy observations; kept here to reproduce it).
#' @return List with each term of the formula and their sum, `log_ml`.
log_ml_niw <- function(Y, X, b, omega, psi, d, beta_at = "posterior") {
  n_var <- ncol(Y)
  n_obs <- nrow(Y)
  omega_inv <- diag(1 / omega, nrow = length(omega))

  # Step 1: the precision matrix in rescaled form.
  # X'X + Omega^-1 can be too ill-conditioned to invert directly (the constant
  # has a diffuse prior, omega = 1e7, and with a tight sum-of-coefficients
  # prior the dummy rows are of order ybar / mu). The algebraically equal form
  #   X'X + Omega^-1 = Omega^-1/2 (I + Omega^1/2 X'X Omega^1/2) Omega^-1/2
  # has a middle matrix M with all eigenvalues >= 1, which is stable to factor.
  # M is factored by a symmetric eigendecomposition of Omega^1/2 X'X Omega^1/2,
  # which is positive semi-definite in exact arithmetic; rounding can make its
  # smallest eigenvalues slightly negative (or break a Cholesky factorization)
  # when its entries are very large, so negative eigenvalues are set to 0.
  omega_half <- diag(sqrt(omega), nrow = length(omega))
  G     <- omega_half %*% crossprod(X) %*% omega_half
  eig   <- eigen((G + t(G)) / 2, symmetric = TRUE)
  ev    <- pmax(eig$values, 0) + 1                      # eigenvalues of M, all >= 1
  M_inv <- eig$vectors %*% diag(1 / ev, nrow = length(ev)) %*% t(eig$vectors)
  log_det_M <- sum(log(ev))

  # Step 2: coefficient estimate and residual cross-product S.
  if (beta_at == "posterior") {
    # B_hat = (X'X + Omega^-1)^-1 (X'Y + Omega^-1 b) = b + Omega^1/2 M^-1 Omega^1/2 X'(Y - X b)
    B_hat <- b + omega_half %*% M_inv %*% omega_half %*% crossprod(X, Y - X %*% b)
    E <- Y - X %*% B_hat
    S <- crossprod(E) + t(B_hat - b) %*% omega_inv %*% (B_hat - b)
  } else if (beta_at == "prior_mean") {
    E <- Y - X %*% b
    S <- crossprod(E)
  } else {
    stop("beta_at must be posterior or prior_mean")
  }

  # Step 3: the terms of the log marginal likelihood.
  #   -n/2 log|Omega| - n/2 log|X'X + Omega^-1| = -n/2 log|M|
  #   d/2 log|Psi| - (N+d)/2 log|Psi + S| = -N/2 log|Psi| - (N+d)/2 log|I + Psi^-1/2 S Psi^-1/2|
  psi_inv_half <- diag(1 / sqrt(psi), nrow = n_var)
  terms <- c(
    pi_term     = -n_var * n_obs / 2 * log(pi),
    gamma_ratio = log_mv_gamma(n_var, (n_obs + d) / 2) - log_mv_gamma(n_var, d / 2),
    precision   = -n_var / 2 * log_det_M,
    psi_det     = -n_obs / 2 * sum(log(psi)),
    fit         = -(n_obs + d) / 2 * log_det(diag(n_var) + psi_inv_half %*% S %*% psi_inv_half)
  )
  list(terms = terms, log_ml = sum(terms), B_hat = if (beta_at == "posterior") B_hat else b)
}


# ---- Hyperpriors ----

#' Shape and scale of a Gamma distribution with a given mode and sd.
#'
#' Solves mode = (k - 1) theta and sd^2 = k theta^2 for k > 1:
#' k = (2 + r + sqrt(r (r + 4))) / 2 with r = mode^2 / sd^2; theta = mode / (k - 1).
#'
#' Implements docs/methodology.md section 4.4.
gamma_shape_scale <- function(mode, sd) {
  r <- mode^2 / sd^2
  k <- (2 + r + sqrt(r * (r + 4))) / 2
  list(shape = k, scale = mode / (k - 1))
}

#' Log inverse-gamma density: a log(beta) - lgamma(a) - (a + 1) log(x) - beta / x.
log_dinvgamma <- function(x, shape, scale) {
  shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - scale / x
}


# ---- Hyperparameter posterior ----

#' Log posterior kernel of the hyperparameters, computed directly.
#'
#' log p(gamma | Y) = log p(Y | gamma) + log p(lambda) + log p(mu) + sum_j log p(psi_j)
#'
#' Implements docs/methodology.md (M10), (M15).
#'
#' @param y Numeric estimation matrix T x n.
#' @param lags Number of lags p.
#' @param lambda,mu,psi Hyperparameter values (psi: length n).
#' @param b_own Prior mean of each own first lag (length n).
#' @param prior_cfg The `prior` section of the run configuration (hyperpriors,
#'   alpha, constant_var, psi shape and scale).
#' @param ybar Centring means of the sum-of-coefficients dummies (length n).
#' @param dummy_beta_at "posterior" (correct) or "prior_mean" (BVAR's shortcut).
#' @param keep NULL, or the equations to use (var_data_matrices()).
#' @return List: log_ml_augmented, log_ml_dummy, log_ml (their difference),
#'   log_hyperprior (by hyperparameter), log_posterior (total).
log_posterior_direct <- function(y, lags, lambda, mu, psi, b_own, prior_cfg, ybar,
                                 dummy_beta_at = "posterior", keep = NULL) {
  n_var <- ncol(y)
  d     <- n_var + 2

  # Step 1: data matrices, prior mean, Omega, dummy observations.
  dm    <- var_data_matrices(y, lags, keep)
  b     <- prior_mean_matrix(b_own, lags)
  omega <- minnesota_omega(lambda, prior_cfg$alpha, psi, lags, prior_cfg$constant_var)
  dum   <- soc_dummy_obs(ybar, mu, lags)

  # Step 2: marginal likelihood of the data augmented with the dummies,
  # divided by that of the dummies alone.
  aug   <- log_ml_niw(rbind(dum$Y, dm$Y), rbind(dum$X, dm$X), b, omega, psi, d)
  dummy <- log_ml_niw(dum$Y, dum$X, b, omega, psi, d, beta_at = dummy_beta_at)
  log_ml <- aug$log_ml - dummy$log_ml

  # Step 3: log hyperprior densities.
  g_lam <- gamma_shape_scale(prior_cfg$lambda$mode, prior_cfg$lambda$sd)
  g_mu  <- gamma_shape_scale(prior_cfg$soc$mode, prior_cfg$soc$sd)
  log_hyperprior <- c(
    lambda = dgamma(lambda, shape = g_lam$shape, scale = g_lam$scale, log = TRUE),
    mu     = dgamma(mu, shape = g_mu$shape, scale = g_mu$scale, log = TRUE),
    psi    = sum(log_dinvgamma(psi, prior_cfg$psi$shape, prior_cfg$psi$scale))
  )

  list(log_ml_augmented = aug$log_ml, log_ml_dummy = dummy$log_ml, log_ml = log_ml,
       log_hyperprior = log_hyperprior,
       log_posterior = log_ml + sum(log_hyperprior))
}
