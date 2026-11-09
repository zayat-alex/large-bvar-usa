# forecast.R: unconditional forecasts from posterior draws of the VAR
# (the baseline of the forecasting engine).
#
# Why it exists: a forecast from a Bayesian VAR is a sample from the
# predictive density, which mixes two kinds of uncertainty: about the
# coefficients (one forecast path per posterior draw of B and Sigma) and
# about future shocks (random innovations fed through the model's dynamics).
# Crump et al. (2021), Section 2.3, use this predictive density. For each
# draw m of (B, Sigma), starting from the last p observed quarters:
#
#   y_{T+h} = c + B_1 y_{T+h-1} + ... + B_p y_{T+h-p} + e_{T+h},
#   e_{T+h} = L u_{T+h},  L L' = Sigma (lower Cholesky factor),  u ~ N(0, I),
#
# for h = 1, ..., H, where a lagged value beyond T is the simulated value of
# the same path. Because each simulated shock enters the next quarter's lags,
# uncertainty accumulates over the horizon. The forecast can start from any
# history (for example the latest data when the coefficients were estimated
# on a shorter sample), since only its last p rows are used.
#
# Notation: n variables, p lags, K = 1 + n p. B (K x n) stacks the
# coefficients as in the estimation code: row 1 is the constant c, and row
# 1 + (s - 1) n + j holds the coefficients on lag s of variable j, so that
# y_t' = x_t' B with x_t = (1, y_{t-1}', ..., y_{t-p}')'.
#
# The formulas and settings are set out in docs/methodology.md (section 9);
# code comments cite its equation labels (M1), (M2), ...
#
# Conditional forecasts (scenarios) are in R/forecast_conditional.R.
#
# Inputs:  coefficient draws (draws x K x n), covariance draws (draws x n x n),
#          a history matrix (rows = quarters, columns = variables in the
#          series order), the horizon, a seed.
# Outputs: forecast draws, an array draws x H x n, in model units.


#' Unconditional forecast paths, one per posterior draw.
#'
#' Implements docs/methodology.md (M18).
#'
#' @param beta Coefficient draws, array draws x K x n.
#' @param sigma Covariance draws, array draws x n x n.
#' @param history Numeric matrix of observed data (quarters x n), in the
#'   variable order of the draws; its last p rows are the starting point.
#' @param lags Number of lags p.
#' @param horizon Number of quarters H to forecast.
#' @param seed Seed of the random number generator (shocks).
#' @param shocks TRUE to simulate shocks (a draw from the predictive density);
#'   FALSE for the path with all future shocks at zero (the mean of the
#'   forecast given the draw's coefficients).
#' @return Array draws x H x n, with the variable names of `history`.
forecast_unconditional <- function(beta, sigma, history, lags, horizon, seed, shocks = TRUE) {
  n_draws <- dim(beta)[1]
  K       <- dim(beta)[2]
  n_var   <- dim(beta)[3]
  if (K != 1 + n_var * lags) stop("beta has ", K, " rows per draw; expected 1 + n * lags = ", 1 + n_var * lags)
  if (ncol(history) != n_var) stop("history has ", ncol(history), " columns; the draws have ", n_var, " variables")
  if (nrow(history) < lags) stop("history needs at least ", lags, " rows")
  if (anyNA(tail(history, lags))) stop("the last ", lags, " rows of history contain missing values")

  out <- array(NA_real_, c(n_draws, horizon, n_var),
               dimnames = list(NULL, paste0("h", seq_len(horizon)), colnames(history)))
  # Starting lags: row s of `start` is y_{T+1-s} (most recent quarter first).
  start <- history[nrow(history) - seq_len(lags) + 1, , drop = FALSE]
  set.seed(seed)

  for (m in seq_len(n_draws)) {
    B <- beta[m, , ]
    L <- if (shocks) t(chol(sigma[m, , ])) else NULL
    lagged <- start                                  # p x n, row s = lag s
    for (h in seq_len(horizon)) {
      # Step 1: regressors x = (1, y_{t-1}', ..., y_{t-p}')'.
      x <- c(1, as.vector(t(lagged)))
      # Step 2: conditional mean y_t = B' x.
      y <- as.vector(crossprod(B, x))
      # Step 3: add a shock with covariance Sigma.
      if (shocks) y <- y + as.vector(L %*% rnorm(n_var))
      out[m, h, ] <- y
      # Step 4: shift the lags: the new value becomes lag 1.
      lagged <- rbind(y, lagged[-lags, , drop = FALSE])
    }
  }
  out
}
