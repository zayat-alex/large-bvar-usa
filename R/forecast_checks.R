# forecast_checks.R: checks of the forecasting code (R/forecast.R) and of the
# BVAR package's forecasting functions.
#
# Why it exists: before the forecasting engine is built, four questions are
# answered with numbers:
#   1. Does the project's simulator produce the right predictive moments?
#      Checked on a VAR with known coefficients, where the mean and the
#      covariance of the h-step forecast have exact formulas.
#   2. How do BVAR::predict()'s forecasts compare when it is given the same
#      draws? Its mean paths should agree; its bands add noise with
#      covariance Sigma^2 to each quarter separately instead of feeding
#      shocks with covariance Sigma through the dynamics (BVAR GitHub issue
#      #86), so they should not.
#   3. Can forecasts start from a history other than the estimation sample
#      (needed when the coefficients are estimated on data up to 2019)?
#   4. How long do BVAR's conditional forecasts (scenarios) take at full size?
# The results are written by scripts/check_forecast.R and shown in
# reports/forecast_checks.qmd.
#
# Inputs:  coefficient and covariance draws, data matrices.
# Outputs: data.frames of comparisons and timings, figures.


# ---- Exact moments of a VAR forecast ----

#' Exact mean and covariance of the h-step forecast of a VAR with known
#' coefficients.
#'
#' Mean: the recursion of forecast_unconditional() with all shocks at zero.
#' Covariance (forecast error variance): with Psi_i the n x n block of F^i
#' (F the companion matrix; Psi_0 = I),
#'   Var(y_{T+h} | data) = sum_{i=0}^{h-1} Psi_i Sigma Psi_i'.
#'
#' @param B Coefficients, K x n (layout as in R/forecast.R).
#' @param Sigma Innovation covariance, n x n.
#' @param history Data matrix whose last p rows start the forecast.
#' @param lags Number of lags p.
#' @param horizon Number of quarters H.
#' @return List: mean (H x n), cov (list of H matrices n x n).
forecast_moments_exact <- function(B, Sigma, history, lags, horizon) {
  n_var <- ncol(B)
  beta1 <- array(B, c(1, dim(B)))
  sigma1 <- array(Sigma, c(1, dim(Sigma)))
  mean_path <- forecast_unconditional(beta1, sigma1, history, lags, horizon, seed = 1,
                                      shocks = FALSE)[1, , ]
  F_mat <- companion_matrix(B, lags)
  cov_list <- vector("list", horizon)
  acc   <- matrix(0, n_var, n_var)
  F_pow <- diag(n_var * lags)                      # F^0
  for (h in seq_len(horizon)) {
    Psi <- F_pow[seq_len(n_var), seq_len(n_var), drop = FALSE]   # Psi_{h-1}
    acc <- acc + Psi %*% Sigma %*% t(Psi)
    cov_list[[h]] <- acc
    F_pow <- F_mat %*% F_pow
  }
  list(mean = matrix(mean_path, horizon, n_var), cov = cov_list)
}

#' Simulated forecast moments against the exact ones, for a known VAR.
#'
#' The same (B, Sigma) is repeated n_sim times as if it were n_sim posterior
#' draws, so the spread of the simulated paths comes from the shocks only and
#' must match the exact forecast error covariance.
#'
#' @return data.frame: horizon, variable, exact and simulated mean and s.d.,
#'   and z = (simulated mean - exact mean) / (Monte Carlo s.e.), and the
#'   ratio of simulated to exact variance with its Monte Carlo s.e.
check_known_var <- function(B, Sigma, history, lags, horizon, n_sim, seed) {
  n_var <- ncol(B)
  beta  <- array(rep(B, each = n_sim), c(n_sim, dim(B)))
  sigma <- array(rep(Sigma, each = n_sim), c(n_sim, dim(Sigma)))
  sim   <- forecast_unconditional(beta, sigma, history, lags, horizon, seed)
  ex    <- forecast_moments_exact(B, Sigma, history, lags, horizon)
  rows <- list()
  for (h in seq_len(horizon)) for (j in seq_len(n_var)) {
    x  <- sim[, h, j]
    sd_exact <- sqrt(ex$cov[[h]][j, j])
    rows[[length(rows) + 1]] <- data.frame(
      horizon = h, variable = colnames(history)[j],
      mean_exact = ex$mean[h, j], mean_sim = mean(x),
      z_mean = (mean(x) - ex$mean[h, j]) / (sd_exact / sqrt(n_sim)),
      sd_exact = sd_exact, sd_sim = sd(x),
      var_ratio = var(x) / sd_exact^2,
      var_ratio_mc_se = sqrt(2 / (n_sim - 1)))
  }
  # Cross-covariance at the last horizon (checks that shocks are correlated through L).
  cov_sim <- cov(sim[, horizon, ])
  list(table = do.call(rbind, rows),
       cov_last_exact = ex$cov[[horizon]], cov_last_sim = cov_sim)
}


# ---- BVAR::predict() given the project's draws ----

#' A minimal BVAR fit object holding the project's draws.
#'
#' BVAR::predict() (version 1.0.5, read in source) uses only these elements:
#' beta (draws x K x n), sigma (draws x n x n), variables, and meta$n_save,
#' K, M, lags, Y, N, where the forecast starts from the last `lags` rows of
#' Y. Building the object directly lets predict() be run on exactly the
#' draws the project's simulator uses.
#'
#' @param beta,sigma Draws.
#' @param history Data matrix (quarters x n); its last p rows start forecasts.
#' @param lags Number of lags p.
#' @return An object of class "bvar".
bvar_shell <- function(beta, sigma, history, lags) {
  n_var <- dim(beta)[3]
  structure(list(beta = beta, sigma = sigma, variables = colnames(history),
                 meta = list(n_save = dim(beta)[1], K = dim(beta)[2], M = n_var, lags = lags,
                             Y = unname(history), N = nrow(history))),
            class = "bvar")
}

#' BVAR::predict() forecast draws (unconditional or conditional).
#'
#' @param shell Output of bvar_shell().
#' @param horizon Quarters.
#' @param seed Seed set immediately before the call.
#' @param newdata Optional history passed as `newdata`.
#' @param cond_path,cond_vars Optional conditions (as in BVAR::bv_fcast()).
#' @return Array draws x H x n.
bvar_predict_draws <- function(shell, horizon, seed, newdata = NULL, cond_path = NULL,
                               cond_vars = NULL) {
  set.seed(seed)
  spec <- BVAR::bv_fcast(horizon = horizon, cond_path = cond_path, cond_vars = cond_vars)
  fc <- if (is.null(newdata)) {
    suppressMessages(predict(shell, spec))
  } else {
    suppressMessages(predict(shell, spec, newdata = newdata))
  }
  fc$fcast
}

#' Spread of forecast draws at selected horizons: own simulator against
#' BVAR::predict(), with the same coefficient draws.
#'
#' @param own,pkg Arrays draws x H x n.
#' @param variables Names of the variables to report.
#' @param horizons Horizons to report.
#' @return data.frame: variable, horizon, mean and s.d. of each, and the
#'   ratio of BVAR's s.d. to the project's.
compare_spread <- function(own, pkg, variables, horizons) {
  rows <- list()
  for (v in variables) for (h in horizons) {
    j <- match(v, dimnames(own)[[3]])
    rows[[length(rows) + 1]] <- data.frame(
      variable = v, horizon = h,
      mean_own = mean(own[, h, j]), mean_bvar = mean(pkg[, h, j]),
      sd_own = sd(own[, h, j]), sd_bvar = sd(pkg[, h, j]),
      sd_ratio_bvar_to_own = sd(pkg[, h, j]) / sd(own[, h, j]))
  }
  do.call(rbind, rows)
}


# ---- Conditional forecasts: cost ----

#' Time BVAR's conditional forecast at full size.
#'
#' A scenario fixes `n_cond_vars` variables for the first `n_cond_q`
#' quarters at the mean (no-shock) path of the first draw plus a constant
#' shift, so the conditions are attainable. predict() first computes the
#' Cholesky impulse responses of every draw, then, per draw, builds the
#' restriction matrix (one row per fixed cell, n H columns) and takes its
#' singular value decomposition (Waggoner and Zha 1999).
#'
#' @param shell Output of bvar_shell() (its number of draws sets the timing sample).
#' @param cond_vars Names of the conditioned variables.
#' @param n_cond_q Number of conditioned quarters.
#' @param horizon Forecast horizon.
#' @param seed Seed.
#' @return List: seconds_irf (impulse responses, all draws), seconds_total,
#'   seconds_per_draw (total / draws), restrictions (fixed cells), draws,
#'   conditioned_cells_match (largest |forecast - condition| over the fixed cells).
time_conditional <- function(shell, cond_vars, n_cond_q, horizon, seed) {
  n_draws <- shell$meta$n_save
  history <- shell$meta$Y
  colnames(history) <- shell$variables
  base <- forecast_unconditional(shell$beta[1, , , drop = FALSE], shell$sigma[1, , , drop = FALSE],
                                 history, shell$meta$lags, horizon, seed, shocks = FALSE)[1, , ]
  path <- base[seq_len(n_cond_q), cond_vars, drop = FALSE] + 1
  t0 <- proc.time()[["elapsed"]]
  irf <- BVAR::irf(shell, BVAR::bv_irf(horizon = horizon, identification = TRUE, fevd = FALSE))
  t_irf <- proc.time()[["elapsed"]] - t0
  shell$irf <- irf
  t0 <- proc.time()[["elapsed"]]
  # Variables are passed by column position: BVAR 1.0.5 finds names with
  # grep() (partial matching), so "pce" would also match "rpce" and "corepce".
  fc <- bvar_predict_draws(shell, horizon, seed, cond_path = unname(path),
                           cond_vars = match(cond_vars, shell$variables))
  t_fc <- proc.time()[["elapsed"]] - t0
  diffs <- abs(fc[, seq_len(n_cond_q), match(cond_vars, shell$variables), drop = FALSE] -
                 aperm(array(path, c(n_cond_q, length(cond_vars), n_draws)), c(3, 1, 2)))
  list(seconds_irf = t_irf, seconds_conditional = t_fc, seconds_total = t_irf + t_fc,
       seconds_per_draw = (t_irf + t_fc) / n_draws, restrictions = n_cond_q * length(cond_vars),
       draws = n_draws, conditioned_cells_match = max(diffs))
}


# ---- Long-horizon summaries ----

#' Quantiles of forecast draws, in model units or as annualized growth.
#'
#' For log series (model units 100 ln x), the annualized quarter-on-quarter
#' growth is 4 * (y_t - y_{t-1}) in percent (log approximation), computed per
#' draw with the last observed quarter as the base of the first change.
#'
#' @param fc Array draws x H x n.
#' @param history Data matrix (its last row is the base of the first change).
#' @param variable Name of the variable.
#' @param growth TRUE for annualized quarterly growth, FALSE for levels.
#' @param probs Quantiles.
#' @return Matrix H x length(probs).
forecast_quantiles <- function(fc, history, variable, growth, probs = c(0.05, 0.16, 0.5, 0.84, 0.95)) {
  j <- match(variable, dimnames(fc)[[3]])
  x <- fc[, , j]
  if (growth) {
    prev <- cbind(history[nrow(history), j], x[, -ncol(x), drop = FALSE])
    x <- 4 * (x - prev)
  }
  t(apply(x, 2, quantile, probs = probs))
}

#' Fan chart of one variable: recent history and forecast quantiles.
#'
#' @param q Output of forecast_quantiles() (H x 5: 5, 16, 50, 84, 95 percent).
#' @param hist_values Recent history in the same units (vector).
#' @param title Panel title.
plot_fan <- function(q, hist_values, title) {
  H <- nrow(q); k <- length(hist_values)
  x_hist <- seq(-k + 1, 0); x_fc <- seq_len(H)
  ylim <- range(c(q, hist_values), na.rm = TRUE)
  plot(x_hist, hist_values, type = "l", xlim = c(-k + 1, H), ylim = ylim, bty = "n",
       col = ink[["primary"]], lwd = 1.5, xlab = "quarters from the last observation", ylab = "",
       main = title, font.main = 1, cex.main = 0.9, col.main = ink[["primary"]],
       col.axis = ink[["secondary"]])
  # Start the bands and the median at the last observation (x = 0), so the
  # forecast is drawn as a continuation of the history.
  last <- hist_values[k]
  x_fc <- c(0, x_fc)
  q <- rbind(rep(last, ncol(q)), q)
  polygon(c(x_fc, rev(x_fc)), c(q[, 1], rev(q[, 5])), col = "#2a78d633", border = NA)
  polygon(c(x_fc, rev(x_fc)), c(q[, 2], rev(q[, 4])), col = "#2a78d666", border = NA)
  lines(x_fc, q[, 3], col = check_colours[["blue"]], lwd = 2)
  abline(v = 0, col = ink[["grid"]])
}
