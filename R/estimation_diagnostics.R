# estimation_diagnostics.R: residual diagnostics of the estimated VAR.
#
# Why it exists: the MCMC report of each run (reports/mcmc_diagnostics.qmd)
# checks the sampler, the hyperparameters against their priors and bounds,
# and the stability of the coefficient draws. What it does not check is the
# fit of the VAR itself: are the one-step-ahead residuals free of large
# outliers and of remaining autocorrelation, and do they co-move as one
# would expect (CPI with core CPI, for example)? A VAR whose residuals are
# strongly autocorrelated has lags missing; isolated huge residuals show
# quarters the model treats as extreme shocks, which then shape Sigma.
#
# Notation: y_t (n variables, model units), x_t = (1, y_{t-1}', ..., y_{t-p}')'
# the regressors, B (K x n) the coefficients, Sigma the residual covariance.
# Residuals e_t = y_t - B' x_t; standardized z_{jt} = e_{jt} / sqrt(Sigma_jj).
# B and Sigma are the posterior means of the run's draws: a point estimate
# for a description of the fit, not a posterior statement. The mean, not the
# element-wise median: residuals are linear in B, so the residual at the mean
# B is the average of the residuals over the draws, whereas element-wise
# medians mix coefficients from different draws. With variables in levels
# (100 x log, around 10^3) and roots near 1, that mixing breaks the
# combinations of coefficients that hold within each draw and leaves large
# constant offsets in the residuals (found when this check was written).
#
# Inputs:  a run's estimation data (covid.rds) and draws (draws.rds).
# Outputs: data.frames and figures used by reports/estimation_diagnostics.qmd.


#' Residuals of the VAR at the posterior mean of (B, Sigma).
#'
#' @param y Estimation data (quarters x n, model units), rows named by quarter.
#' @param draws The run's draws.rds (beta: draws x K x n, sigma: draws x n x n, lags).
#' @return List: e (residuals, quarters x n), z (standardized residuals),
#'   B, Sigma (posterior means), quarters (labels of the residual rows).
var_residuals <- function(y, draws) {
  lags <- draws$lags
  B     <- apply(draws$beta, c(2, 3), mean)
  Sigma <- apply(draws$sigma, c(2, 3), mean)
  dm <- var_data_matrices(y, lags)
  e <- dm$Y - dm$X %*% B
  quarters <- rownames(y)[(lags + 1):nrow(y)]
  dimnames(e) <- list(quarters, colnames(y))
  z <- sweep(e, 2, sqrt(diag(Sigma)), "/")
  list(e = e, z = z, B = B, Sigma = Sigma, quarters = quarters)
}

#' Standardized residuals beyond a threshold.
#'
#' @param z Standardized residuals (quarters x n).
#' @param threshold Absolute value flagged (spec 9.3: 5).
#' @return data.frame: quarter, variable, z, sorted by |z|.
residual_outliers <- function(z, threshold = 5) {
  idx <- which(abs(z) > threshold, arr.ind = TRUE)
  if (!nrow(idx)) return(data.frame(quarter = character(0), variable = character(0), z = numeric(0)))
  out <- data.frame(quarter = rownames(z)[idx[, 1]], variable = colnames(z)[idx[, 2]],
                    z = z[idx])
  out[order(-abs(out$z)), , drop = FALSE]
}

#' Autocorrelation of each variable's residuals.
#'
#' A well-specified VAR leaves residuals close to white noise. Reported per
#' variable: the autocorrelations at lags 1 to max_lag, the largest in
#' absolute value, and the Ljung-Box test of no autocorrelation up to
#' max_lag (stats::Box.test). With 31 variables, a few p-values below 0.05
#' are expected by chance; the test also assumes estimated residuals behave
#' like data, which they only approximately do.
#'
#' @param e Residuals (quarters x n).
#' @param max_lag Largest lag (spec: 8).
#' @return data.frame: variable, acf_1 .. acf_<max_lag>, max_abs_acf, ljung_box_p.
residual_autocorrelation <- function(e, max_lag = 8) {
  rows <- lapply(colnames(e), function(v) {
    a <- stats::acf(e[, v], lag.max = max_lag, plot = FALSE)$acf[-1]
    data.frame(variable = v, t(setNames(a, paste0("acf_", seq_len(max_lag)))),
               max_abs_acf = max(abs(a)),
               ljung_box_p = stats::Box.test(e[, v], lag = max_lag, type = "Ljung-Box")$p.value)
  })
  do.call(rbind, rows)
}

#' Largest residual correlations between pairs of variables.
#'
#' @param e Residuals (quarters x n).
#' @param k Number of pairs.
#' @return data.frame: variable_1, variable_2, correlation, the k largest |correlation|.
top_residual_correlations <- function(e, k = 15) {
  r <- stats::cor(e)
  idx <- which(upper.tri(r), arr.ind = TRUE)
  out <- data.frame(variable_1 = colnames(r)[idx[, 1]], variable_2 = colnames(r)[idx[, 2]],
                    correlation = r[idx])
  head(out[order(-abs(out$correlation)), , drop = FALSE], k)
}

#' Rows of a residual matrix before a quarter (e.g. to leave out 2020-21).
#'
#' @param m Matrix with rows named by quarter ("2019Q4").
#' @param before First quarter left out.
#' @return The rows before it.
rows_before <- function(m, before) {
  m[parse_quarter(rownames(m)) < parse_quarter(before), , drop = FALSE]
}


# ---- Plots (base graphics; colours of R/estimate_checks.R) ----

#' Standardized residuals over time, one row per variable.
#'
#' Cells beyond +/-3 are coloured (orange positive, blue negative; darker
#' beyond +/-5), so outlier quarters show as vertical streaks.
#'
#' @param z Standardized residuals (quarters x n).
#' @param title Plot title.
plot_residual_grid <- function(z, title) {
  cls <- cut(z, c(-Inf, -5, -3, 3, 5, Inf), labels = FALSE)
  cols <- c("#2a78d6", "#2a78d666", "#ffffff", "#eb683466", "#eb6834")
  m <- matrix(cls, nrow(z))
  n <- ncol(z)
  op <- par(mar = c(4, 8, 2.5, 1), col.axis = ink[["secondary"]])
  on.exit(par(op))
  image(seq_len(nrow(z)), seq_len(n), m[, n:1, drop = FALSE], col = cols, breaks = 0.5:5.5,
        xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = title, font.main = 1,
        col.main = ink[["primary"]])
  years <- substr(rownames(z), 1, 4)
  at <- which(grepl("Q1$", rownames(z)) & as.integer(years) %% 5 == 0)
  axis(1, at = at, labels = years[at], cex.axis = 0.7)
  axis(2, at = seq_len(n), labels = rev(colnames(z)), las = 1, cex.axis = 0.6, tick = FALSE)
  box(col = ink[["grid"]])
  legend("bottom", inset = c(0, -0.13), xpd = TRUE, horiz = TRUE, bty = "n", cex = 0.7,
         fill = cols[-3], legend = c("< -5", "-5 to -3", "3 to 5", "> 5"), border = NA)
}

#' Residual correlation matrix as a heat map.
#'
#' @param e Residuals (quarters x n).
#' @param title Plot title.
plot_residual_correlation <- function(e, title) {
  r <- stats::cor(e)
  n <- ncol(r)
  pal <- grDevices::colorRampPalette(c("#2a78d6", "#ffffff", "#eb6834"))(41)
  op <- par(mar = c(6, 6, 2.5, 1), col.axis = ink[["secondary"]])
  on.exit(par(op))
  image(seq_len(n), seq_len(n), r[, n:1], col = pal, zlim = c(-1, 1), xaxt = "n", yaxt = "n",
        xlab = "", ylab = "", main = title, font.main = 1, col.main = ink[["primary"]])
  axis(1, at = seq_len(n), labels = colnames(r), las = 2, cex.axis = 0.55, tick = FALSE)
  axis(2, at = seq_len(n), labels = rev(colnames(r)), las = 1, cex.axis = 0.55, tick = FALSE)
}
