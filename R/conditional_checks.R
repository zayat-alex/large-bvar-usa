# conditional_checks.R: checks of the two conditional forecast methods
# (R/forecast_conditional.R).
#
# Why it exists: scenarios are computed by two algorithms that target the
# same conditional distribution: direct Gaussian conditioning (Waggoner and
# Zha 1999) and the Kalman simulation smoother used by Crump et al. (2021),
# Section 2.3. Before either is used, this file checks that they agree, and
# that both agree with a third computation written from the textbook formula
# for a conditional normal distribution, in which the joint covariance of the
# future path is built from powers of the companion matrix (a different
# route from both methods). The comparisons are:
#   1. conditional means, exactly;
#   2. conditional covariances within each quarter, exactly;
#   3. simulated paths with the same seed: both methods are the conditional
#      mean plus A (I - P) z for the same standard normal z, so they should
#      agree path by path up to rounding;
#   4. simulated paths against the exact moments (many paths of one draw);
#   5. production-size draws (one path per posterior draw): quantiles of the
#      two methods, with a second seed as the size of simulation noise.
# The results are written by scripts/check_conditional.R and shown in
# reports/conditional_checks.qmd.
#
# Inputs:  coefficient and covariance draws, a history matrix, conditions
#          matrices (H x n, model units, NA = free).
# Outputs: data.frames of comparisons and timings.


# ---- Textbook reference ----

#' Exact joint mean and covariance of the stacked future path of a VAR.
#'
#' With F the companion matrix (M3) and J = [I_n 0 ... 0], the response of
#' y_{T+h} to the innovation of quarter T+k is Psi_{h-k} = J F^{h-k} J', so
#'   Cov(y_{T+h}, y_{T+g}) = sum_{k=1}^{min(h,g)} Psi_{h-k} Sigma Psi_{g-k}'.
#' The mean is the no-shock path, computed by iterating the companion form.
#'
#' @param B,Sigma One draw (K x n, n x n).
#' @param history Data matrix; its last p rows start the forecast.
#' @param lags,horizon Lags p and horizon H.
#' @return List: mean (nH, stacked quarter by quarter), cov (nH x nH).
forecast_joint_moments_exact <- function(B, Sigma, history, lags, horizon) {
  n_var <- ncol(B)
  F_mat <- companion_matrix(B, lags)
  J <- cbind(diag(n_var), matrix(0, n_var, n_var * (lags - 1)))
  # Mean: z_{t+1} = c_tilde + F z_t on the stacked lags z_t = (y_t', ..., y_{t-p+1}')'.
  z <- as.vector(t(history[nrow(history) - seq_len(lags) + 1, , drop = FALSE]))
  c_tilde <- c(B[1, ], rep(0, n_var * (lags - 1)))
  mean_stacked <- numeric(0)
  for (h in seq_len(horizon)) {
    z <- c_tilde + as.vector(F_mat %*% z)
    mean_stacked <- c(mean_stacked, z[seq_len(n_var)])
  }
  # Psi_j = J F^j J'.
  Psi <- vector("list", horizon)
  F_pow <- diag(n_var * lags)
  for (j in seq_len(horizon)) {
    Psi[[j]] <- J %*% F_pow %*% t(J)
    F_pow <- F_mat %*% F_pow
  }
  cov <- matrix(0, n_var * horizon, n_var * horizon)
  for (h in seq_len(horizon)) for (g in seq_len(horizon)) {
    acc <- matrix(0, n_var, n_var)
    for (k in seq_len(min(h, g))) acc <- acc + Psi[[h - k + 1]] %*% Sigma %*% t(Psi[[g - k + 1]])
    cov[(h - 1) * n_var + seq_len(n_var), (g - 1) * n_var + seq_len(n_var)] <- acc
  }
  list(mean = mean_stacked, cov = cov)
}

#' Textbook conditional normal: mean and covariance given fixed cells.
#'
#' For y ~ N(mu, V) partitioned into fixed cells (1) and free cells (2):
#'   E[y_2 | y_1 = r] = mu_2 + V_21 V_11^{-1} (r - mu_1),
#'   Var[y_2 | y_1]   = V_22 - V_21 V_11^{-1} V_12;
#' fixed cells have their values and zero variance.
#'
#' @param joint Output of forecast_joint_moments_exact().
#' @param conditions H x n matrix, NA = free.
#' @return List: mean (H x n), cov (nH x nH), var (list of H matrices n x n).
condition_textbook <- function(joint, conditions) {
  horizon <- nrow(conditions)
  n_var <- ncol(conditions)
  r_all <- as.vector(t(conditions))
  f1 <- which(!is.na(r_all))
  mean_c <- joint$mean
  cov_c <- joint$cov
  if (length(f1) > 0) {
    V11_inv_V12 <- solve(joint$cov[f1, f1, drop = FALSE], joint$cov[f1, , drop = FALSE])
    mean_c <- joint$mean + as.vector(t(V11_inv_V12) %*% (r_all[f1] - joint$mean[f1]))
    cov_c <- joint$cov - joint$cov[, f1, drop = FALSE] %*% V11_inv_V12
  }
  list(mean = matrix(mean_c, horizon, n_var, byrow = TRUE, dimnames = list(NULL, colnames(conditions))),
       cov = cov_c,
       var = lapply(seq_len(horizon), function(h) {
         idx <- (h - 1) * n_var + seq_len(n_var)
         cov_c[idx, idx]
       }))
}


# ---- Test scenarios ----

#' Test scenarios built around a baseline path.
#'
#' Each scenario fixes some cells at the baseline plus a deviation (model
#' units: log points for log series, percentage points for rates). They are
#' test cases for the methods, not scenarios for publication:
#'   gdp_h1        rgdp + 1 at T+1 (Crump et al. 2021, Section 3.1);
#'   gdp_h8_12     rgdp + 1 at T+8 to T+12 (Section 3.1, the medium-run case);
#'   many_cells_h1 gs2 + 0.25 at T+1 with `hold` at the baseline at T+1 (many
#'                 fixed cells in one quarter, the shape of the recursive
#'                 policy exercise of Section 3.2);
#'   mixed         wti + 10 and gs2, gs10 + 0.5 at T+1 to T+8, sp500 - 5 at
#'                 T+2, T+4, T+6 only (skipped quarters);
#'   heavy         ten variables + 0.5 in all quarters.
#'
#' @param base Baseline path, H x n (for example the median no-shock path).
#' @param hold Variables held at the baseline in many_cells_h1.
#' @return Named list of H x n conditions matrices (NA = free).
make_test_scenarios <- function(base, hold) {
  free <- base
  free[] <- NA_real_
  s <- list()
  s$gdp_h1 <- free
  s$gdp_h1[1, "rgdp"] <- base[1, "rgdp"] + 1
  s$gdp_h8_12 <- free
  s$gdp_h8_12[8:12, "rgdp"] <- base[8:12, "rgdp"] + 1
  s$many_cells_h1 <- free
  s$many_cells_h1[1, hold] <- base[1, hold]
  s$many_cells_h1[1, "gs2"] <- base[1, "gs2"] + 0.25
  s$mixed <- free
  s$mixed[1:8, "wti"] <- base[1:8, "wti"] + 10
  s$mixed[1:8, c("gs2", "gs10")] <- base[1:8, c("gs2", "gs10")] + 0.5
  s$mixed[c(2, 4, 6), "sp500"] <- base[c(2, 4, 6), "sp500"] - 5
  heavy_vars <- c("rgdp", "pce", "corepce", "urate", "gs2", "gs10", "baa", "wti", "sp500", "usd")
  s$heavy <- free
  s$heavy[, heavy_vars] <- base[, heavy_vars] + 0.5
  s
}


# ---- Comparisons 1-2: exact moments ----

#' Conditional means and per-quarter variances of the two methods, draw by
#' draw, and against the textbook formula on some draws.
#'
#' Variance differences are scaled by the unconditional variance of the same
#' cell, so they are comparable across variables and horizons. The table
#' also records two properties of a conditional distribution: the variance
#' of a fixed cell is zero, and no conditional variance exceeds the
#' unconditional one.
#'
#' @param beta,sigma Draws.
#' @param history,lags,horizon As in forecast_conditional().
#' @param scenarios Named list of conditions matrices.
#' @param draw_index Draws to compare.
#' @param textbook_index Draws (a subset) also compared with the textbook.
#' @return data.frame, one row per draw and scenario.
compare_moments <- function(beta, sigma, history, lags, horizon, scenarios, draw_index, textbook_index) {
  free <- scenarios[[1]]
  free[] <- NA_real_
  rows <- list()
  for (i in draw_index) {
    B <- beta[i, , ]
    S <- sigma[i, , ]
    var_unc <- sapply(conditional_moments_gaussian(B, S, history, lags, horizon, free)$var, diag)  # n x H
    joint <- if (i %in% textbook_index) forecast_joint_moments_exact(B, S, history, lags, horizon) else NULL
    for (nm in names(scenarios)) {
      cond <- scenarios[[nm]]
      fixed <- which(!is.na(t(cond)))                 # positions in n x H, as var_unc
      g <- conditional_moments_gaussian(B, S, history, lags, horizon, cond)
      k <- conditional_moments_kalman(B, S, history, lags, horizon, cond)
      vg <- sapply(g$var, diag)
      vk <- sapply(k$var, diag)
      row <- data.frame(
        scenario = nm, draw = i, fixed_cells = length(fixed),
        mean_max_abs_diff = max(abs(g$mean - k$mean)),
        var_max_rel_diff = max(abs(vg - vk) / var_unc),
        fixed_var_max_rel = max(c(vg[fixed], vk[fixed]) / var_unc[fixed]),
        var_above_unc_max_rel = max((pmax(vg, vk) - var_unc) / var_unc),
        textbook_mean_max_abs_diff = NA_real_, textbook_var_max_rel_diff = NA_real_)
      if (!is.null(joint)) {
        tb <- condition_textbook(joint, cond)
        vt <- sapply(tb$var, diag)
        row$textbook_mean_max_abs_diff <- max(abs(c(g$mean - tb$mean, k$mean - tb$mean)))
        row$textbook_var_max_rel_diff <- max(abs(cbind(vg - vt, vk - vt)) / cbind(var_unc, var_unc))
      }
      rows[[length(rows) + 1]] <- row
    }
  }
  do.call(rbind, rows)
}


# ---- Comparisons 3 and 5: paths ----

#' Paths of the two methods with the same seed, and quantiles at production size.
#'
#' Runs the Kalman and Gaussian methods with the same seed (comparison 3:
#' same paths up to rounding), and the Gaussian method with a second seed,
#' whose differences from the first show the size of simulation noise with
#' one path per draw (comparison 5).
#'
#' @param beta,sigma Draws (all draws used).
#' @param history,lags,horizon,conditions As in forecast_conditional().
#' @param seed,seed_alt The common seed and the second seed.
#' @param variables,horizons Cells of the quantile table.
#' @return List: paths (kalman, gaussian, gaussian_alt), max_abs_diff,
#'   max_abs_level, fixed_cells_max_error, ms_per_draw, quantiles (data.frame
#'   in change from the last observation: q05, q50, q95 by method; the
#'   difference of medians from the Gaussian run, as a share of its 90
#'   percent band width).
compare_paths <- function(beta, sigma, history, lags, horizon, conditions, seed, seed_alt,
                          variables, horizons) {
  n_draws <- dim(beta)[1]
  timed <- function(method, s) {
    t0 <- proc.time()[["elapsed"]]
    x <- forecast_conditional(beta, sigma, history, lags, horizon, conditions, s, method = method)
    list(x = x, ms = 1000 * (proc.time()[["elapsed"]] - t0) / n_draws)
  }
  k  <- timed("kalman", seed)
  g  <- timed("gaussian", seed)
  g2 <- timed("gaussian", seed_alt)
  fixed <- which(!is.na(conditions))
  fixed_err <- 0
  if (length(fixed)) {
    for (p in list(k$x, g$x)) {
      for (i in seq_len(n_draws)) fixed_err <- max(fixed_err, abs(p[i, , ][fixed] - conditions[fixed]))
    }
  }
  last <- history[nrow(history), ]
  rows <- list()
  for (v in variables) for (h in horizons) {
    j <- match(v, colnames(history))
    q <- lapply(list(kalman = k$x, gaussian = g$x, gaussian_alt = g2$x),
                function(p) quantile(p[, h, j] - last[j], c(0.05, 0.5, 0.95), names = FALSE))
    width <- q$gaussian[3] - q$gaussian[1]
    rows[[length(rows) + 1]] <- data.frame(
      variable = v, horizon = h, fixed = !is.na(conditions[h, j]),
      q05_kalman = q$kalman[1], q50_kalman = q$kalman[2], q95_kalman = q$kalman[3],
      q05_gaussian = q$gaussian[1], q50_gaussian = q$gaussian[2], q95_gaussian = q$gaussian[3],
      median_diff_kalman_share = if (width > 0) (q$kalman[2] - q$gaussian[2]) / width else NA_real_,
      median_diff_seed_share = if (width > 0) (q$gaussian_alt[2] - q$gaussian[2]) / width else NA_real_)
  }
  list(paths = list(kalman = k$x, gaussian = g$x, gaussian_alt = g2$x),
       max_abs_diff = max(abs(k$x - g$x)), max_abs_level = max(abs(g$x)),
       fixed_cells_max_error = fixed_err,
       ms_per_draw = c(kalman = k$ms, gaussian = g$ms),
       quantiles = do.call(rbind, rows))
}


# ---- Comparison 4: simulated paths against exact moments ----

#' Many paths of one draw from each method against the textbook moments.
#'
#' For every free cell: z = (simulated mean - exact mean) / (exact s.d. /
#' sqrt(paths)), and the ratio of simulated to exact variance. With N free
#' cells and independent z, the median of the largest |z| is about
#' qnorm(1 - 0.25 / N); it is reported beside the result. The cells of one
#' path are strongly correlated (a variable at neighbouring horizons), so
#' these per-cell summaries are rough guides. The joint test accounts for
#' the correlation: with mbar the vector of simulated means of the free
#' cells, mu and V their exact conditional mean and covariance, and N paths,
#'   N (mbar - mu)' V^{-1} (mbar - mu) ~ chi-squared(number of free cells)
#' if the simulated paths have the right distribution (approximately, by the
#' central limit theorem; exactly for Gaussian paths).
#'
#' @param B,Sigma One draw.
#' @param history,lags,horizon,conditions As in forecast_conditional().
#' @param n_paths Named integer vector: paths per method (kalman, gaussian).
#' @param seed Seed.
#' @return data.frame, one row per method.
check_simulated_vs_exact <- function(B, Sigma, history, lags, horizon, conditions, n_paths, seed) {
  tb <- condition_textbook(forecast_joint_moments_exact(B, Sigma, history, lags, horizon), conditions)
  rows <- list()
  for (m in names(n_paths)) {
    N <- n_paths[[m]]
    t0 <- proc.time()[["elapsed"]]
    p <- forecast_conditional(array(B, c(1, dim(B))), array(Sigma, c(1, dim(Sigma))), history, lags,
                              horizon, conditions, seed, method = m, n_paths = N)
    secs <- proc.time()[["elapsed"]] - t0
    rows[[length(rows) + 1]] <- cbind(method = m, paths_vs_exact_stats(p, tb, conditions),
                                      ms_per_path = 1000 * secs / N)
  }
  do.call(rbind, rows)
}

#' Statistics of simulated paths of one draw against its exact conditional moments.
#'
#' The statistics of check_simulated_vs_exact(), for paths from any source
#' (also BVAR::predict(), R/engine_checks.R).
#'
#' @param p Paths, array N x H x n, all from the same draw and scenario.
#' @param tb Exact conditional moments, output of condition_textbook().
#' @param conditions The scenario's conditions matrix (H x n, NA = free).
#' @return One-row data.frame: paths, free_cells, chi2, chi2_p_value,
#'   max_abs_z, typical_max_abs_z, share_abs_z_above_2, var_ratio_min,
#'   var_ratio_max, var_ratio_max_abs_dev_in_se.
paths_vs_exact_stats <- function(p, tb, conditions) {
  N <- dim(p)[1]
  horizon <- nrow(conditions)
  free <- which(is.na(conditions))                      # positions in the H x n matrix
  # The same cells in the stacked (quarter by quarter) order of tb$cov.
  free_stacked <- which(is.na(as.vector(t(conditions))))
  V_chol <- chol(tb$cov[free_stacked, free_stacked])
  sd_exact <- sqrt(matrix(diag(tb$cov), horizon, ncol(conditions), byrow = TRUE)[free])
  sim_mean <- apply(p, c(2, 3), mean)
  z <- (sim_mean[free] - tb$mean[free]) / (sd_exact / sqrt(N))
  ratio <- apply(p, c(2, 3), var)[free] / sd_exact^2
  # Joint test: e = mbar - mu in stacked order; chi2 = N e' V^{-1} e.
  e <- (as.vector(t(sim_mean)) - as.vector(t(tb$mean)))[free_stacked]
  chi2 <- N * sum(backsolve(V_chol, e, transpose = TRUE)^2)
  data.frame(paths = N, free_cells = length(free),
             chi2 = chi2, chi2_p_value = pchisq(chi2, df = length(free), lower.tail = FALSE),
             max_abs_z = max(abs(z)), typical_max_abs_z = qnorm(1 - 0.25 / length(free)),
             share_abs_z_above_2 = mean(abs(z) > 2),
             var_ratio_min = min(ratio), var_ratio_max = max(ratio),
             var_ratio_max_abs_dev_in_se = max(abs(ratio - 1)) / sqrt(2 / (N - 1)))
}


# ---- Figure ----

#' Fan charts of one scenario under both methods, against the baseline.
#'
#' One panel per variable, in change from the last observation (model
#' units). Grey: the baseline's 90 percent band. Blue: the scenario's 68 and
#' 90 percent bands and median from the Gaussian method. Orange dashed: the
#' Kalman method's median and 90 percent band edges (same seed), which should
#' lie on the blue ones. Dots: the fixed cells.
#'
#' @param baseline Unconditional forecast draws (draws x H x n).
#' @param kalman,gaussian Conditional forecast draws of the two methods.
#' @param conditions The scenario's conditions matrix.
#' @param history Data matrix (its last row is the base of the changes).
#' @param variables Variables to show.
#' @param units Named character vector of unit labels per variable.
plot_conditional_fans <- function(baseline, kalman, gaussian, conditions, history, variables, units) {
  H <- dim(gaussian)[2]
  x <- seq_len(H)
  last <- history[nrow(history), ]
  qs <- function(p, j) t(apply(p[, , j] - last[j], 2, quantile, c(0.05, 0.16, 0.5, 0.84, 0.95)))
  for (v in variables) {
    j <- match(v, colnames(history))
    qb <- qs(baseline, j)
    qg <- qs(gaussian, j)
    qk <- qs(kalman, j)
    fx <- which(!is.na(conditions[, j]))
    ylim <- range(c(qb, qg, qk, conditions[fx, j] - last[j], 0))
    plot(x, qg[, 3], type = "n", ylim = ylim, bty = "n", xlab = "quarters ahead", ylab = units[[v]],
         main = v, font.main = 1, cex.main = 0.95, col.main = ink[["primary"]],
         col.axis = ink[["secondary"]], col.lab = ink[["secondary"]])
    abline(h = 0, col = ink[["grid"]])
    polygon(c(x, rev(x)), c(qb[, 1], rev(qb[, 5])), col = "#52514e22", border = NA)
    polygon(c(x, rev(x)), c(qg[, 1], rev(qg[, 5])), col = "#2a78d633", border = NA)
    polygon(c(x, rev(x)), c(qg[, 2], rev(qg[, 4])), col = "#2a78d666", border = NA)
    lines(x, qg[, 3], col = check_colours[["blue"]], lwd = 2)
    for (col_q in c(1, 3, 5)) lines(x, qk[, col_q], col = check_colours[["orange"]], lwd = 1.5, lty = 2)
    if (length(fx)) points(fx, conditions[fx, j] - last[j], pch = 19, cex = 0.9, col = ink[["primary"]])
  }
}

#' Legend for plot_conditional_fans(), drawn across the bottom of the figure.
#'
#' Call after the panels, with the figure's outer bottom margin (par(oma))
#' at least 2 lines.
legend_conditional_fans <- function() {
  op <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  on.exit(par(op))
  plot.new()
  legend("bottom", horiz = TRUE, bty = "n", cex = 0.9, text.col = ink[["secondary"]],
         legend = c("baseline 90%", "scenario 90% / 68% (Gaussian)", "scenario median (Gaussian)",
                    "Kalman median and 90% edges", "fixed cells"),
         fill = c("#52514e22", "#2a78d655", NA, NA, NA), border = NA,
         col = c(NA, NA, check_colours[["blue"]], check_colours[["orange"]], ink[["primary"]]),
         lty = c(NA, NA, 1, 2, NA), lwd = c(NA, NA, 2, 1.5, NA), pch = c(NA, NA, NA, NA, 19))
}
