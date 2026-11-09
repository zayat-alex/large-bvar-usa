# engine_checks.R: properties of the scenario engine and the cross-check
# against the BVAR package's conditional forecasts.
#
# Why it exists: a conditional forecast has properties that hold exactly for
# every posterior draw, whatever the data (a list of them is the basis of
# the engine's tests). They are checked here on the model's own draws, for
# both conditional forecast methods (R/forecast_conditional.R):
#   null scenario   fixing cells at the draw's own no-shock path returns the
#                   no-shock path everywhere (a scenario that changes nothing
#                   has no effect);
#   linearity       scaling the deviation of the fixed cells from the
#                   no-shock path by alpha scales the effect on every
#                   conditional mean by alpha;
#   superposition   the effect of two deviations together equals the sum of
#                   their effects, provided each is computed with the other's
#                   cells fixed at the no-shock path (the conditioning set
#                   must be the same; fixing a cell, even at its baseline
#                   value, changes what the model infers);
#   ordering        the results do not depend on the order of the variables
#                   in the model: with all n shocks free, the conditional
#                   distribution is the same whatever Cholesky ordering is
#                   used internally.
# The fixed cells being met, zero variance at fixed cells and conditional
# variances no larger than unconditional ones are checked in
# R/conditional_checks.R.
#
# Cross-check: BVAR::predict() with cond_path implements the Waggoner and
# Zha (1999) method independently (BVAR 1.0.5, read in source: cond_fcast()
# draws structural shocks eta = R+ r + V2 z from the Cholesky impulse
# responses). Its paths for one draw are compared with the exact conditional
# moments, and its quantiles at production size with the project's. Its
# handling of a single fixed cell is examined separately (it builds the
# inverse of the singular values with diag(1 / d), which for a scalar d is an
# identity matrix of size floor(1 / d), not the 1 x 1 matrix 1 / d).
#
# Used by the tests (tests/testthat/test-forecast-conditional.R, test-irf.R).
#
# Inputs:  coefficient and covariance draws, a history matrix, conditions
#          matrices (H x n, model units, NA = free).
# Outputs: data.frames of property checks and comparisons.


# ---- Helpers ----

#' Conditional mean and per-quarter covariances by either method.
#'
#' @param method "gaussian" or "kalman".
#' @return List: mean (H x n), var (list of H matrices n x n).
conditional_moments <- function(B, Sigma, history, lags, horizon, conditions, method) {
  switch(method,
         gaussian = conditional_moments_gaussian(B, Sigma, history, lags, horizon, conditions),
         kalman = conditional_moments_kalman(B, Sigma, history, lags, horizon, conditions),
         stop("unknown method: ", method))
}

#' Fix the cells of a template at a path plus a deviation.
#'
#' @param template Conditions matrix whose non-NA cells are the cells to fix.
#' @param path H x n path (for example a draw's no-shock path).
#' @param deviation H x n deviations added at the fixed cells (0 for none).
#' @return Conditions matrix: path + deviation at the template's cells, NA elsewhere.
fix_cells_at <- function(template, path, deviation = 0) {
  out <- template
  fixed <- !is.na(template)
  out[fixed] <- (path + deviation)[fixed]
  out
}

#' Reorder the variables of a VAR.
#'
#' Row 1 of B (the constant) stays; within each lag block the rows follow
#' `perm`; the columns (equations) follow `perm`.
#'
#' @param B,Sigma One draw (K x n, n x n).
#' @param history,conditions Matrices with n columns.
#' @param perm Permutation of 1..n.
#' @param lags Number of lags p.
#' @return List: B, Sigma, history, conditions, all in the new order.
permute_variables <- function(B, Sigma, history, conditions, perm, lags) {
  n_var <- ncol(B)
  rows <- c(1, unlist(lapply(seq_len(lags), function(s) 1 + (s - 1) * n_var + perm)))
  list(B = B[rows, perm], Sigma = Sigma[perm, perm], history = history[, perm, drop = FALSE],
       conditions = conditions[, perm, drop = FALSE])
}


# ---- Properties ----

#' The engine's exact properties for one draw, one scenario and one method.
#'
#' The scenario is used as a template: its fixed cells and its deviation
#' `delta` from the baseline it was built around. For the draw, the
#' deviation is applied to the draw's own no-shock path m, so the property
#' is exact for the draw. Differences are reported as the largest absolute
#' difference and as a share of the largest effect involved.
#'
#' @param B,Sigma One draw.
#' @param history,lags,horizon As in forecast_conditional().
#' @param template Conditions matrix of the scenario.
#' @param delta H x n deviations of the scenario from its baseline (0 where free).
#' @param partner Second template for superposition (with its delta), or NULL.
#' @param partner_delta Deviations of the partner.
#' @param alphas Scale factors for linearity.
#' @param perms Named list of permutations for the ordering check.
#' @param method "gaussian" or "kalman".
#' @return data.frame: property, detail, max_abs_diff, scale, rel_diff.
engine_properties <- function(B, Sigma, history, lags, horizon, template, delta, partner, partner_delta,
                              alphas, perms, method) {
  m <- forecast_unconditional(array(B, c(1, dim(B))), array(Sigma, c(1, dim(Sigma))), history, lags,
                              horizon, seed = 1, shocks = FALSE)[1, , ]
  effect <- function(cond) conditional_moments(B, Sigma, history, lags, horizon, cond, method)$mean - m
  row <- function(property, detail, diff, scale) {
    data.frame(property = property, detail = detail, max_abs_diff = max(abs(diff)),
               scale = scale, rel_diff = max(abs(diff)) / max(scale, .Machine$double.eps))
  }
  rows <- list()

  # Null scenario: cells fixed at the draw's no-shock path.
  rows[[1]] <- row("null scenario", "cells at the no-shock path", effect(fix_cells_at(template, m)),
                   max(abs(m)))

  # Linearity: effect of alpha * delta = alpha * effect of delta.
  e1 <- effect(fix_cells_at(template, m, delta))
  for (a in alphas) {
    ea <- effect(fix_cells_at(template, m, a * delta))
    rows[[length(rows) + 1]] <- row("linearity", paste0("alpha = ", a), ea - a * e1, max(abs(a * e1)))
  }

  # Superposition, on the union of the two sets of fixed cells.
  if (!is.null(partner)) {
    union <- template
    union[!is.na(partner)] <- 0
    union[!is.na(template)] <- 0
    d1 <- ifelse(is.na(template), 0, delta)
    d2 <- ifelse(is.na(partner), 0, partner_delta)
    e_a  <- effect(fix_cells_at(union, m, d1))
    e_b  <- effect(fix_cells_at(union, m, d2))
    e_ab <- effect(fix_cells_at(union, m, d1 + d2))
    rows[[length(rows) + 1]] <- row("superposition", "effect(a + b) - effect(a) - effect(b)",
                                    e_ab - e_a - e_b, max(abs(e_ab)))
  }

  # Ordering: reorder the variables, compute, put back in the original order.
  cond <- fix_cells_at(template, m, delta)
  base_m <- conditional_moments(B, Sigma, history, lags, horizon, cond, method)
  for (nm in names(perms)) {
    perm <- perms[[nm]]
    pv <- permute_variables(B, Sigma, history, cond, perm, lags)
    pm <- conditional_moments(pv$B, pv$Sigma, pv$history, lags, horizon, pv$conditions, method)
    back <- order(perm)
    rows[[length(rows) + 1]] <- row("ordering", paste0(nm, ": mean"), pm$mean[, back] - base_m$mean,
                                    max(abs(base_m$mean)))
    var_diff <- unlist(lapply(seq_len(horizon), function(h) pm$var[[h]][back, back] - base_m$var[[h]]))
    rows[[length(rows) + 1]] <- row("ordering", paste0(nm, ": variance"), var_diff,
                                    max(unlist(lapply(base_m$var, diag))))
  }
  do.call(rbind, rows)
}


# ---- BVAR cross-check ----

#' BVAR's Cholesky impulse responses for one draw.
#'
#' BVAR::irf() requires at least 10 draws (its n_thin check), so the draw is
#' repeated 10 times and the first copy kept.
#'
#' @param B,Sigma One draw.
#' @param history,lags,horizon As in forecast_conditional().
#' @return The bvar_irf object for 10 copies of the draw; element irf is
#'   10 x n (response) x H (horizon, 1 = impact) x n (shock).
bvar_irf_one_draw <- function(B, Sigma, history, lags, horizon) {
  shell <- bvar_shell(array(rep(B, each = 10), c(10, dim(B))), array(rep(Sigma, each = 10), c(10, dim(Sigma))),
                      history, lags)
  BVAR::irf(shell, BVAR::bv_irf(horizon = horizon, identification = TRUE, fevd = FALSE))
}

#' BVAR's Cholesky impulse responses against the project's Psi_{h-1} L.
#'
#' @return Largest absolute difference over responses, horizons and shocks.
bvar_irf_agreement <- function(B, Sigma, history, lags, horizon) {
  ir <- bvar_irf_one_draw(B, Sigma, history, lags, horizon)
  Psi <- ma_coefficients(B, lags, horizon)
  L <- t(chol(Sigma))
  own <- array(NA_real_, c(ncol(B), horizon, ncol(B)))
  for (h in seq_len(horizon)) own[, h, ] <- Psi[[h]] %*% L
  max(abs(ir$irf[1, , , ] - own))
}

#' BVAR::predict() conditional paths for one draw, many paths.
#'
#' The draw is repeated n_paths times; its impulse responses are computed
#' once and repeated, so BVAR does not recompute them for every copy.
#' Conditioned variables are passed by column position (BVAR 1.0.5 finds
#' names with grep(), so "pce" would also match "rpce" and "corepce").
#'
#' @return Array n_paths x H x n.
bvar_conditional_one_draw <- function(B, Sigma, history, lags, horizon, conditions, n_paths, seed) {
  ir <- bvar_irf_one_draw(B, Sigma, history, lags, horizon)
  ir$irf <- ir$irf[rep(1, n_paths), , , , drop = FALSE]
  shell <- bvar_shell(array(rep(B, each = n_paths), c(n_paths, dim(B))),
                      array(rep(Sigma, each = n_paths), c(n_paths, dim(Sigma))), history, lags)
  shell$irf <- ir
  vars <- which(colSums(!is.na(conditions)) > 0)
  p <- bvar_predict_draws(shell, horizon, seed, cond_path = unname(conditions[, vars, drop = FALSE]),
                          cond_vars = vars)
  dimnames(p) <- list(NULL, paste0("h", seq_len(horizon)), colnames(history))
  p
}

#' BVAR with a single fixed cell.
#'
#' For one fixed cell (variable j, quarter h) the restriction matrix R has
#' one row, whose norm d is the unconditional standard deviation of that
#' cell. BVAR 1.0.5 builds diag(1 / d), which is an identity matrix of size
#' floor(1 / d): 1 x 1 (= 1, not 1 / d) when 1 <= 1 / d < 2, so the fixed
#' cell is moved by d r instead of r; a 0 x 0 or larger matrix otherwise,
#' which makes the matrix product fail. This function records what happens.
#'
#' @param B,Sigma One draw.
#' @param history,lags,horizon As in forecast_conditional().
#' @param cells data.frame: variable, quarter, deviation (from the draw's
#'   no-shock path).
#' @return data.frame per cell: d, 1/d, requested deviation r, BVAR's
#'   deviation at the fixed cell (NA if it failed), d r, error message.
bvar_single_cell <- function(B, Sigma, history, lags, horizon, cells) {
  m <- forecast_unconditional(array(B, c(1, dim(B))), array(Sigma, c(1, dim(Sigma))), history, lags,
                              horizon, seed = 1, shocks = FALSE)[1, , ]
  sd_unc <- sqrt(matrix(diag(forecast_joint_moments_exact(B, Sigma, history, lags, horizon)$cov),
                        horizon, ncol(B), byrow = TRUE))
  rows <- list()
  for (k in seq_len(nrow(cells))) {
    v <- cells$variable[k]; h <- cells$quarter[k]; r <- cells$deviation[k]
    j <- match(v, colnames(history))
    cond <- matrix(NA_real_, horizon, ncol(B), dimnames = list(NULL, colnames(history)))
    cond[h, j] <- m[h, j] + r
    res <- tryCatch(bvar_conditional_one_draw(B, Sigma, history, lags, horizon, cond, n_paths = 10, seed = 1),
                    error = function(e) conditionMessage(e))
    d <- sd_unc[h, j]
    rows[[k]] <- data.frame(variable = v, quarter = h, d = d, one_over_d = 1 / d, requested = r,
                            bvar_deviation = if (is.character(res)) NA_real_ else mean(res[, h, j]) - m[h, j],
                            d_times_r = d * r,
                            error = if (is.character(res)) res else "")
  }
  do.call(rbind, rows)
}

#' BVAR and the project at production size: quantiles of one scenario.
#'
#' All draws, one path each: BVAR::predict() (impulse responses computed for
#' every draw, as in production use), the project's Gaussian method, and the
#' Gaussian method with a second seed (the size of simulation noise).
#'
#' @return List: quantiles (data.frame in change from the last observation,
#'   with the difference of medians as a share of the 90 percent band width,
#'   BVAR vs project and seed vs seed), fixed_cells_max_error (BVAR), seconds.
bvar_production_compare <- function(beta, sigma, history, lags, horizon, conditions, seed, seed_alt,
                                    variables, horizons) {
  t0 <- proc.time()[["elapsed"]]
  shell <- bvar_shell(beta, sigma, history, lags)
  shell$irf <- BVAR::irf(shell, BVAR::bv_irf(horizon = horizon, identification = TRUE, fevd = FALSE))
  vars <- which(colSums(!is.na(conditions)) > 0)
  pb <- bvar_predict_draws(shell, horizon, seed, cond_path = unname(conditions[, vars, drop = FALSE]),
                           cond_vars = vars)
  t_bvar <- proc.time()[["elapsed"]] - t0
  pg  <- forecast_conditional(beta, sigma, history, lags, horizon, conditions, seed, method = "gaussian")
  pg2 <- forecast_conditional(beta, sigma, history, lags, horizon, conditions, seed_alt, method = "gaussian")
  fixed <- which(!is.na(conditions))
  fixed_err <- max(apply(pb, 1, function(path) max(abs(path[fixed] - conditions[fixed]))))
  last <- history[nrow(history), ]
  rows <- list()
  for (v in variables) for (h in horizons) {
    j <- match(v, colnames(history))
    q <- lapply(list(bvar = pb, gaussian = pg, gaussian_alt = pg2),
                function(p) quantile(p[, h, j] - last[j], c(0.05, 0.5, 0.95), names = FALSE))
    width <- q$gaussian[3] - q$gaussian[1]
    rows[[length(rows) + 1]] <- data.frame(
      variable = v, horizon = h, fixed = !is.na(conditions[h, j]),
      q05_bvar = q$bvar[1], q50_bvar = q$bvar[2], q95_bvar = q$bvar[3],
      q05_gaussian = q$gaussian[1], q50_gaussian = q$gaussian[2], q95_gaussian = q$gaussian[3],
      median_diff_bvar_share = if (width > 0) (q$bvar[2] - q$gaussian[2]) / width else NA_real_,
      median_diff_seed_share = if (width > 0) (q$gaussian_alt[2] - q$gaussian[2]) / width else NA_real_,
      width_ratio_bvar = if (width > 0) (q$bvar[3] - q$bvar[1]) / width else NA_real_,
      width_ratio_seed = if (width > 0) (q$gaussian_alt[3] - q$gaussian_alt[1]) / width else NA_real_)
  }
  list(quantiles = do.call(rbind, rows), fixed_cells_max_error = fixed_err, seconds_bvar = t_bvar)
}
