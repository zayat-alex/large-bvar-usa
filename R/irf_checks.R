# irf_checks.R: checks of the impulse responses (R/irf.R).
#
# Why it exists: the impulse responses are computed with the closed-form
# Cholesky formula, but Crump et al. (2021), Section 2.3, define them as a
# conditional forecast minus the baseline. The two must agree exactly for
# every posterior draw, and that is checked here:
#   equivalence   the generalized IRF of variable j equals the effect on the
#                 conditional mean of setting y_{j,T+1} to its baseline + 1
#                 (everything else free); the recursive IRF equals the
#                 effect of setting the slow variables to their baseline and
#                 the policy variable to baseline + 1 in T + 1 (fast
#                 variables free). The conditional means come from the
#                 scenario engine (R/forecast_conditional.R), which
#                 shares no code with R/irf.R except the moving-average
#                 coefficients.
#   invariance    the recursive IRF does not depend on the order of the
#                 variables within the slow block or within the fast block.
#   BVAR          BVAR::irf() (BVAR 1.0.5) computes Cholesky responses to
#                 one-standard-deviation shocks in the column order of the
#                 data. With the draw reordered so that the identification
#                 order is the data order, and the response divided by its
#                 impact on the shocked variable, it must equal the project's.
# It also summarises the responses over all draws: quantiles for the figures
# and how much the bands widen at long horizons (every draw has a largest
# companion root above 1; docs/methodology.md, section 3).
#
# Used by the tests (tests/testthat/test-irf.R).
#
# Inputs:  coefficient and covariance draws, a history matrix (the engine
#          needs one to condition; the responses do not depend on it).
# Outputs: data.frames of comparisons and summaries.


#' Effect of conditions on the first forecast quarter, by the scenario engine.
#'
#' Fixes the listed variables in T + 1 at the draw's no-shock path plus the
#' given deviations, leaves every other cell free, and returns the exact
#' conditional mean minus the no-shock path (Crump et al. 2021, Section 2.3,
#' the object D) over H + 1 quarters.
#'
#' @param B,Sigma One draw.
#' @param history Data matrix; its last p rows start the forecast.
#' @param lags,horizon Lags p and last IRF horizon H (the forecast runs over
#'   H + 1 quarters, T + 1 to T + H + 1).
#' @param deviations Named numeric vector: deviation from the no-shock path
#'   in T + 1 of each fixed variable.
#' @return Matrix (H + 1) x n, rows h0..hH.
irf_by_conditioning <- function(B, Sigma, history, lags, horizon, deviations) {
  n_q <- horizon + 1
  m <- forecast_unconditional(array(B, c(1, dim(B))), array(Sigma, c(1, dim(Sigma))), history, lags,
                              n_q, seed = 1, shocks = FALSE)[1, , ]
  cond <- matrix(NA_real_, n_q, ncol(history), dimnames = list(NULL, colnames(history)))
  j <- match(names(deviations), colnames(history))
  if (anyNA(j)) stop("unknown variable(s) in deviations")
  cond[1, j] <- m[1, j] + deviations
  effect <- conditional_moments_gaussian(B, Sigma, history, lags, n_q, cond)$mean - m
  dimnames(effect) <- list(paste0("h", 0:horizon), colnames(history))
  effect
}

#' Row of a comparison table.
irf_compare_row <- function(check, identification, diff, scale) {
  data.frame(check = check, identification = identification, max_abs_diff = max(abs(diff)),
             scale = scale, rel_diff = max(abs(diff)) / max(scale, .Machine$double.eps))
}

#' Closed-form IRFs against the scenario engine, one draw.
#'
#' @param B,Sigma One draw.
#' @param history,lags,horizon As in irf_by_conditioning().
#' @param girf_shocks Variables for generalized IRFs.
#' @param recursive List slow, policy, fast.
#' @return data.frame: check, identification, max_abs_diff, scale (largest
#'   absolute response), rel_diff; plus a row with the largest impact
#'   response of a slow variable in the recursive IRF (which must be 0).
irf_equivalence <- function(B, Sigma, history, lags, horizon, girf_shocks, recursive) {
  variables <- colnames(history)
  rows <- list()
  g <- girf(B, Sigma, lags, horizon, variables, girf_shocks)
  for (k in seq_along(girf_shocks)) {
    e <- irf_by_conditioning(B, Sigma, history, lags, horizon, setNames(1, girf_shocks[k]))
    rows[[k]] <- irf_compare_row("GIRF = conditioning", paste0("girf_", girf_shocks[k]),
                                 g[, , k] - e, max(abs(g[, , k])))
  }
  r <- recursive_irf(B, Sigma, lags, horizon, variables, recursive$slow, recursive$policy, recursive$fast)
  dev <- setNames(c(rep(0, length(recursive$slow)), 1), c(recursive$slow, recursive$policy))
  e <- irf_by_conditioning(B, Sigma, history, lags, horizon, dev)
  id <- paste0("recursive_", recursive$policy)
  rows[[length(rows) + 1]] <- irf_compare_row("recursive = conditioning", id, r - e, max(abs(r)))
  rows[[length(rows) + 1]] <- irf_compare_row("slow variables on impact", id, r[1, recursive$slow], 1)
  do.call(rbind, rows)
}

#' Recursive IRF under reorderings within the blocks, one draw.
#'
#' @param B,Sigma One draw.
#' @param lags,horizon Lags p and last horizon H.
#' @param variables Variable names, in model order.
#' @param recursive List slow, policy, fast.
#' @param perms Named list of functions; each takes a character vector and
#'   returns it reordered (applied to the slow and the fast block).
#' @return data.frame as irf_compare_row(), one row per reordering.
irf_block_invariance <- function(B, Sigma, lags, horizon, variables, recursive, perms) {
  base <- recursive_irf(B, Sigma, lags, horizon, variables, recursive$slow, recursive$policy, recursive$fast)
  rows <- lapply(names(perms), function(nm) {
    f <- perms[[nm]]
    r <- recursive_irf(B, Sigma, lags, horizon, variables, f(recursive$slow), recursive$policy,
                       f(recursive$fast))
    irf_compare_row("order within blocks", nm, r - base, max(abs(base)))
  })
  do.call(rbind, rows)
}

#' BVAR's Cholesky IRF for an identification order, rescaled to unit impact.
#'
#' The draw is reordered so that `order` is the data order
#' (permute_variables()); BVAR::irf() then gives responses to
#' one-standard-deviation shocks in that order. The response to the shock of
#' `shock` is divided by its impact on `shock` and put back in the model order.
#'
#' @param B,Sigma One draw.
#' @param history,lags Data matrix and lags (BVAR needs a data matrix).
#' @param horizon Last horizon H (BVAR is asked for H + 1 periods, 1 = impact).
#' @param order Identification order (all variable names).
#' @param shock Shocked variable.
#' @return Matrix (H + 1) x n in model order.
bvar_irf_unit <- function(B, Sigma, history, lags, horizon, order, shock) {
  variables <- colnames(history)
  perm <- match(order, variables)
  dummy <- matrix(NA_real_, 1, length(variables), dimnames = list(NULL, variables))
  pv <- permute_variables(B, Sigma, history, dummy, perm, lags)
  ir <- bvar_irf_one_draw(pv$B, pv$Sigma, pv$history, lags, horizon + 1)
  k <- match(shock, order)
  resp <- ir$irf[1, , , k]                              # response (permuted order) x horizon
  resp <- t(resp / resp[k, 1])                          # (H + 1) x n, unit impact on the shock
  out <- resp[, order(perm), drop = FALSE]              # back to the model order
  dimnames(out) <- list(paste0("h", 0:horizon), variables)
  out
}

#' Project IRFs against BVAR::irf(), one draw.
#'
#' @param B,Sigma One draw.
#' @param history,lags,horizon As in bvar_irf_unit().
#' @param girf_shocks,recursive As in irf_equivalence().
#' @return data.frame as irf_compare_row().
irf_bvar_agreement <- function(B, Sigma, history, lags, horizon, girf_shocks, recursive) {
  variables <- colnames(history)
  g <- girf(B, Sigma, lags, horizon, variables, girf_shocks)
  rows <- lapply(seq_along(girf_shocks), function(k) {
    j <- girf_shocks[k]
    b <- bvar_irf_unit(B, Sigma, history, lags, horizon, c(j, setdiff(variables, j)), j)
    irf_compare_row("project = BVAR irf()", paste0("girf_", j), g[, , k] - b, max(abs(g[, , k])))
  })
  order_r <- c(recursive$slow, recursive$policy, recursive$fast)
  r <- recursive_irf(B, Sigma, lags, horizon, variables, recursive$slow, recursive$policy, recursive$fast)
  b <- bvar_irf_unit(B, Sigma, history, lags, horizon, order_r, recursive$policy)
  rows[[length(rows) + 1]] <- irf_compare_row("project = BVAR irf()", paste0("recursive_", recursive$policy),
                                              r - b, max(abs(r)))
  do.call(rbind, rows)
}

#' How the bands of an IRF widen with the horizon.
#'
#' For each response: the 90 percent band width at selected horizons; the
#' share of draws whose absolute response at the last horizon exceeds that at
#' an earlier horizon; and the rank correlation, over draws, between the
#' largest companion root and the absolute response at the last horizon. With
#' roots above 1, responses of some draws grow without bound; these numbers
#' show how much that matters within the horizon.
#'
#' @param irf Array draws x (H + 1) x n.
#' @param max_root Largest companion modulus of each draw.
#' @param horizons Horizons for the band widths (the last one is H).
#' @param compare_h Earlier horizon for the share.
#' @param variables Responses to report.
#' @return data.frame per response.
irf_band_growth <- function(irf, max_root, horizons, compare_h, variables) {
  H <- dim(irf)[2] - 1
  rows <- lapply(variables, function(v) {
    x <- irf[, , v]
    w <- sapply(horizons, function(h) diff(quantile(x[, h + 1], c(0.05, 0.95), names = FALSE)))
    out <- data.frame(variable = v)
    for (k in seq_along(horizons)) out[[paste0("width90_h", horizons[k])]] <- w[k]
    out[[paste0("share_abs_h", H, "_gt_h", compare_h)]] <- mean(abs(x[, H + 1]) > abs(x[, compare_h + 1]))
    out$rank_cor_root <- cor(max_root, abs(x[, H + 1]), method = "spearman")
    out
  })
  do.call(rbind, rows)
}

#' Fan chart of impulse responses, one panel per response.
#'
#' Blue: posterior median, 68 and 90 percent bands over draws (parameter
#' uncertainty). Horizon 0 is the impact quarter.
#'
#' @param irf Array draws x (H + 1) x n.
#' @param variables Responses to show.
#' @param units Named character vector of unit labels per variable.
plot_irf_fans <- function(irf, variables, units) {
  H <- dim(irf)[2] - 1
  x <- 0:H
  for (v in variables) {
    q <- t(apply(irf[, , v], 2, quantile, c(0.05, 0.16, 0.5, 0.84, 0.95)))
    plot(x, q[, 3], type = "n", ylim = range(c(q, 0)), bty = "n", xlab = "quarters after impact",
         ylab = units[[v]], main = v, font.main = 1, cex.main = 0.95, col.main = ink[["primary"]],
         col.axis = ink[["secondary"]], col.lab = ink[["secondary"]])
    abline(h = 0, col = ink[["grid"]])
    polygon(c(x, rev(x)), c(q[, 1], rev(q[, 5])), col = "#2a78d633", border = NA)
    polygon(c(x, rev(x)), c(q[, 2], rev(q[, 4])), col = "#2a78d666", border = NA)
    lines(x, q[, 3], col = check_colours[["blue"]], lwd = 2)
  }
}

#' Legend for plot_irf_fans(), across the bottom of the figure.
legend_irf_fans <- function() {
  op <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  on.exit(par(op))
  plot.new()
  legend("bottom", horiz = TRUE, bty = "n", cex = 0.9, text.col = ink[["secondary"]],
         legend = c("90% band", "68% band", "median"),
         fill = c("#2a78d633", "#2a78d666", NA), border = NA,
         col = c(NA, NA, check_colours[["blue"]]), lty = c(NA, NA, 1), lwd = c(NA, NA, 2))
}
