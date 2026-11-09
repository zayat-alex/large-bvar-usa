# irf_review.R: figures and tables for the review of the impulse responses
# (spec 9.5): signs and rough magnitudes of the responses, both COVID
# methods side by side.
#
# Why it exists: the impulse responses are checked for correctness against
# exact results in tests/testthat/test-irf.R. This review asks a different
# question: do the responses make economic sense, and how do they compare
# with Crump et al. (2021), Figures 6, 8 and 9? Each production run's engine
# output (outputs/<run_id>/engine/irf.rds) already holds the quantiles of
# every response; nothing is re-estimated.
#
# Units: unit shocks (1 log point, about 1 percent, for log series; 1
# percentage point for rates) in the first forecast quarter (horizon 0).
# Responses in model units: 100 x log difference from the no-shock path
# (about percent) for log series, percentage points for rates.
#
# Inputs:  irf.rds of each run (quantiles: horizon, variable, q05 ... q95).
# Outputs: figures (scripts/check_irf_review.R).


#' One shock's responses, all COVID methods, one panel per response.
#'
#' The first method is drawn with shaded 68 and 90 percent bands (blue), each
#' further method with lines (orange, then green: median solid, 68 percent
#' dashed, 90 percent dotted), so they can be compared directly.
#'
#' @param q Named list (one element per method) of quantile data.frames for
#'   one shock (irf.rds$quantiles[[shock]]).
#' @param responses Variables to show.
#' @param max_h Last horizon shown.
#' @param title Figure title (the shock).
#' @param scale Factor applied to every response (1 = unit shock).
plot_irf_methods <- function(q, responses, max_h, title, scale = 1) {
  n <- length(responses)
  n_col <- 4
  op <- par(mfrow = c(ceiling(n / n_col), n_col), mar = c(3, 3, 2.2, 0.8), oma = c(2.5, 0, 2.5, 0),
            col.axis = ink[["secondary"]], mgp = c(1.8, 0.5, 0))
  on.exit(par(op))
  for (v in responses) {
    d <- lapply(q, function(x) {
      x <- x[x$variable == v & x$horizon <= max_h, ]
      x[, grep("^q", names(x))] <- x[, grep("^q", names(x))] * scale
      x
    })
    ylim <- range(unlist(lapply(d, function(x) x[, c("q05", "q95")])), 0)
    h <- d[[1]]$horizon
    plot(NA, xlim = c(0, max_h), ylim = ylim, bty = "n", xlab = "quarters after impact", ylab = "",
         main = v, font.main = 1, cex.main = 1, col.main = ink[["primary"]], cex.axis = 0.8)
    abline(h = 0, col = ink[["secondary"]], lwd = 0.8)
    x <- d[[1]]
    polygon(c(h, rev(h)), c(x$q05, rev(x$q95)), col = "#2a78d633", border = NA)
    polygon(c(h, rev(h)), c(x$q16, rev(x$q84)), col = "#2a78d666", border = NA)
    lines(h, x$q50, col = check_colours[["blue"]], lwd = 2)
    line_cols <- c(check_colours[["orange"]], check_colours[["aqua"]])
    for (k in seq_along(d)[-1]) {
      x <- d[[k]]; cl <- line_cols[k - 1]
      lines(h, x$q50, col = cl, lwd = 2)
      lines(h, x$q16, col = cl, lty = 2)
      lines(h, x$q84, col = cl, lty = 2)
      lines(h, x$q05, col = cl, lty = 3)
      lines(h, x$q95, col = cl, lty = 3)
    }
  }
  mtext(title, outer = TRUE, line = 0.8, cex = 1.1, col = ink[["primary"]])
  others <- if (length(q) > 1) paste0(c("orange", "green")[seq_len(length(q) - 1)], " ", names(q)[-1],
                                      collapse = ", ")
  mtext(paste0("Blue (shaded 68% and 90%): ", names(q)[1],
               if (length(q) > 1) paste0(".  Lines (median, 68% dashed, 90% dotted): ", others), "."),
        side = 1, outer = TRUE, line = 1, cex = 0.8, col = ink[["secondary"]])
}


# ---- Credit-spread widening (a conditional forecast) ----
#
# The Baa yield usually moves with Treasury yields, so a generalized Baa
# shock is mostly a rise in all rates. A widening of the Baa spread is
# isolated by conditioning, in the first forecast quarter, on Baa up by
# `size` and on the held yield (gs10) at its baseline. The expected
# forecast errors given those two values are
#   a = Sigma[, S] Sigma[S, S]^{-1} r,   S = (shocked, held),  r = (size, 0),
# the conditional mean of a normal vector (as for the generalized IRF,
# docs/methodology.md (M28), with two conditions instead of one), and the
# response is propagated as r_h = Psi_h a (M27). This equals the difference
# between the conditional and the unconditional forecast means of a
# scenario with those two cells fixed (checked against the Gaussian
# conditioning code). It is the combination of shocks that most likely
# produces such a widening, not a structural credit shock.

#' Impact vector of a spread widening: one variable up, another held.
#'
#' @param Sigma Covariance of one draw (n x n, model order).
#' @param variables Variable names, in model order.
#' @param shocked Variable moved by `size` (e.g. "baa").
#' @param held Variable held at its baseline (e.g. "gs10").
#' @param size Move of `shocked`, model units.
#' @return Named impact vector (length n).
spread_impact <- function(Sigma, variables, shocked, held, size = 1) {
  S <- match(c(shocked, held), variables)
  if (anyNA(S)) stop("unknown variable in ", shocked, " / ", held)
  a <- Sigma[, S, drop = FALSE] %*% solve(Sigma[S, S], c(size, 0))
  setNames(as.vector(a), variables)
}

#' Spread-widening responses for every posterior draw.
#'
#' @param beta,sigma Draws (draws x K x n, draws x n x n).
#' @param lags,horizon Lags p and last horizon H.
#' @param variables Variable names, in model order.
#' @param shocked,held,size As in spread_impact().
#' @return Array draws x (H + 1) x n.
spread_irf_draws <- function(beta, sigma, lags, horizon, variables, shocked, held, size = 1) {
  out <- array(NA_real_, c(dim(beta)[1], horizon + 1, length(variables)),
               dimnames = list(NULL, paste0("h", 0:horizon), variables))
  for (i in seq_len(dim(beta)[1])) {
    a <- spread_impact(sigma[i, , ], variables, shocked, held, size)
    out[i, , ] <- propagate_impact(beta[i, , ], a, lags, horizon, variables)
  }
  out
}

#' Cross-check: spread response = conditional minus unconditional mean.
#'
#' For one draw, the scenario fixes `shocked` at its no-shock forecast plus
#' `size` and `held` at its no-shock forecast in the first quarter; the
#' difference of the exact conditional and unconditional means
#' (conditional_moments_gaussian()) must equal the propagated impact.
#'
#' @return Largest absolute difference over horizons 0..H-1 and variables.
spread_vs_conditional <- function(B, Sigma, history, lags, horizon, shocked, held, size = 1) {
  vars <- colnames(history)
  free <- matrix(NA_real_, horizon, length(vars), dimnames = list(NULL, vars))
  base <- conditional_moments_gaussian(B, Sigma, history, lags, horizon, free)$mean
  cond <- free
  cond[1, shocked] <- base[1, shocked] + size
  cond[1, held]    <- base[1, held]
  scen <- conditional_moments_gaussian(B, Sigma, history, lags, horizon, cond)$mean
  irf <- propagate_impact(B, spread_impact(Sigma, vars, shocked, held, size), lags, horizon - 1, vars)
  max(abs((scen - base) - irf))
}


# ---- Expectations of spec 9.5 ----

#' Result of each expected response, from the response quantiles.
#'
#' For each row: within the horizons of the window, the horizon at which the
#' median moves furthest in the expected direction; the median and the 68
#' percent band there (scaled to the shock size of the row). Verdict: "as
#' expected" if the sign matches and the 68 percent band excludes zero;
#' "weak" if the sign matches but the band includes zero; "opposite" if the
#' median never moves in the expected direction.
#'
#' @param q Quantile data.frame of one shock (horizon, variable, q05 ... q95).
#' @param rows data.frame: variable, sign (+1 or -1), h_from, h_to, scale.
#' @return rows with horizon, median, q16, q84 (scaled), verdict.
expectation_results <- function(q, rows) {
  out <- lapply(seq_len(nrow(rows)), function(k) {
    r <- rows[k, ]
    d <- q[q$variable == r$variable & q$horizon >= r$h_from & q$horizon <= r$h_to, ]
    s <- r$scale
    m <- d$q50 * s
    lo <- pmin(d$q16 * s, d$q84 * s); hi <- pmax(d$q16 * s, d$q84 * s)
    j <- which.max(r$sign * m)
    sign_ok <- r$sign * m[j] > 0
    excl <- (lo[j] > 0) || (hi[j] < 0)
    data.frame(r, horizon = d$horizon[j], median = m[j], q16 = lo[j], q84 = hi[j],
               verdict = if (!sign_ok) "opposite" else if (excl) "as expected" else "weak")
  })
  do.call(rbind, out)
}

#' Quantiles of the effect on YoY growth at chosen horizons.
#'
#' @param draws Response draws (draws x (H + 1) x n).
#' @param variables Variables to report.
#' @param is_log Logical, one per model variable (config/series.csv transform).
#' @param horizons Horizons to report.
#' @param scale Shock size factor.
#' @param probs Quantiles.
#' @return data.frame: variable, horizon, measure, q16, q50, q84.
irf_growth_table <- function(draws, variables, is_log, horizons, scale = 1, probs = c(0.16, 0.5, 0.84)) {
  g <- irf_annual_measure(draws, is_log) * scale
  rows <- list()
  for (v in variables) {
    j <- match(v, dimnames(draws)[[3]])
    for (h in horizons) {
      x <- quantile(g[, h + 1, j], probs, names = FALSE)
      rows[[length(rows) + 1]] <- data.frame(variable = v, horizon = h,
                                             measure = if (is_log[j]) "YoY growth, pp" else "level, pp",
                                             q16 = min(x[1], x[3]), q50 = x[2], q84 = max(x[1], x[3]))
    }
  }
  do.call(rbind, rows)
}

#' Our statistic for one row of the paper comparison.
#'
#' @param q Quantile data.frame of one shock.
#' @param variable,stat,h_from,h_to,scale Row fields: stat is "max", "min" or
#'   "at" (the median at h_from).
#' @return List: value (median, scaled), horizon.
paper_stat <- function(q, variable, stat, h_from, h_to, scale) {
  d <- q[q$variable == variable & q$horizon >= h_from & q$horizon <= h_to, ]
  m <- d$q50 * scale
  j <- switch(stat, max = which.max(m), min = which.min(m), at = 1L)
  list(value = m[j], horizon = d$horizon[j])
}
