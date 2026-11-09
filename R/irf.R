# irf.R: impulse responses from posterior draws of the VAR.
#
# Why it exists: an impulse response shows how the model's forecast of every
# variable changes, quarter by quarter, after a surprise in one variable.
# Crump et al. (2021), Section 2.3, define it as a scenario minus the
# baseline: for one posterior draw, the difference between the forecast mean
# given a condition on the first forecast quarter and the forecast mean
# without it. For the two identifications used in the paper this difference
# has a closed form, a Cholesky impulse response, which is what this file
# computes. The scenario engine (R/forecast_conditional.R) gives the
# same numbers by conditioning; tests/testthat/test-irf.R checks that.
#
# Notation (as in R/forecast.R): n variables, p lags, B (K x n) the
# coefficients, Sigma the covariance of the one-step forecast errors e_t,
# Psi_h the moving-average coefficients (Psi_0 = I; response of y_{t+h} to a
# unit forecast error e_t). A shock is summarised by its impact vector a:
# the forecast errors it causes in the first forecast quarter T + 1. The
# response h quarters after impact is
#   r_h = Psi_h a,   h = 0, 1, ..., H                                      (M27)
# h = 0 is the impact quarter T + 1. Since the VAR is linear, r_h does not
# depend on the history or on the forecast origin, only on (B, Sigma).
#
# Impact vectors:
#   Cholesky, ordering pi, shock to variable j (M27): order the variables by
#     pi, take the lower Cholesky factor L_pi of Sigma in that order
#     (L_pi L_pi' = Sigma_pi), and use the column of L_pi that belongs to j,
#     rescaled so that j itself moves by exactly 1 (a unit shock: 1 log point,
#     about 1 percent, for log series; 1 percentage point for rates).
#     Variables ordered before j do not move on impact.
#   Generalized (Koop, Pesaran and Potter 1996), variable j (M28):
#     a = Sigma e_j / Sigma_jj: every variable moves by its regression
#     coefficient on j's forecast error. This equals the Cholesky shock with j
#     ordered first, and the conditional forecast with y_{j,T+1} set to its
#     baseline + 1 and everything else free (Crump et al. 2021, Section 2.3).
#     It has no structural meaning: it is the mix of shocks that usually
#     moves j (Section 3.1).
#   Recursive (M29): ordering (slow block, policy variable, fast block),
#     shock to the policy variable. Slow variables do not move on impact,
#     fast ones may. It equals the conditional forecast with the slow
#     variables at their baseline and the policy variable at baseline + 1 in
#     T + 1, fast variables free (Section 2.3; Section 3.2, Figure 9). The
#     order within each block does not matter.
#
# The bands of an impulse response come from the posterior draws: one
# response per draw of (B, Sigma), no simulated shocks, as in the paper
# (Section 2.3: "for different draws of theta"). They show uncertainty
# about the parameters.
#
# The formulas are set out in docs/methodology.md (section 10); code comments
# cite its equation labels (M27)-(M29).
#
# Inputs:  coefficient draws (draws x K x n), covariance draws (draws x n x n),
#          the variable names, the horizon, the shocks (from config/run.yml,
#          section irf).
# Outputs: response draws, arrays draws x (H + 1) x n, in model units per
#          unit shock; quantile tables.


#' Responses to an impact vector over horizons 0..H.
#'
#' Implements docs/methodology.md (M27), r_h = Psi_h a.
#'
#' @param B Coefficients of one draw, K x n.
#' @param impact Impact vector a (length n, model order).
#' @param lags Number of lags p.
#' @param horizon Last horizon H (quarters after impact).
#' @param variables Variable names, in model order.
#' @return Matrix (H + 1) x n; row h + 1 is the response h quarters after impact.
propagate_impact <- function(B, impact, lags, horizon, variables) {
  Psi <- ma_coefficients(B, lags, horizon + 1)          # Psi_0, ..., Psi_H
  out <- matrix(NA_real_, horizon + 1, length(impact),
                dimnames = list(paste0("h", 0:horizon), variables))
  for (h in 0:horizon) out[h + 1, ] <- Psi[[h + 1]] %*% impact
  out
}

#' Impact vector of a Cholesky shock with unit impact.
#'
#' Implements docs/methodology.md (M27). The lower Cholesky factor of Sigma is
#' computed with the variables in the order `order`; the column of `shock`
#' is mapped back to the model order and divided by its own entry, so the
#' shocked variable moves by 1 and the variables ordered before it by 0.
#'
#' @param Sigma Covariance of one draw, n x n (model order).
#' @param variables Variable names, in model order.
#' @param order All variable names, in the identification order.
#' @param shock Name of the shocked variable.
#' @return Numeric vector of length n (model order), named.
cholesky_impact <- function(Sigma, variables, order, shock) {
  if (!setequal(order, variables) || length(order) != length(variables)) {
    stop("order must contain every variable exactly once")
  }
  if (!shock %in% variables) stop("unknown shock variable: ", shock)
  idx <- match(order, variables)                        # position in the model of each ordered variable
  L <- t(chol(Sigma[idx, idx]))                         # L L' = Sigma in the identification order
  k <- match(shock, order)
  a <- numeric(length(variables))
  a[idx] <- L[, k]                                      # back to the model order
  a <- a / a[match(shock, variables)]                   # unit impact on the shocked variable
  setNames(a, variables)
}

#' Cholesky impulse response with unit impact, one draw.
#'
#' Implements docs/methodology.md (M27).
#'
#' @param B,Sigma One draw (K x n, n x n).
#' @param lags,horizon Lags p and last horizon H.
#' @param variables Variable names, in model order.
#' @param order Identification order (all variables).
#' @param shock Shocked variable.
#' @return Matrix (H + 1) x n.
cholesky_irf <- function(B, Sigma, lags, horizon, variables, order, shock) {
  propagate_impact(B, cholesky_impact(Sigma, variables, order, shock), lags, horizon, variables)
}

#' Generalized impulse responses, one draw.
#'
#' Implements docs/methodology.md (M28): impact Sigma e_j / Sigma_jj for each
#' shocked variable j.
#'
#' @param B,Sigma One draw.
#' @param lags,horizon Lags p and last horizon H.
#' @param variables Variable names, in model order.
#' @param shocks Names of the shocked variables.
#' @return Array (H + 1) x n (response) x length(shocks) (shock).
girf <- function(B, Sigma, lags, horizon, variables, shocks) {
  unknown <- setdiff(shocks, variables)
  if (length(unknown)) stop("unknown shock variable(s): ", paste(unknown, collapse = ", "))
  Psi <- ma_coefficients(B, lags, horizon + 1)
  j <- match(shocks, variables)
  impact <- Sigma[, j, drop = FALSE] %*% diag(1 / diag(Sigma)[j], length(j))   # n x shocks
  out <- array(NA_real_, c(horizon + 1, length(variables), length(shocks)),
               dimnames = list(paste0("h", 0:horizon), variables, shocks))
  for (h in 0:horizon) out[h + 1, , ] <- Psi[[h + 1]] %*% impact
  out
}

#' Check the blocks of a recursive identification.
#'
#' The slow block, the policy variable and the fast block must not overlap
#' and must together contain every variable of the model once.
#'
#' @param slow,fast Character vectors of variable names.
#' @param policy Name of the policy variable.
#' @param variables Variable names of the model.
#' @return TRUE invisibly; stops with a message otherwise.
check_recursive_blocks <- function(slow, policy, fast, variables) {
  if (length(policy) != 1) stop("give exactly one policy variable")
  all_named <- c(slow, policy, fast)
  dup <- unique(all_named[duplicated(all_named)])
  if (length(dup)) stop("variable(s) in more than one block: ", paste(dup, collapse = ", "))
  unknown <- setdiff(all_named, variables)
  if (length(unknown)) stop("unknown variable(s) in the blocks: ", paste(unknown, collapse = ", "))
  missing_v <- setdiff(variables, all_named)
  if (length(missing_v)) stop("variable(s) in no block: ", paste(missing_v, collapse = ", "))
  invisible(TRUE)
}

#' Recursive impulse response to the policy variable, one draw.
#'
#' Implements docs/methodology.md (M29): Cholesky ordering (slow, policy,
#' fast), unit shock to the policy variable.
#'
#' @param B,Sigma One draw.
#' @param lags,horizon Lags p and last horizon H.
#' @param variables Variable names, in model order.
#' @param slow,policy,fast The blocks (see check_recursive_blocks()).
#' @return Matrix (H + 1) x n.
recursive_irf <- function(B, Sigma, lags, horizon, variables, slow, policy, fast) {
  check_recursive_blocks(slow, policy, fast, variables)
  cholesky_irf(B, Sigma, lags, horizon, variables, order = c(slow, policy, fast), shock = policy)
}

#' Impulse responses for every posterior draw.
#'
#' One response per draw of (B, Sigma), no simulated shocks: the spread over
#' draws is parameter uncertainty (Crump et al. 2021, Section 2.3).
#'
#' @param beta,sigma Draws (draws x K x n, draws x n x n).
#' @param lags,horizon Lags p and last horizon H.
#' @param variables Variable names, in model order.
#' @param girf_shocks Variables to shock one at a time (generalized IRFs).
#' @param recursive NULL, or a list with elements slow, policy, fast.
#' @return Named list of arrays draws x (H + 1) x n: "girf_<variable>" for
#'   each generalized IRF, "recursive_<policy>" for the recursive shock.
irf_draws <- function(beta, sigma, lags, horizon, variables, girf_shocks, recursive = NULL) {
  n_draws <- dim(beta)[1]
  if (dim(beta)[3] != length(variables)) stop("the draws have ", dim(beta)[3], " variables; ",
                                              length(variables), " names given")
  if (!is.null(recursive)) check_recursive_blocks(recursive$slow, recursive$policy, recursive$fast, variables)
  ids <- c(if (length(girf_shocks)) paste0("girf_", girf_shocks),
           if (!is.null(recursive)) paste0("recursive_", recursive$policy))
  out <- lapply(setNames(ids, ids), function(id) {
    array(NA_real_, c(n_draws, horizon + 1, length(variables)),
          dimnames = list(NULL, paste0("h", 0:horizon), variables))
  })
  for (i in seq_len(n_draws)) {
    B <- beta[i, , ]
    S <- sigma[i, , ]
    if (length(girf_shocks)) {
      g <- girf(B, S, lags, horizon, variables, girf_shocks)
      for (k in seq_along(girf_shocks)) out[[paste0("girf_", girf_shocks[k])]][i, , ] <- g[, , k]
    }
    if (!is.null(recursive)) {
      out[[paste0("recursive_", recursive$policy)]][i, , ] <-
        recursive_irf(B, S, lags, horizon, variables, recursive$slow, recursive$policy, recursive$fast)
    }
  }
  out
}

#' Posterior quantiles of an impulse response, per horizon and response.
#'
#' @param irf Array draws x (H + 1) x n (one element of irf_draws()).
#' @param probs Probabilities (config: quantiles).
#' @return data.frame: horizon (0..H), variable, one column per probability
#'   (q05, q16, ...), in model units per unit shock.
irf_quantiles <- function(irf, probs) {
  horizon <- dim(irf)[2] - 1
  variables <- dimnames(irf)[[3]]
  rows <- list()
  for (j in seq_along(variables)) {
    q <- t(apply(irf[, , j, drop = FALSE], 2, quantile, probs, names = FALSE))   # (H + 1) x probs
    colnames(q) <- sprintf("q%02d", round(100 * probs))
    rows[[j]] <- data.frame(horizon = 0:horizon, variable = variables[j], q)
  }
  do.call(rbind, rows)
}

#' IRF settings from the run configuration, checked.
#'
#' @param cfg Parsed config/run.yml.
#' @param variables Variable names of the model.
#' @return List: horizon, girf_shocks, recursive (slow, policy, fast).
irf_settings <- function(cfg, variables) {
  s <- cfg$irf
  if (is.null(s$horizon) || s$horizon < 1) stop("irf$horizon must be a positive number of quarters")
  girf_shocks <- unlist(s$girf_shocks)
  unknown <- setdiff(girf_shocks, variables)
  if (length(unknown)) stop("irf$girf_shocks: unknown variable(s): ", paste(unknown, collapse = ", "))
  rec <- s$recursive
  recursive <- if (is.null(rec)) NULL else list(slow = unlist(rec$slow), policy = rec$policy,
                                                fast = unlist(rec$fast))
  if (!is.null(recursive)) check_recursive_blocks(recursive$slow, recursive$policy, recursive$fast, variables)
  list(horizon = s$horizon, girf_shocks = girf_shocks, recursive = recursive)
}
