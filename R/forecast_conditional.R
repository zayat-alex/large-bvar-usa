# forecast_conditional.R: conditional forecasts (scenarios) from posterior
# draws of the VAR.
#
# Why it exists: a scenario fixes future values of some variables in some
# quarters and asks what the model implies for everything else (Crump et al.
# 2021, Section 2.3). For each posterior draw of (B, Sigma) the forecast is
# drawn from the predictive density given the fixed cells. Two equivalent
# methods compute that conditional distribution, chosen with
# forecast_conditional(method = ...):
#   1. direct conditioning of a Gaussian vector (Waggoner and Zha 1999);
#   2. the Kalman filter and simulation smoother, the paper's algorithm
#      (Banbura, Giannone and Lenza 2015).
# They share no computational code, so that each checks the other; with the
# same seed they give the same paths up to rounding
# (reports/conditional_checks.qmd).
#
# Contents: the input check and the function that dispatches to a method;
# then the Gaussian method; then the Kalman method. Each method's section
# starts with its derivation.
#
# Notation as in R/forecast.R (unconditional forecasts): n variables, p
# lags, K = 1 + n p, B (K x n) with the constant in row 1.
#
# Inputs:  coefficient draws (draws x K x n), covariance draws (draws x n x n),
#          a history matrix (its last p rows start the forecast), the horizon,
#          a conditions matrix (H x n, model units, NA = free), a seed.
# Outputs: conditional forecast draws ((draws x paths) x H x n), or the
#          exact conditional moments of one draw (for checks).


# ==== Choosing the method ====

#' Check a conditions matrix against the history and the horizon.
#'
#' A conditions matrix has one row per forecast quarter (row h is T + h) and
#' one column per variable, named and ordered as the history; a number fixes
#' that cell (in model units), NA leaves it free. The rules of a scenario
#' file (which quarters and variables may be fixed) are checked when the
#' file is read; this check only guards the forecasting functions.
#'
#' @param conditions Matrix H x n.
#' @param history Data matrix whose columns name the variables.
#' @param horizon Number of quarters H.
#' @return TRUE invisibly; stops with a message otherwise.
check_conditions <- function(conditions, history, horizon) {
  if (!is.matrix(conditions) || !is.numeric(conditions)) stop("conditions must be a numeric matrix")
  if (nrow(conditions) != horizon) {
    stop("conditions has ", nrow(conditions), " rows; expected one per forecast quarter (", horizon, ")")
  }
  if (!identical(colnames(conditions), colnames(history))) {
    stop("the columns of conditions must be the variables of history, in the same order")
  }
  if (any(is.infinite(conditions) | is.nan(conditions))) stop("conditions contains Inf or NaN")
  invisible(TRUE)
}

#' Conditional forecast paths (a scenario), one per posterior draw.
#'
#' Draws from the predictive density given fixed future values of some
#' variables in some quarters (Crump et al. 2021, Section 2.3). Both methods
#' compute the same conditional distribution; they differ only in the
#' algorithm, and with the same seed they give the same paths up to rounding
#' (reports/conditional_checks.qmd).
#'
#' @param beta,sigma Draws (draws x K x n, draws x n x n).
#' @param history Data matrix; its last p rows start the forecast.
#' @param lags,horizon Lags p and horizon H.
#' @param conditions H x n matrix in model units, NA = free cell (see
#'   check_conditions()).
#' @param seed Seed of the random number generator.
#' @param method "kalman" (Kalman simulation smoother, the paper's
#'   algorithm) or "gaussian" (direct conditioning, Waggoner and Zha 1999).
#' @param shocks TRUE for draws; FALSE for each draw's conditional mean.
#' @param n_paths Paths per posterior draw.
#' @return Array (draws x n_paths) x H x n.
forecast_conditional <- function(beta, sigma, history, lags, horizon, conditions, seed, method,
                                 shocks = TRUE, n_paths = 1) {
  if (missing(method)) stop("give the method: \"kalman\" or \"gaussian\"")
  f <- switch(method,
              kalman = forecast_conditional_kalman,
              gaussian = forecast_conditional_gaussian,
              stop("unknown conditional forecast method: ", method, " (use \"kalman\" or \"gaussian\")"))
  f(beta, sigma, history, lags, horizon, conditions, seed, shocks = shocks, n_paths = n_paths)
}


# ==== Method 1: direct Gaussian conditioning ====
#
# Conditional forecasts by direct conditioning of a Gaussian vector (the
# Waggoner and Zha 1999 form).
#
# Why it exists: a scenario fixes future values of some variables in some
# quarters and asks what the model implies for everything else. Given one
# posterior draw of (B, Sigma), the future path over H quarters is a linear
# function of the future shocks, so it is jointly Gaussian; fixing some of
# its cells and drawing the rest is conditioning a Gaussian vector on some of
# its elements, which has a closed form. This section computes that closed form.
# It is one of two implementations of the same conditional distribution; the
# other is the Kalman simulation smoother of Crump et al. (2021), Section 2.3
# (the Kalman section below). The two share no computational code, so
# that each checks the other.
#
# The idea, for one draw (notation as in R/forecast.R):
#   Stack the future path y = (y_{T+1}', ..., y_{T+H}')' (n H x 1). Then
#     y = m + A u,   u ~ N(0, I_{nH}),                                  (M21)
#   where m is the no-shock path, and A = Theta (I_H (x) L) with L L' = Sigma
#   and Theta block lower triangular, block (h, k) = Psi_{h-k}: Psi_j is the
#   response of y after j quarters to a unit innovation (Psi_0 = I,
#   Psi_j = sum_{s=1}^{min(j,p)} B_s Psi_{j-s}).
#   A scenario fixes q cells: R y = r, with R selecting them. In terms of
#   the standardized shocks this is G u = d, with G = R A (the rows of A of
#   the fixed cells) and d = r - R m. Then (M22)
#     E[y | R y = r]   = m + A G'(G G')^{-1} d,
#     Var[y | R y = r] = A (I - P) A',   P = G'(G G')^{-1} G,
#   and a draw from the conditional distribution is (M23)
#     u = G'(G G')^{-1} d + (I - P) z,   z ~ N(0, I),   y = m + A u.
#   Since A is block lower triangular with L (invertible) on its diagonal, G
#   has full row rank whenever the fixed cells are distinct, so G G' can be
#   inverted for any scenario that fixes each cell at most once.
#
# The formulas are set out in docs/methodology.md (section 9); code comments
# cite its equation labels (M21), (M22), ...
#
# Inputs:  coefficient draws (draws x K x n), covariance draws (draws x n x n),
#          a history matrix (its last p rows start the forecast), the horizon,
#          a conditions matrix (H x n, model units, NA = free), a seed.
# Outputs: conditional forecast draws (draws x H x n), or the exact
#          conditional moments of one draw (for checks).


#' Moving-average coefficients Psi_0, ..., Psi_{H-1} of a VAR.
#'
#' Psi_j is the response of y_{t+j} to a unit innovation e_t: Psi_0 = I and
#' Psi_j = sum_{s=1}^{min(j, p)} B_s Psi_{j-s}, where (B_s)_ij is the effect
#' of variable j at lag s on variable i, read from row 1 + (s - 1) n + j,
#' column i of B.
#'
#' @param B Coefficients, K x n.
#' @param lags Number of lags p.
#' @param horizon Number of quarters H.
#' @return List of H matrices n x n; element j + 1 is Psi_j.
ma_coefficients <- function(B, lags, horizon) {
  n_var <- ncol(B)
  B_lag <- lapply(seq_len(lags), function(s) t(B[1 + (s - 1) * n_var + seq_len(n_var), , drop = FALSE]))
  Psi <- vector("list", horizon)
  Psi[[1]] <- diag(n_var)
  for (j in seq_len(horizon - 1)) {
    acc <- matrix(0, n_var, n_var)
    for (s in seq_len(min(j, lags))) acc <- acc + B_lag[[s]] %*% Psi[[j - s + 1]]
    Psi[[j + 1]] <- acc
  }
  Psi
}

#' The loading matrix A of the stacked future path on standardized shocks.
#'
#' Implements docs/methodology.md (M21): A = Theta (I_H (x) L), so block
#' (h, k) of A, for k <= h, is Psi_{h-k} L; blocks above the diagonal are 0
#' (a shock in quarter k cannot move earlier quarters).
#'
#' @param Psi Output of ma_coefficients().
#' @param L Lower Cholesky factor of Sigma (L L' = Sigma).
#' @return Matrix nH x nH; rows and columns ordered quarter by quarter
#'   (all n variables of T+1, then of T+2, ...).
stacked_shock_loading <- function(Psi, L) {
  n_var <- nrow(L)
  horizon <- length(Psi)
  A <- matrix(0, n_var * horizon, n_var * horizon)
  for (h in seq_len(horizon)) for (k in seq_len(h)) {
    rows <- (h - 1) * n_var + seq_len(n_var)
    cols <- (k - 1) * n_var + seq_len(n_var)
    A[rows, cols] <- Psi[[h - k + 1]] %*% L
  }
  A
}

#' The pieces of the conditional distribution for one draw.
#'
#' Computes m, A and, if there are conditions, the quantities of (M22)-(M23)
#' that do not depend on z: the fixed cells, G, the Cholesky factor of G G',
#' and the conditional mean of the standardized shocks, G'(G G')^{-1} d.
#'
#' @param B,Sigma One draw (K x n, n x n).
#' @param history Data matrix; its last p rows start the forecast.
#' @param lags,horizon Lags p and horizon H.
#' @param conditions H x n matrix, NA = free cell.
#' @return List: m (stacked no-shock path, nH), A, fixed (positions of the
#'   fixed cells in the stacked vector), G, G_chol, u_mean.
gaussian_conditioning_setup <- function(B, Sigma, history, lags, horizon, conditions) {
  # Step 1: the no-shock path m: y_t = B' x_t with all future shocks at zero.
  # (Computed here rather than with forecast_unconditional(), which sets the
  # seed and would restart the random numbers of the calling loop.)
  n_var  <- ncol(B)
  lagged <- history[nrow(history) - seq_len(lags) + 1, , drop = FALSE]   # row s = lag s
  m <- numeric(0)
  for (h in seq_len(horizon)) {
    y <- as.vector(crossprod(B, c(1, as.vector(t(lagged)))))
    m <- c(m, y)                                   # stacked quarter by quarter
    lagged <- rbind(y, lagged[-lags, , drop = FALSE])
  }
  # Step 2: the loading A on standardized shocks (M21).
  A <- stacked_shock_loading(ma_coefficients(B, lags, horizon), t(chol(Sigma)))
  # Step 3: the fixed cells, in the same stacked order (t() puts quarter first).
  stacked_conditions <- as.vector(t(conditions))
  fixed <- which(!is.na(stacked_conditions))
  out <- list(m = m, A = A, fixed = fixed, G = NULL, G_chol = NULL, u_mean = rep(0, length(m)))
  if (length(fixed) == 0) return(out)
  # Step 4: G u = d with G = rows of A, d = r - R m; u_mean = G'(G G')^{-1} d (M22).
  G <- A[fixed, , drop = FALSE]
  d <- stacked_conditions[fixed] - m[fixed]
  G_chol <- chol(tcrossprod(G))                    # upper triangular, G G' = R'R
  out$G <- G
  out$G_chol <- G_chol
  out$u_mean <- as.vector(crossprod(G, chol2inv(G_chol) %*% d))
  out
}

#' Conditional forecast paths, one or more per posterior draw (Gaussian method).
#'
#' For each draw: the setup of gaussian_conditioning_setup(), then for each
#' path u = u_mean + (I - P) z, y = m + A u (M23). With shocks = FALSE it
#' returns each draw's conditional mean m + A u_mean (M22). With no fixed
#' cells it returns the unconditional forecast; with the same seed and one
#' path per draw it uses the same random numbers as forecast_unconditional().
#'
#' Implements docs/methodology.md (M21)-(M23).
#'
#' @param beta Coefficient draws, array draws x K x n.
#' @param sigma Covariance draws, array draws x n x n.
#' @param history Numeric matrix of observed data (quarters x n), in the
#'   variable order of the draws; its last p rows start the forecast.
#' @param lags Number of lags p.
#' @param horizon Number of quarters H.
#' @param conditions H x n matrix in model units, column names as in
#'   `history`; NA marks a free cell.
#' @param seed Seed of the random number generator.
#' @param shocks TRUE for draws from the conditional predictive density;
#'   FALSE for each draw's conditional mean.
#' @param n_paths Paths per posterior draw (more than 1 only for checks; the
#'   setup is computed once per draw and reused).
#' @return Array (draws x n_paths) x H x n, paths of draw 1 first.
forecast_conditional_gaussian <- function(beta, sigma, history, lags, horizon, conditions, seed,
                                          shocks = TRUE, n_paths = 1) {
  check_conditions(conditions, history, horizon)
  n_draws <- dim(beta)[1]
  n_var   <- dim(beta)[3]
  if (!shocks) n_paths <- 1
  out <- array(NA_real_, c(n_draws * n_paths, horizon, n_var),
               dimnames = list(NULL, paste0("h", seq_len(horizon)), colnames(history)))
  set.seed(seed)
  row <- 0
  for (i in seq_len(n_draws)) {
    s <- gaussian_conditioning_setup(beta[i, , ], sigma[i, , ], history, lags, horizon, conditions)
    for (k in seq_len(n_paths)) {
      u <- s$u_mean
      if (shocks) {
        z <- rnorm(n_var * horizon)
        # (I - P) z = z - G'(G G')^{-1} G z: remove the part of z that would move the fixed cells.
        if (!is.null(s$G)) z <- z - as.vector(crossprod(s$G, chol2inv(s$G_chol) %*% (s$G %*% z)))
        u <- u + z
      }
      y <- s$m + as.vector(s$A %*% u)
      row <- row + 1
      out[row, , ] <- matrix(y, horizon, n_var, byrow = TRUE)
    }
  }
  out
}

#' Exact conditional mean and covariance of one draw (Gaussian method).
#'
#' Implements docs/methodology.md (M22). Used by checks; the covariance is
#' the full nH x nH matrix A (I - P) A'.
#'
#' @param B,Sigma One draw (K x n, n x n).
#' @param history,lags,horizon,conditions As in forecast_conditional_gaussian().
#' @return List: mean (H x n), cov (nH x nH, stacked quarter by quarter),
#'   var (list of H matrices n x n: the covariance within each quarter).
conditional_moments_gaussian <- function(B, Sigma, history, lags, horizon, conditions) {
  check_conditions(conditions, history, horizon)
  n_var <- ncol(B)
  s <- gaussian_conditioning_setup(B, Sigma, history, lags, horizon, conditions)
  mean_stacked <- s$m + as.vector(s$A %*% s$u_mean)
  # A (I - P) = A - (A G') (G G')^{-1} G; since I - P is symmetric and
  # idempotent, A (I - P) A' = [A (I - P)] [A (I - P)]'.
  A_free <- s$A
  if (!is.null(s$G)) A_free <- s$A - (s$A %*% t(s$G)) %*% chol2inv(s$G_chol) %*% s$G
  cov <- tcrossprod(A_free)
  var <- lapply(seq_len(horizon), function(h) {
    idx <- (h - 1) * n_var + seq_len(n_var)
    cov[idx, idx]
  })
  list(mean = matrix(mean_stacked, horizon, n_var, byrow = TRUE,
                     dimnames = list(paste0("h", seq_len(horizon)), colnames(history))),
       cov = cov, var = var)
}


# ==== Method 2: Kalman filter and simulation smoother ====
#
# Conditional forecasts with the Kalman filter and the simulation smoother,
# as in Crump et al. (2021), Section 2.3, following Banbura, Giannone and
# Lenza (2015).
#
# Why it exists: this is the paper's algorithm for scenarios. The VAR is
# written as a state-space model; the future quarters are appended to the
# sample, the cells a scenario fixes are "observed" and all other future
# cells are missing. The Kalman filter and smoother then give the mean of
# the future path given the fixed cells, and the simulation smoother of
# Durbin and Koopman (2002) gives draws from its conditional distribution.
# The same conditional distribution is computed in closed form in
# the Gaussian section above; the two sections share no computational
# code, so that each checks the other.
#
# State space (paper eqs. 6-8), for one posterior draw:
#   y_t     = C_t S_t                       (measurement: C_t selects the
#                                            observed cells of quarter t)
#   S_{t+1} = G S_{t} + w_t,  w_t ~ N(0, Q)                              (M24)
#   S_t = (y_t', y_{t-1}', ..., y_{t-p+1}', c')'    (n p + n elements)
#   G   = [ B_1 B_2 ... B_p  I_n ]     (first n rows: the VAR equation)
#         [ I_n(p-1)    0     0  ]     (shift the lags down by one)
#         [ 0    ...    0    I_n ]     (the constant stays constant)
#   Q   = Sigma in the top-left n x n block, 0 elsewhere.
# There is no measurement error: a fixed cell is matched exactly.
#
# Forecast quarters are t = T+1, ..., T+H, written t = 1, ..., H below. The
# state at T is known (the last p observations and the draw's constant), so
# the filter starts at a_1 = G S_T with P_1 = Q.
#
# Kalman filter and state smoother (Durbin and Koopman 2012, Sections 4.3
# and 4.4), with Z_t the rows of the identity for the cells observed at t
# (M25):
#   v_t = y_t - Z_t a_t,  F_t = Z_t P_t Z_t',  K_t = G P_t Z_t' F_t^{-1},
#   L_t = G - K_t Z_t,    a_{t+1} = G a_t + K_t v_t,
#   P_{t+1} = G P_t L_t' + Q;
#   backwards from r_H = 0, N_H = 0:
#   r_{t-1} = Z_t' F_t^{-1} v_t + L_t' r_t,
#   N_{t-1} = Z_t' F_t^{-1} Z_t + L_t' N_t L_t,
#   E[S_t | fixed cells] = a_t + P_t r_{t-1},
#   Var[S_t | fixed cells] = P_t - P_t N_{t-1} P_t.
# In a quarter with no fixed cells, Z_t is empty: K_t = 0, L_t = G.
#
# Simulation smoother by mean correction (Durbin and Koopman 2002) (M26):
#   1. simulate an unconditional path S^+ from the known state at T, and
#      read off its values y^+ at the fixed cells;
#   2. smooth twice: S_hat from the scenario's values, S_hat^+ from y^+;
#   3. the draw is S_tilde = S_hat - S_hat^+ + S^+.
# The covariances P_t, F_t, K_t, L_t do not depend on the data, so they are
# computed once per posterior draw and reused for both smoothing passes and
# for all paths of that draw.
#
# The formulas are set out in docs/methodology.md (section 9); code comments
# cite its equation labels (M24), (M25), (M26).
#
# Inputs:  coefficient draws (draws x K x n), covariance draws (draws x n x n),
#          a history matrix (its last p rows start the forecast), the horizon,
#          a conditions matrix (H x n, model units, NA = free), a seed.
# Outputs: conditional forecast draws (draws x H x n), or the conditional
#          means and per-quarter covariances of one draw (for checks).


#' State-space matrices of the VAR (paper eqs. 6-8).
#'
#' Implements docs/methodology.md (M24).
#'
#' @param B,Sigma One draw (K x n, n x n).
#' @param history Data matrix; its last p rows give the known state at T.
#' @param lags Number of lags p.
#' @return List: G (transition, m x m with m = n p + n), Q (state shock
#'   covariance, m x m), S_T (known state at T), n_var.
var_state_space <- function(B, Sigma, history, lags) {
  n_var <- ncol(B)
  n_lag_states <- n_var * lags
  m <- n_lag_states + n_var
  G <- matrix(0, m, m)
  # First n rows: y_{t+1} = B_1 y_t + ... + B_p y_{t-p+1} + c, with
  # (B_s)_ij = B[1 + (s - 1) n + j, i]; the rows 2..K of B, transposed, are
  # [B_1 ... B_p].
  G[seq_len(n_var), seq_len(n_lag_states)] <- t(B[-1, , drop = FALSE])
  G[seq_len(n_var), n_lag_states + seq_len(n_var)] <- diag(n_var)
  # Lags shift down: the new lag s + 1 is the old lag s.
  if (lags > 1) {
    for (s in seq_len(lags - 1)) {
      G[s * n_var + seq_len(n_var), (s - 1) * n_var + seq_len(n_var)] <- diag(n_var)
    }
  }
  # The constant carries over unchanged.
  G[n_lag_states + seq_len(n_var), n_lag_states + seq_len(n_var)] <- diag(n_var)
  Q <- matrix(0, m, m)
  Q[seq_len(n_var), seq_len(n_var)] <- Sigma
  # Known state at T: y_T, y_{T-1}, ..., y_{T-p+1}, then the constant c (row 1 of B).
  lagged <- history[nrow(history) - seq_len(lags) + 1, , drop = FALSE]
  S_T <- c(as.vector(t(lagged)), B[1, ])
  list(G = G, Q = Q, S_T = S_T, n_var = n_var)
}

#' The data-independent part of the Kalman filter: P_t, F_t^{-1}, K_t, L_t.
#'
#' Implements the covariance recursions of docs/methodology.md (M25),
#' starting from P_1 = Q (the state at T is known).
#'
#' @param ss Output of var_state_space().
#' @param observed List of H integer vectors: the variables fixed in each
#'   quarter (positions 1..n, which are also their positions in the state).
#' @return List of H elements, each a list: obs, P, F_inv, K, L.
kalman_covariances <- function(ss, observed) {
  G <- ss$G
  P <- ss$Q
  out <- vector("list", length(observed))
  for (t in seq_along(observed)) {
    obs <- observed[[t]]
    GP <- G %*% P
    if (length(obs) > 0) {
      F_inv <- chol2inv(chol(P[obs, obs, drop = FALSE]))      # F_t = Z P Z'
      K <- GP[, obs, drop = FALSE] %*% F_inv                  # K_t = G P Z' F^{-1}
      L <- G
      L[, obs] <- L[, obs] - K                                # L_t = G - K Z
    } else {
      F_inv <- NULL
      K <- NULL
      L <- G
    }
    out[[t]] <- list(obs = obs, P = P, F_inv = F_inv, K = K, L = L)
    P <- GP %*% t(L) + ss$Q                                   # P_{t+1} = G P L' + Q
  }
  out
}

#' Smoothed state means given values at the observed cells.
#'
#' Runs the mean part of the Kalman filter forwards and the state smoother
#' backwards (docs/methodology.md (M25)), using the covariances of
#' kalman_covariances().
#'
#' @param ss Output of var_state_space().
#' @param kc Output of kalman_covariances().
#' @param values List of H numeric vectors: the values of the observed cells
#'   in each quarter, in the order of kc[[t]]$obs.
#' @return Matrix H x m: E[S_t | observed cells], t = 1..H.
kalman_smoothed_means <- function(ss, kc, values) {
  H <- length(kc)
  m <- length(ss$S_T)
  a <- matrix(NA_real_, H, m)
  v <- vector("list", H)
  # Forward: a_1 = G S_T; a_{t+1} = G a_t + K_t v_t.
  a_t <- as.vector(ss$G %*% ss$S_T)
  for (t in seq_len(H)) {
    a[t, ] <- a_t
    obs <- kc[[t]]$obs
    a_next <- as.vector(ss$G %*% a_t)
    if (length(obs) > 0) {
      v[[t]] <- values[[t]] - a_t[obs]
      a_next <- a_next + as.vector(kc[[t]]$K %*% v[[t]])
    }
    a_t <- a_next
  }
  # Backward: r_{t-1} = Z' F^{-1} v_t + L_t' r_t;  E[S_t | .] = a_t + P_t r_{t-1}.
  smoothed <- matrix(NA_real_, H, m)
  r <- rep(0, m)
  for (t in rev(seq_len(H))) {
    obs <- kc[[t]]$obs
    r_prev <- as.vector(crossprod(kc[[t]]$L, r))
    if (length(obs) > 0) r_prev[obs] <- r_prev[obs] + as.vector(kc[[t]]$F_inv %*% v[[t]])
    smoothed[t, ] <- a[t, ] + as.vector(kc[[t]]$P %*% r_prev)
    r <- r_prev
  }
  smoothed
}

#' Smoothed state covariances, Var[S_t | observed cells].
#'
#' The N_t recursion of docs/methodology.md (M25). Only used by checks: the
#' draws do not need it.
#'
#' @param kc Output of kalman_covariances().
#' @return List of H matrices m x m.
kalman_smoothed_variances <- function(kc) {
  H <- length(kc)
  m <- nrow(kc[[1]]$P)
  N <- matrix(0, m, m)
  out <- vector("list", H)
  for (t in rev(seq_len(H))) {
    obs <- kc[[t]]$obs
    N_prev <- crossprod(kc[[t]]$L, N %*% kc[[t]]$L)            # L' N L
    if (length(obs) > 0) N_prev[obs, obs] <- N_prev[obs, obs] + kc[[t]]$F_inv   # + Z' F^{-1} Z
    P <- kc[[t]]$P
    out[[t]] <- P - P %*% N_prev %*% P
    N <- N_prev
  }
  out
}

#' Fixed cells of a conditions matrix, quarter by quarter.
#'
#' @param conditions H x n matrix, NA = free.
#' @return List: observed (list of H integer vectors of variable positions),
#'   values (list of H numeric vectors).
kalman_observed_cells <- function(conditions) {
  observed <- lapply(seq_len(nrow(conditions)), function(t) which(!is.na(conditions[t, ])))
  values <- lapply(seq_len(nrow(conditions)), function(t) unname(conditions[t, observed[[t]]]))
  list(observed = observed, values = values)
}

#' Conditional forecast paths, one or more per posterior draw (Kalman
#' simulation smoother).
#'
#' For each draw: the state space (M24), the filter covariances (M25), the
#' smoothed means given the scenario's values; then for each path the mean
#' correction of Durbin and Koopman (2002) (M26). With shocks = FALSE it
#' returns each draw's conditional mean (the smoothed means). With no fixed
#' cells the smoothed means are the no-shock path and every draw is the
#' simulated unconditional path; with the same seed and one path per draw
#' the random numbers are those of forecast_unconditional().
#'
#' Implements docs/methodology.md (M24)-(M26).
#'
#' @param beta Coefficient draws, array draws x K x n.
#' @param sigma Covariance draws, array draws x n x n.
#' @param history Numeric matrix of observed data (quarters x n), in the
#'   variable order of the draws; its last p rows start the forecast.
#' @param lags Number of lags p.
#' @param horizon Number of quarters H.
#' @param conditions H x n matrix in model units, column names as in
#'   `history`; NA marks a free cell.
#' @param seed Seed of the random number generator.
#' @param shocks TRUE for draws from the conditional predictive density;
#'   FALSE for each draw's conditional mean.
#' @param n_paths Paths per posterior draw (more than 1 only for checks; the
#'   filter covariances are computed once per draw and reused).
#' @return Array (draws x n_paths) x H x n, paths of draw 1 first.
forecast_conditional_kalman <- function(beta, sigma, history, lags, horizon, conditions, seed,
                                        shocks = TRUE, n_paths = 1) {
  check_conditions(conditions, history, horizon)
  n_draws <- dim(beta)[1]
  n_var   <- dim(beta)[3]
  if (!shocks) n_paths <- 1
  cells <- kalman_observed_cells(conditions)
  out <- array(NA_real_, c(n_draws * n_paths, horizon, n_var),
               dimnames = list(NULL, paste0("h", seq_len(horizon)), colnames(history)))
  set.seed(seed)
  row <- 0
  for (i in seq_len(n_draws)) {
    ss <- var_state_space(beta[i, , ], sigma[i, , ], history, lags)
    kc <- kalman_covariances(ss, cells$observed)
    S_hat <- kalman_smoothed_means(ss, kc, cells$values)
    L_sigma <- if (shocks) t(chol(sigma[i, , ])) else NULL
    for (k in seq_len(n_paths)) {
      S_draw <- S_hat
      if (shocks) {
        # Step 1: an unconditional path S^+ from the known state at T.
        S_plus <- matrix(NA_real_, horizon, length(ss$S_T))
        state <- ss$S_T
        for (t in seq_len(horizon)) {
          w <- rep(0, length(state))
          w[seq_len(n_var)] <- as.vector(L_sigma %*% rnorm(n_var))
          state <- as.vector(ss$G %*% state) + w
          S_plus[t, ] <- state
        }
        values_plus <- lapply(seq_len(horizon), function(t) S_plus[t, cells$observed[[t]]])
        # Steps 2-3: S_tilde = S_hat - S_hat^+ + S^+.
        S_draw <- S_hat - kalman_smoothed_means(ss, kc, values_plus) + S_plus
      }
      row <- row + 1
      out[row, , ] <- S_draw[, seq_len(n_var), drop = FALSE]
    }
  }
  out
}

#' Conditional means and per-quarter covariances of one draw (Kalman method).
#'
#' Implements docs/methodology.md (M25). Used by checks.
#'
#' @param B,Sigma One draw (K x n, n x n).
#' @param history,lags,horizon,conditions As in forecast_conditional_kalman().
#' @return List: mean (H x n), var (list of H matrices n x n: the covariance
#'   of y_t given the fixed cells, from the top-left block of the smoothed
#'   state covariance).
conditional_moments_kalman <- function(B, Sigma, history, lags, horizon, conditions) {
  check_conditions(conditions, history, horizon)
  n_var <- ncol(B)
  cells <- kalman_observed_cells(conditions)
  ss <- var_state_space(B, Sigma, history, lags)
  kc <- kalman_covariances(ss, cells$observed)
  S_hat <- kalman_smoothed_means(ss, kc, cells$values)
  V <- kalman_smoothed_variances(kc)
  list(mean = matrix(S_hat[, seq_len(n_var)], horizon, n_var,
                     dimnames = list(paste0("h", seq_len(horizon)), colnames(history))),
       var = lapply(V, function(v) v[seq_len(n_var), seq_len(n_var)]))
}
