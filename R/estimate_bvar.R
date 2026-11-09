# estimate_bvar.R: estimation of the Bayesian VAR (estimation stage): the
# hyperparameter posterior, its mode, the Metropolis sampler, and the draws of
# the coefficients and the innovation covariance.
#
# Why it exists: the model's hyperparameters (lambda, mu, psi) are estimated
# from their posterior, and the coefficients B and covariance Sigma are drawn
# given each hyperparameter draw, following Giannone, Lenza and Primiceri
# (2015), Appendix B (the method Crump et al. 2021, Section 2.2, use):
#   1. find the mode of the hyperparameter posterior;
#   2. propose new hyperparameters from a normal distribution centred on the
#      current draw, with covariance c * W, W = inverse of minus the Hessian of
#      the log posterior at the mode, c tuned for an acceptance rate of about
#      20-30 percent;
#   3. accept or reject (Metropolis);
#   4. for each retained hyperparameter draw, draw Sigma and B from their
#      normal-inverse-Wishart posterior.
# The hyperparameter posterior is the one computed in R/estimate_ml.R
# (log_posterior_direct), where the sum-of-coefficients dummy observations are
# treated with the textbook formula. The hyperparameters are sampled on the log
# scale, theta = log(h); the target density of theta includes the Jacobian of
# that change of variables, so the posterior of h is unchanged.
#
# Notation: h = (lambda, soc, psi1, ..., psin) hyperparameters (soc = mu);
# theta = log(h); n variables, p lags, K = 1 + n p regressors.
#
# The formulas and settings are set out in docs/methodology.md (sections 5-7);
# code comments cite its equation labels (M1), (M2), ...
#
# Inputs:  the estimation matrix (est_data_to_matrix()), the prior and model
#          configuration, the series configuration.
# Outputs: a list with the hyperparameter draws, coefficient draws
#          (beta: draws x K x n), covariance draws (sigma: draws x n x n),
#          the mode search results, the proposal and sampler diagnostics.


# ---- Setup ----

#' Collect everything the estimation needs in one object.
#'
#' @param y Numeric estimation matrix T x n (est_data_to_matrix()).
#' @param series Series configuration (gives the prior mean of each variable).
#' @param prior_cfg The `prior` section of the run configuration.
#' @param lags Number of lags p.
#' @param keep NULL (every quarter's equation is used), or a logical vector
#'   with one entry per row of y: FALSE for the quarters whose equation is
#'   left out (COVID treatment `exclude_window`; var_data_matrices()).
#' @return List: y, lags, prior_cfg, b_own, ybar (centring of the
#'   sum-of-coefficients dummies), start, lower, upper (named vectors on the
#'   original scale, in the order lambda, soc, psi1..psin), psi_table; and
#'   keep, only when quarters are excluded.
estimation_setup <- function(y, series, prior_cfg, lags, keep = NULL) {
  # Step 1: prior mean of the own first lag and centring of the dummies.
  # The centring uses the first p observations, before any excluded quarter.
  b_own <- prior_mean_vector(series)
  ybar  <- soc_centring_means(y, lags, prior_cfg$soc_ybar)

  # Step 2: starting values and bounds of psi (as for the prior object), on
  # the equations used in the estimation.
  psi_table <- psi_start_table(y, lags, prior_cfg$psi$start, prior_cfg$psi$bound_factor, keep)
  psi_names <- paste0("psi", seq_len(ncol(y)))

  # Step 3: starting values and bounds of all hyperparameters.
  start <- c(lambda = prior_cfg$lambda$mode, soc = prior_cfg$soc$mode,
             setNames(psi_table$start, psi_names))
  lower <- c(lambda = prior_cfg$lambda$min, soc = prior_cfg$soc$min,
             setNames(psi_table$min, psi_names))
  upper <- c(lambda = prior_cfg$lambda$max, soc = prior_cfg$soc$max,
             setNames(psi_table$max, psi_names))

  setup <- list(y = y, lags = lags, prior_cfg = prior_cfg, b_own = b_own, ybar = ybar,
                start = start, lower = lower, upper = upper, psi_table = psi_table,
                all_names = names(start), fixed = NULL)
  # Stored only when quarters are excluded, so that the setup of a run
  # without exclusions is unchanged.
  if (!is.null(keep)) {
    if (length(keep) != nrow(y) || anyNA(keep)) stop("keep must have one TRUE/FALSE per row of y")
    setup$keep <- as.logical(keep)
  }
  setup
}

#' Hold some hyperparameters fixed; the others are estimated.
#'
#' Used for checks (for example, estimating only lambda and mu with psi held
#' at given values). The estimated hyperparameters keep their starting values
#' and bounds; the fixed ones are stored in `fixed` and added back whenever
#' the full vector is needed (complete_hyper()).
#'
#' @param setup Output of estimation_setup().
#' @param fixed Named numeric vector of the hyperparameters to fix and their values.
#' @return The setup with start, lower, upper restricted to the estimated
#'   hyperparameters and `fixed` set.
fix_hyperparameters <- function(setup, fixed) {
  if (!all(names(fixed) %in% setup$all_names)) stop("unknown hyperparameter in `fixed`")
  free <- setdiff(setup$all_names, names(fixed))
  setup$start <- setup$start[free]
  setup$lower <- setup$lower[free]
  setup$upper <- setup$upper[free]
  setup$fixed <- fixed
  setup
}

#' Full hyperparameter vector from the estimated ones plus the fixed ones.
#'
#' @param h_free Named vector of the estimated hyperparameters.
#' @param setup Output of estimation_setup() (possibly with fixed values).
#' @return Named vector in the order lambda, soc, psi1..psin.
complete_hyper <- function(h_free, setup) {
  c(h_free, setup$fixed)[setup$all_names]
}


# ---- Target density ----

#' Log posterior of the hyperparameters on the original scale.
#'
#' @param h Named vector (lambda, soc, psi1..psin).
#' @param setup Output of estimation_setup().
#' @return Scalar: log p(h | Y) up to a constant; -Inf outside the bounds
#'   (the hyperpriors are truncated to the bounds).
log_post_hyper <- function(h, setup) {
  if (any(h < setup$lower | h > setup$upper)) return(-Inf)
  log_post_inside(complete_hyper(h, setup), setup)
}

#' Log posterior of the hyperparameters, without the bounds check.
#'
#' @param h Full named vector (lambda, soc, psi1..psin), assumed inside the bounds.
#' @param setup Output of estimation_setup().
#' @return Scalar: log p(h | Y) up to a constant.
log_post_inside <- function(h, setup) {
  psi <- h[grep("^psi[0-9]+$", names(h))]
  log_posterior_direct(setup$y, setup$lags, h[["lambda"]], h[["soc"]], psi, setup$b_own,
                       setup$prior_cfg, setup$ybar, dummy_beta_at = "posterior",
                       keep = setup$keep)$log_posterior
}

#' Log posterior of theta = log(h), the scale on which the sampler moves.
#'
#' Change of variables: p(theta) = p(h) * |dh / dtheta| = p(h) * prod(h), so
#' log p(theta) = log p(h) + sum(theta).
#' The bounds are checked on the log scale: exp(log(x)) can differ from x in
#' the last digit, so a point exactly on a bound in theta would otherwise be
#' wrongly rejected on the original scale.
#'
#' Implements docs/methodology.md (M19).
#'
#' @param theta Named vector, log of the estimated hyperparameters (the
#'   fixed ones, if any, are added back by complete_hyper()).
#' @param setup Output of estimation_setup().
#' @return Scalar; -Inf outside the bounds. The Jacobian covers only the
#'   estimated hyperparameters.
log_post_theta <- function(theta, setup) {
  if (any(theta < log(setup$lower) | theta > log(setup$upper))) return(-Inf)
  log_post_inside(complete_hyper(exp(theta), setup), setup) + sum(theta)
}


# ---- Step 1: mode ----

#' Objective of the mode search: log p(theta), with theta clipped to the log bounds.
#'
#' L-BFGS-B keeps its iterates inside the bounds, but at a bound its
#' finite-difference gradient can evaluate a point that lies outside by one
#' rounding step (seen at full size: 4.4e-16 above log(50) for mu). There
#' log_post_theta() returns -Inf and optim() stops with an error. Clipping
#' changes nothing inside the bounds; it only maps such a point onto the bound.
#' The Metropolis sampler uses log_post_theta() itself, so proposals outside
#' the bounds are still rejected.
#'
#' Implements docs/methodology.md section 7, step 1.
#'
#' @param theta Named vector on the log scale.
#' @param setup Output of estimation_setup().
#' @return Scalar log p(theta) at the clipped point.
mode_objective <- function(theta, setup) {
  clipped <- pmin(pmax(theta, log(setup$lower)), log(setup$upper))
  log_post_theta(clipped, setup)
}

#' Find the mode of the theta posterior from several starting points.
#'
#' Each start is optimized with L-BFGS-B within the (log) bounds. Reporting
#' every start shows whether they reach the same mode; the best one is used.
#'
#' Implements docs/methodology.md section 7, step 1.
#'
#' @param setup Output of estimation_setup().
#' @param start_factors Numeric vector: each start multiplies the configured
#'   starting values (all hyperparameters) by one of these factors.
#' @param maxit Iteration limit of the optimizer.
#' @return List: table (one row per start: factor, convergence, message, log
#'   posterior, and the mode on the original scale), theta_mode (best start,
#'   named, log scale), best (row index).
find_hyper_mode <- function(setup, start_factors = c(1, 0.5, 2), maxit = 2000) {
  lower_t <- log(setup$lower)
  upper_t <- log(setup$upper)
  rows <- list()
  thetas <- list()
  for (i in seq_along(start_factors)) {
    # Step 1: starting point, kept inside the bounds.
    h0 <- pmin(pmax(setup$start * start_factors[i], setup$lower * 1.01), setup$upper / 1.01)
    theta0 <- log(h0)

    # Step 2: maximize log p(theta).
    opt <- optim(par = theta0, fn = mode_objective, setup = setup, method = "L-BFGS-B",
                 lower = lower_t, upper = upper_t,
                 control = list(fnscale = -1, maxit = maxit))
    names(opt$par) <- names(setup$start)
    thetas[[i]] <- opt$par

    # Step 3: record the result on the original scale.
    h <- exp(opt$par)
    rows[[i]] <- data.frame(start_factor = start_factors[i], convergence = opt$convergence,
                            message = opt$message, log_post_theta = opt$value,
                            t(h), check.names = FALSE)
  }
  table <- do.call(rbind, rows)
  best  <- which.max(table$log_post_theta)
  list(table = table, theta_mode = thetas[[best]], best = best)
}


# ---- Step 2: proposal ----

#' Proposal covariance W: inverse of minus the Hessian at the mode.
#'
#' The Hessian of log p(theta) is computed numerically (stats::optimHess).
#' At a maximum it is negative definite, so -H is positive definite and its
#' inverse is the covariance of the normal approximation to the posterior.
#'
#' Implements docs/methodology.md section 7, step 2.
#'
#' @param theta_mode Named mode on the log scale.
#' @param setup Output of estimation_setup().
#' @return List: W (covariance matrix), hessian, eigenvalues of -H (all must
#'   be positive), and the implied approximate posterior s.d. of each theta.
hyper_proposal <- function(theta_mode, setup) {
  H <- optimHess(theta_mode, log_post_theta, setup = setup)
  neg_H_eigen <- eigen(-H, symmetric = TRUE, only.values = TRUE)$values
  if (any(neg_H_eigen <= 0)) {
    stop("the Hessian at the mode is not negative definite: the mode may be at a bound ",
         "or not a maximum")
  }
  W <- solve(-H)
  W <- (W + t(W)) / 2                       # remove rounding asymmetry
  dimnames(W) <- list(names(theta_mode), names(theta_mode))
  list(W = W, hessian = H, neg_hessian_eigenvalues = neg_H_eigen,
       approx_sd = sqrt(diag(W)))
}


# ---- Step 3: Metropolis sampler ----

#' Random-walk Metropolis sampler for any log density.
#'
#' Proposal: theta* = theta + L z, z ~ N(0, I), L L' = c W. Accept with
#' probability min(1, exp(log_target(theta*) - log_target(theta))).
#' The scale c starts at 2.38^2 / k (the standard choice for a normal target
#' in k dimensions, acceptance about 23 percent) and is adjusted only during
#' the first `adapt_share` of the burn-in, in windows of `window` iterations: times
#' 0.8 if the window's acceptance is below `target[1]`, times 1.25 if above
#' `target[2]`. After that c is fixed, so the retained draws come from a
#' sampler with a fixed proposal (adapting it further would invalidate them).
#'
#' Implements docs/methodology.md section 7, step 4.
#'
#' @param log_target function(theta) returning the log density (up to a constant).
#' @param theta_start Named starting vector.
#' @param W Proposal covariance (k x k).
#' @param n_draw Total iterations, including burn-in.
#' @param n_burn Burn-in iterations (discarded).
#' @param n_thin Keep every n_thin-th iteration after the burn-in.
#' @param seed Seed of the random number generator.
#' @param target Acceptance band used for tuning during the burn-in.
#' @param adapt_share Share of the burn-in during which c is tuned.
#' @param window Length of a tuning window (iterations).
#' @return List: theta (retained draws, matrix), log_target (per retained
#'   draw), acceptance (after burn-in), scale (final c), scale_path (c and
#'   acceptance per tuning window), burn_theta and burn_log_target (every
#'   n_thin-th burn-in iteration, discarded from the posterior but kept so
#'   that trace plots show how the chain moved away from its start), settings.
run_metropolis <- function(log_target, theta_start, W, n_draw, n_burn, n_thin = 1, seed,
                           target = c(0.20, 0.35), adapt_share = 0.5, window = 100) {
  if (n_burn >= n_draw) stop("n_burn must be smaller than n_draw")
  set.seed(seed)
  k <- length(theta_start)
  scale <- 2.38^2 / k
  L <- t(chol(scale * W))
  n_adapt <- floor(n_burn * adapt_share)

  # Storage for the retained draws.
  n_keep <- floor((n_draw - n_burn) / n_thin)
  theta_store <- matrix(NA_real_, n_keep, k, dimnames = list(NULL, names(theta_start)))
  lp_store    <- rep(NA_real_, n_keep)
  scale_path  <- list()
  n_burn_keep <- floor(n_burn / n_thin)
  burn_theta  <- matrix(NA_real_, n_burn_keep, k, dimnames = list(NULL, names(theta_start)))
  burn_lp     <- rep(NA_real_, n_burn_keep)

  # Current state.
  theta <- theta_start
  lp    <- log_target(theta)
  if (!is.finite(lp)) stop("the starting point has zero posterior density")
  accepted_window <- 0
  accepted_after_burn <- 0
  keep <- 0

  for (iter in seq_len(n_draw)) {
    # Step 1: propose.
    theta_prop <- theta + as.vector(L %*% rnorm(k))
    names(theta_prop) <- names(theta)
    lp_prop <- log_target(theta_prop)

    # Step 2: accept or reject.
    accept <- is.finite(lp_prop) && log(runif(1)) < lp_prop - lp
    if (accept) {
      theta <- theta_prop
      lp    <- lp_prop
    }

    # Step 3: tune the scale during the first part of the burn-in.
    if (iter <= n_adapt) {
      accepted_window <- accepted_window + accept
      if (iter %% window == 0) {
        rate <- accepted_window / window
        scale_path[[length(scale_path) + 1]] <- c(iteration = iter, scale = scale, acceptance = rate)
        if (rate < target[1]) scale <- scale * 0.8
        if (rate > target[2]) scale <- scale * 1.25
        L <- t(chol(scale * W))
        accepted_window <- 0
      }
    }

    # Step 4: during the burn-in, record every n_thin-th iteration for the
    # trace plots; after it, count acceptances and store every n_thin-th draw.
    if (iter <= n_burn && iter %% n_thin == 0) {
      burn_theta[iter / n_thin, ] <- theta
      burn_lp[iter / n_thin] <- lp
    }
    if (iter > n_burn) {
      accepted_after_burn <- accepted_after_burn + accept
      if ((iter - n_burn) %% n_thin == 0 && keep < n_keep) {
        keep <- keep + 1
        theta_store[keep, ] <- theta
        lp_store[keep] <- lp
      }
    }
  }

  list(theta = theta_store, log_target = lp_store,
       acceptance = accepted_after_burn / (n_draw - n_burn),
       scale = scale,
       scale_path = if (length(scale_path)) do.call(rbind, scale_path) else NULL,
       burn_theta = burn_theta, burn_log_target = burn_lp,
       settings = list(n_draw = n_draw, n_burn = n_burn, n_thin = n_thin, seed = seed,
                       target = target, adapt_share = adapt_share, window = window,
                       n_adapt = n_adapt))
}


# ---- Step 4: coefficients and covariance given the hyperparameters ----

#' Posterior quantities of B and Sigma given the hyperparameters.
#'
#' The data are augmented with the sum-of-coefficients dummy observations;
#' then (Giannone, Lenza and Primiceri 2015, Appendix A):
#'   B_hat = (X'X + Omega^-1)^-1 (X'Y + Omega^-1 b)
#'   S     = (Y - X B_hat)'(Y - X B_hat) + (B_hat - b)' Omega^-1 (B_hat - b)
#'   Sigma | h, Y     ~ IW(Psi + S, N + d)        (N includes the dummy rows)
#'   vec(B) | Sigma   ~ N(vec(B_hat), Sigma (x) (X'X + Omega^-1)^-1)
#' (X'X + Omega^-1)^-1 is computed in the rescaled form
#' Omega^1/2 (I + Omega^1/2 X'X Omega^1/2)^-1 Omega^1/2 (see log_ml_niw()).
#'
#' Implements docs/methodology.md (M12), (M13), (M16), (M17).
#'
#' @param h Named hyperparameters (lambda, soc, psi1..psin).
#' @param setup Output of estimation_setup().
#' @return List: B_hat (K x n), scale_matrix (Psi + S), dof (N + d),
#'   P (K x K, with P P' = (X'X + Omega^-1)^-1).
niw_posterior <- function(h, setup) {
  y <- setup$y; lags <- setup$lags; pc <- setup$prior_cfg
  n_var <- ncol(y)
  psi <- h[grep("^psi[0-9]+$", names(h))]

  # Step 1: data with the dummy observations on top (excluded quarters'
  # equations left out, if any).
  dm  <- var_data_matrices(y, lags, setup$keep)
  dum <- soc_dummy_obs(setup$ybar, h[["soc"]], lags)
  Y   <- rbind(dum$Y, dm$Y)
  X   <- rbind(dum$X, dm$X)

  # Step 2: prior mean and Omega.
  b     <- prior_mean_matrix(setup$b_own, lags)
  omega <- minnesota_omega(h[["lambda"]], pc$alpha, psi, lags, pc$constant_var)
  omega_half <- diag(sqrt(omega), nrow = length(omega))

  # Step 3: (X'X + Omega^-1)^-1 = P P', P = Omega^1/2 V D^-1/2, from the
  # eigendecomposition V D V' of I + Omega^1/2 X'X Omega^1/2.
  G   <- omega_half %*% crossprod(X) %*% omega_half
  eig <- eigen((G + t(G)) / 2, symmetric = TRUE)
  ev  <- pmax(eig$values, 0) + 1
  P   <- omega_half %*% eig$vectors %*% diag(1 / sqrt(ev), nrow = length(ev))

  # Step 4: posterior mean of B and the residual cross-product S.
  B_hat <- b + P %*% t(P) %*% crossprod(X, Y - X %*% b)
  E     <- Y - X %*% B_hat
  omega_inv <- diag(1 / omega, nrow = length(omega))
  S     <- crossprod(E) + t(B_hat - b) %*% omega_inv %*% (B_hat - b)

  list(B_hat = B_hat, scale_matrix = diag(psi, nrow = n_var) + S,
       dof = nrow(Y) + n_var + 2, P = P)
}

#' Draw Sigma and B given the hyperparameters.
#'
#' Sigma ~ IW(scale_matrix, dof): draw W ~ Wishart(dof, scale_matrix^-1) and
#' set Sigma = W^-1. B = B_hat + P Z L', with Z a K x n matrix of independent
#' standard normals and L L' = Sigma, so that vec(B) has covariance
#' Sigma (x) P P'.
#'
#' Implements docs/methodology.md (M16), (M17).
#'
#' @param post Output of niw_posterior().
#' @param n_draws Number of draws.
#' @return List: beta (n_draws x K x n), sigma (n_draws x n x n).
draw_coefficients <- function(post, n_draws) {
  K <- nrow(post$B_hat); n_var <- ncol(post$B_hat)
  beta  <- array(NA_real_, c(n_draws, K, n_var))
  sigma <- array(NA_real_, c(n_draws, n_var, n_var))
  scale_inv <- solve(post$scale_matrix)
  scale_inv <- (scale_inv + t(scale_inv)) / 2
  for (i in seq_len(n_draws)) {
    # Step 1: Sigma from the inverse Wishart.
    W_draw <- stats::rWishart(1, df = post$dof, Sigma = scale_inv)[, , 1]
    Sigma  <- solve(W_draw)
    Sigma  <- (Sigma + t(Sigma)) / 2
    # Step 2: B given Sigma.
    L <- t(chol(Sigma))
    Z <- matrix(rnorm(K * n_var), K, n_var)
    beta[i, , ]  <- post$B_hat + post$P %*% Z %*% t(L)
    sigma[i, , ] <- Sigma
  }
  list(beta = beta, sigma = sigma)
}

#' One draw of Sigma and B for each row of a matrix of hyperparameter draws.
#'
#' Step 4 of the algorithm: for hyperparameter draw h_m, draw
#' (Sigma_m, B_m) from p(Sigma, B | h_m, Y) (niw_posterior(), then
#' draw_coefficients()). Each pair (h_m, Sigma_m, B_m) is then a draw from the
#' joint posterior p(h, Sigma, B | Y).
#'
#' Implements docs/methodology.md section 7, step 6.
#'
#' @param hyper Matrix of hyperparameter draws on the original scale, one row
#'   per draw, columns named lambda, soc, psi1..psin (all hyperparameters,
#'   complete_hyper()).
#' @param setup Output of estimation_setup().
#' @param seed Seed, set once before the first draw.
#' @return List: beta (draws x K x n), sigma (draws x n x n).
draw_posterior_coefficients <- function(hyper, setup, seed) {
  n_keep <- nrow(hyper)
  K <- 1 + ncol(setup$y) * setup$lags; n_var <- ncol(setup$y)
  beta  <- array(NA_real_, c(n_keep, K, n_var))
  sigma <- array(NA_real_, c(n_keep, n_var, n_var))
  set.seed(seed)
  for (i in seq_len(n_keep)) {
    post <- niw_posterior(hyper[i, ], setup)
    d <- draw_coefficients(post, 1)
    beta[i, , ]  <- d$beta[1, , ]
    sigma[i, , ] <- d$sigma[1, , ]
  }
  list(beta = beta, sigma = sigma)
}


# ---- All steps ----

#' Estimate the BVAR: mode, proposal, Metropolis, coefficient draws.
#'
#' @param setup Output of estimation_setup().
#' @param n_draw,n_burn,n_thin Metropolis settings (iterations, burn-in, thinning).
#' @param seed Seed (the Metropolis step and the coefficient draws use it in turn).
#' @param start_factors Starting points of the mode search (see find_hyper_mode()).
#' @return List: hyper (retained draws on the original scale, matrix, all
#'   hyperparameters; fixed ones constant), estimated (names of the estimated ones),
#'   log_post (per draw), beta (draws x K x n), sigma (draws x n x n),
#'   mode (find_hyper_mode() result), proposal (hyper_proposal() result),
#'   sampler (acceptance, scale, scale_path, settings), variables, timing.
estimate_bvar <- function(setup, n_draw, n_burn, n_thin = 1, seed,
                          start_factors = c(1, 0.5, 2)) {
  timing <- list()

  # Step 1: mode of the hyperparameter posterior.
  t0 <- Sys.time()
  mode <- find_hyper_mode(setup, start_factors)
  timing$mode <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Step 2: proposal covariance from the Hessian at the mode.
  t0 <- Sys.time()
  proposal <- hyper_proposal(mode$theta_mode, setup)
  timing$hessian <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Step 3: Metropolis on theta = log(h), started at the mode.
  t0 <- Sys.time()
  target <- function(theta) log_post_theta(theta, setup)
  mh <- run_metropolis(target, mode$theta_mode, proposal$W, n_draw, n_burn, n_thin, seed)
  timing$metropolis <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Step 4: Sigma and B for every retained hyperparameter draw (fixed
  # hyperparameters, if any, are added back as constant columns).
  t0 <- Sys.time()
  hyper <- t(apply(exp(mh$theta), 1, complete_hyper, setup = setup))
  coef  <- draw_posterior_coefficients(hyper, setup, seed + 1)
  timing$coefficients <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  list(hyper = hyper, estimated = names(setup$start),
       log_post = mh$log_target - rowSums(mh$theta), beta = coef$beta, sigma = coef$sigma,
       mode = mode, proposal = proposal,
       sampler = mh[c("acceptance", "scale", "scale_path", "settings")],
       variables = colnames(setup$y), timing = timing)
}
