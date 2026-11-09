# estimate_run.R: the estimation stage as it runs in the pipeline (stage 3):
# MCMC settings, several Metropolis chains, the pooled posterior draws, the
# convergence summary, the run identifier and the run's output files.
#
# Why it exists: R/estimate_bvar.R implements the estimator (Giannone, Lenza
# and Primiceri 2015, Appendix B) step by step. The pipeline (_targets.R)
# runs those steps as separate targets, so that an expensive step is redone
# only when its own inputs change: the mode search and the Hessian do not
# depend on the chain lengths, so changing the number of draws reruns only
# the chains. This file holds the functions those targets call:
#   mcmc_settings()      the MCMC settings of the chosen profile, checked;
#   covid_for_method()   the COVID treatment (R/covid.R) for one method;
#   timed_call()         run one estimation step and record its run time;
#   chain_start()        the starting point of one chain (the mode, or an
#                        overdispersed point around it);
#   run_chains()         several Metropolis chains, in parallel;
#   pool_and_draw()      pool the chains, keep an evenly spaced subset of the
#                        hyperparameter draws, draw Sigma and B for each;
#   summarize_mcmc()     acceptance, effective sample size, R-hat, Geweke,
#                        draws near a bound, companion-matrix eigenvalues;
#   check_convergence()  the run's diagnostics against the convergence targets;
#   derived_quantities() quantities the engine uses, per kept draw;
#   make_run_id()        the run identifier;
#   write_run_files()    outputs/<run_id>/: draws, estimation details,
#                        summary, meta.json.
#
# Two models are estimated side by side, one per COVID method (config
# covid$methods): "none" (parameters estimated through the latest quarter)
# and "precovid_params" (parameters estimated through covid$cutoff). Each is a
# separate run with its own run_id and output folder.
#
# Notation as in R/estimate_bvar.R: h = (lambda, soc, psi1..psin)
# hyperparameters, theta = log(h); n variables, p lags, K = 1 + n p.
#
# The formulas and settings are set out in docs/methodology.md (sections 7-8);
# code comments cite its equation labels (M1), (M2), ...
#
# Inputs:  the run configuration, the series configuration, model_data from
#          the data stage.
# Outputs: R objects for the pipeline's targets, and the files under
#          <output_root>/<run_id>/ (write_run_files()).


# ---- Settings ----

#' MCMC settings of the configured profile, checked before anything runs.
#'
#' The run configuration has two profiles: "smoke" (a short run that checks
#' the pipeline end to end; its draws are not for use) and "full" (the
#' production run); a third, "pilot", is a medium-length run used to measure
#' how fast the chains mix before the full run is sized. Seeds: chain k uses
#' seed + k; its starting point (chain_start()) uses seed + 2000 + k; the
#' coefficient draws use seed + 1000. Every chain therefore has its own
#' random-number streams, and the results do not depend on how the chains are
#' spread over processor cores.
#'
#' @param cfg Run configuration (list from config/run.yml).
#' @return List: profile, n_draw, n_burn, n_thin, n_chains, n_engine_draws,
#'   n_keep_per_chain, seed, chain_seeds, start_seeds, coefficient_seed,
#'   chain_start (type, scale), tuning (target, adapt_share, window). Stops if
#'   a setting of the profile is missing or invalid.
mcmc_settings <- function(cfg) {
  m <- cfg$mcmc
  profile <- m$profile
  if (is.null(profile) || !profile %in% names(m$profiles)) {
    stop("mcmc$profile must be one of: ", paste(names(m$profiles), collapse = ", "))
  }
  p <- m$profiles[[profile]]

  # Step 1: every count must be set (the full profile is null until the
  # draw counts are chosen from full-size diagnostics).
  keys <- c("n_draw", "n_burn", "n_thin", "n_chains", "n_engine_draws")
  missing <- keys[vapply(keys, function(k) is.null(p[[k]]), logical(1))]
  if (length(missing)) {
    stop("mcmc profile '", profile, "' has no value for: ", paste(missing, collapse = ", "),
         " (set them in config/run.yml)")
  }
  for (k in keys) {
    x <- p[[k]]
    if (!is.numeric(x) || length(x) != 1 || x < 1 || x != round(x)) {
      stop("mcmc profile '", profile, "': ", k, " must be a positive whole number")
    }
  }

  # Step 2: consistency of the counts.
  if (p$n_burn >= p$n_draw) stop("mcmc profile '", profile, "': n_burn must be smaller than n_draw")
  n_keep <- floor((p$n_draw - p$n_burn) / p$n_thin)
  if (p$n_engine_draws > n_keep * p$n_chains) {
    stop("mcmc profile '", profile, "': n_engine_draws (", p$n_engine_draws,
         ") exceeds the retained draws of all chains (", n_keep * p$n_chains, ")")
  }

  # Step 3: tuning of the proposal scale (run_metropolis()).
  tu <- m$tuning
  if (length(tu$target) != 2 || tu$target[1] >= tu$target[2]) {
    stop("mcmc$tuning$target must be two increasing acceptance rates")
  }

  # Step 4: where the chains start (chain_start()).
  cs <- m$chain_start
  if (is.null(cs$type) || !cs$type %in% c("mode", "overdispersed")) {
    stop("mcmc$chain_start$type must be 'mode' or 'overdispersed'")
  }
  if (cs$type == "overdispersed" && (!is.numeric(cs$scale) || length(cs$scale) != 1 || cs$scale <= 0)) {
    stop("mcmc$chain_start$scale must be a positive number")
  }
  chain_start_cfg <- list(type = cs$type, scale = if (cs$type == "mode") 0 else cs$scale)

  list(profile = profile, n_draw = p$n_draw, n_burn = p$n_burn, n_thin = p$n_thin,
       n_chains = p$n_chains, n_engine_draws = p$n_engine_draws, n_keep_per_chain = n_keep,
       seed = m$seed, chain_seeds = m$seed + seq_len(p$n_chains),
       start_seeds = m$seed + 2000 + seq_len(p$n_chains), coefficient_seed = m$seed + 1000,
       chain_start = chain_start_cfg,
       tuning = list(target = unlist(tu$target), adapt_share = tu$adapt_share, window = tu$window))
}


# ---- Pipeline steps ----

#' COVID treatment for one method.
#'
#' @param model_data data.frame from the data stage.
#' @param cfg Run configuration.
#' @param method One of the methods in cfg$covid$methods.
#' @return Output of covid_treatment() (R/covid.R) with cfg$covid$method = method.
covid_for_method <- function(model_data, cfg, method) {
  cfg$covid$method <- method
  covid_treatment(model_data, cfg)
}

#' Run one step and record how long it took.
#'
#' @param f Function returning a list (for example find_hyper_mode).
#' @param ... Arguments passed to f.
#' @return f(...) with an element `seconds` added (wall-clock time).
timed_call <- function(f, ...) {
  t0 <- proc.time()[["elapsed"]]
  out <- f(...)
  out$seconds <- proc.time()[["elapsed"]] - t0
  out
}

#' Starting point of one chain.
#'
#' Why overdispersed starts: R-hat (Gelman and Rubin 1992) compares the
#' spread between chains with the spread within them. Chains that all start
#' at the same point can agree with each other before any of them has
#' explored the posterior, so R-hat close to 1 would prove little. Started
#' from points more spread out than the posterior, chains that agree must
#' each have forgotten their start.
#'
#' The overdispersed start is theta_mode + scale * L z, z ~ N(0, I),
#' L L' = W (the normal approximation to the posterior at the mode,
#' hyper_proposal()), so with scale 2 each start lies about 2 approximate
#' posterior standard deviations from the mode in every direction. A
#' coordinate that would fall outside the search bounds is moved to 1
#' percent inside the bound (as the starts of find_hyper_mode()).
#'
#' Implements docs/methodology.md (M20).
#'
#' @param theta_mode Named mode on the log scale.
#' @param W Proposal covariance (hyper_proposal()).
#' @param setup Output of estimation_setup() (bounds).
#' @param type "mode" (start at the mode) or "overdispersed".
#' @param scale Multiple of the approximate posterior s.d. (overdispersed only).
#' @param seed Seed of the start's random numbers (separate from the chain's).
#' @return List: theta (named start on the log scale), clipped (names of the
#'   coordinates moved inside a bound; empty if none).
chain_start <- function(theta_mode, W, setup, type, scale, seed) {
  if (type == "mode") return(list(theta = theta_mode, clipped = character(0)))

  # Step 1: a point around the mode, spread as scale^2 W.
  set.seed(seed)
  z <- rnorm(length(theta_mode))
  theta <- theta_mode + scale * as.vector(t(chol(W)) %*% z)
  names(theta) <- names(theta_mode)

  # Step 2: keep it inside the bounds (log scale, 1 percent margin).
  lo <- log(setup$lower * 1.01)
  hi <- log(setup$upper / 1.01)
  clipped <- names(theta)[theta < lo | theta > hi]
  theta <- pmin(pmax(theta, lo), hi)
  list(theta = theta, clipped = clipped)
}

#' One Metropolis chain on theta = log(h).
#'
#' @param k Index of the chain (its seeds are settings$chain_seeds[k] and
#'   settings$start_seeds[k]).
#' @param setup Output of estimation_setup().
#' @param theta_mode Named mode on the log scale (find_hyper_mode()).
#' @param W Proposal covariance (hyper_proposal()).
#' @param settings Output of mcmc_settings().
#' @return Output of run_metropolis(), with `seconds`, `start` (the starting
#'   point, log scale) and `start_clipped` (chain_start()) added.
run_one_chain <- function(k, setup, theta_mode, W, settings) {
  cs <- settings$chain_start
  start <- chain_start(theta_mode, W, setup, cs$type, cs$scale, settings$start_seeds[k])
  target <- function(theta) log_post_theta(theta, setup)
  tu <- settings$tuning
  out <- timed_call(run_metropolis, target, start$theta, W, settings$n_draw, settings$n_burn,
                    settings$n_thin, settings$chain_seeds[k], target = tu$target,
                    adapt_share = tu$adapt_share, window = tu$window)
  out$start <- start$theta
  out$start_clipped <- start$clipped
  out
}

#' Run several Metropolis chains, in parallel when more than one core is given.
#'
#' Every chain uses the same proposal c W; the chains differ in their seed
#' and, with overdispersed starts, in their starting point (mcmc_settings(),
#' chain_start()). Comparing them (R-hat, in summarize_mcmc()) shows whether
#' each one has explored the same posterior.
#' Parallel runs use a socket cluster (base R `parallel`, which works on
#' Windows): each worker is a separate R process that loads the project's
#' functions from `r_dir`. Because each chain sets its own seed, the results
#' are identical to running the chains one after the other.
#'
#' @param setup Output of estimation_setup().
#' @param theta_mode Named mode on the log scale.
#' @param W Proposal covariance.
#' @param settings Output of mcmc_settings().
#' @param n_cores Number of R processes to use (1 = run in this process).
#' @param r_dir Folder with the project's R files (sourced by each worker).
#' @return List of run_metropolis() outputs, one per chain, each with
#'   `chain` (index) and `seconds`.
run_chains <- function(setup, theta_mode, W, settings, n_cores = 1, r_dir = "R") {
  idx <- seq_along(settings$chain_seeds)
  n_workers <- min(n_cores, length(idx))

  if (n_workers <= 1) {
    chains <- lapply(idx, run_one_chain, setup = setup, theta_mode = theta_mode, W = W,
                     settings = settings)
  } else {
    cl <- parallel::makePSOCKcluster(n_workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    # Each worker uses this project's package library and functions.
    r_files <- normalizePath(list.files(r_dir, pattern = "\\.R$", full.names = TRUE))
    parallel::clusterCall(cl, function(lib, files) {
      .libPaths(lib)
      for (f in files) source(f)
      NULL
    }, .libPaths(), r_files)
    chains <- parallel::clusterApply(cl, idx, run_one_chain, setup = setup,
                                     theta_mode = theta_mode, W = W, settings = settings)
  }
  for (k in seq_along(chains)) chains[[k]]$chain <- k
  chains
}

#' Pool the chains and draw Sigma and B for an evenly spaced subset.
#'
#' The retained hyperparameter draws of all chains are stacked (chain 1
#' first). n_engine_draws of them, evenly spaced over the stacked draws so
#' that every chain contributes about equally, are kept, and one (Sigma, B) is
#' drawn for each (draw_posterior_coefficients()). Keeping evenly spaced draws
#' is a further thinning of each chain; it reduces the autocorrelation
#' between kept draws.
#'
#' Implements docs/methodology.md section 7, step 6.
#'
#' @param chains Output of run_chains().
#' @param setup Output of estimation_setup().
#' @param settings Output of mcmc_settings() (n_engine_draws, coefficient_seed).
#' @return List: hyper (n_engine_draws x all hyperparameters, original
#'   scale), chain and draw (source of each kept row), log_post (log p(h | Y)
#'   per kept row, original scale), beta (draws x K x n), sigma (draws x n x n),
#'   variables, lags, seconds.
pool_and_draw <- function(chains, setup, settings) {
  t0 <- proc.time()[["elapsed"]]

  # Step 1: stack the chains.
  theta <- do.call(rbind, lapply(chains, function(ch) ch$theta))
  lp    <- unlist(lapply(chains, function(ch) ch$log_target))
  chain <- unlist(lapply(chains, function(ch) rep(ch$chain, nrow(ch$theta))))
  draw  <- unlist(lapply(chains, function(ch) seq_len(nrow(ch$theta))))

  # Step 2: evenly spaced rows.
  pick <- round(seq(1, nrow(theta), length.out = settings$n_engine_draws))

  # Step 3: hyperparameters on the original scale and the coefficient draws.
  # log p(h | Y) = log p(theta | Y) - sum(theta) (the Jacobian, log_post_theta()).
  hyper <- t(apply(exp(theta[pick, , drop = FALSE]), 1, complete_hyper, setup = setup))
  coef  <- draw_posterior_coefficients(hyper, setup, settings$coefficient_seed)

  list(hyper = hyper, chain = chain[pick], draw = draw[pick],
       log_post = lp[pick] - rowSums(theta[pick, , drop = FALSE]),
       beta = coef$beta, sigma = coef$sigma, variables = colnames(setup$y), lags = setup$lags,
       seconds = proc.time()[["elapsed"]] - t0)
}


# ---- Convergence summary ----

#' Largest modulus of the eigenvalues of the companion matrix of one draw.
#'
#' Companion form (Crump et al. 2021, eq. 1 written as a VAR(1) in the
#' stacked vector z_t = (y_t', ..., y_{t-p+1}')'):
#'   F = [ B_1 B_2 ... B_p ; I_{n(p-1)} 0 ].
#' Row 1 + (s-1) n + j of B (K x n) holds the coefficients on lag s of
#' variable j in every equation, so B_s = t(B[rows of lag s, ]).
#' A modulus above 1 means the draw implies explosive dynamics; with a prior
#' centred on unit roots, values close to 1 are expected.
#'
#' @param B Coefficient matrix of one draw, K x n (constant in row 1).
#' @param lags Number of lags p.
#' @return Scalar: the largest eigenvalue modulus.
companion_max_modulus <- function(B, lags) {
  max(Mod(eigen(companion_matrix(B, lags), only.values = TRUE)$values))
}

#' Companion matrix F of a VAR(p) from its coefficient matrix (see
#' companion_max_modulus() for the layout).
#'
#' Implements docs/methodology.md (M3).
#'
#' @param B Coefficient matrix, K x n (constant in row 1).
#' @param lags Number of lags p.
#' @return Matrix np x np.
companion_matrix <- function(B, lags) {
  n_var <- ncol(B)
  F_mat <- matrix(0, n_var * lags, n_var * lags)
  for (s in seq_len(lags)) {
    rows <- 1 + (s - 1) * n_var + seq_len(n_var)
    F_mat[seq_len(n_var), (s - 1) * n_var + seq_len(n_var)] <- t(B[rows, ])
  }
  if (lags > 1) F_mat[(n_var + 1):(n_var * lags), 1:(n_var * (lags - 1))] <- diag(n_var * (lags - 1))
  F_mat
}

#' Convergence and sanity summary of a run.
#'
#' Computed on theta = log(h), the scale on which the sampler moves:
#'   - acceptance rate and final proposal scale c of each chain;
#'   - effective sample size (ESS) of each hyperparameter, summed over chains
#'     (coda::effectiveSize): the number of independent draws the
#'     autocorrelated chain is worth;
#'   - R-hat (coda::gelman.diag, Gelman and Rubin 1992), with two or more
#'     chains: compares the variance between chains with the variance within
#'     them; values close to 1 mean the chains agree;
#'   - Geweke z (coda::geweke.diag): compares the mean of the first 10 percent
#'     of a chain with the mean of its last 50 percent; the largest |z| over
#'     chains is reported;
#'   - share of draws within 1 percent of a search bound.
#' The same chain statistics are computed for log p(theta | Y), the log
#' posterior along the chain (row "log_post"), which summarizes all
#' hyperparameters at once; the multivariate R-hat (Brooks and Gelman 1998)
#' does the same for the joint vector.
#' On the kept coefficient draws: whether every Sigma is positive definite,
#' and the largest companion-matrix eigenvalue modulus of each draw
#' (companion_max_modulus(); above 1 = explosive).
#' With `history`, also derived_quantities(); with `convergence`, the
#' check_convergence() table. Neither changes the draws.
#'
#' Implements docs/methodology.md section 8.2.
#'
#' @param chains Output of run_chains().
#' @param draws Output of pool_and_draw().
#' @param setup Output of estimation_setup().
#' @param history Optional numeric matrix (quarters x n) to forecast from
#'   (the forecast history, origin_data), for derived_quantities().
#' @param convergence Optional list: ess_min, rhat_max (config mcmc$convergence).
#' @return List: chains (data.frame per chain, including the log posterior at
#'   the start and the number of clipped start coordinates), hyper
#'   (data.frame per hyperparameter plus a "log_post" row: q05, median, q95 on
#'   the original scale, ess, ess_per_1000 (per 1,000 retained iterations of
#'   one chain), rhat, rhat_upper, geweke_max_abs_z, share_near_bound),
#'   rhat_multivariate, sigma_pd (share), max_modulus (vector, one per kept
#'   draw), max_modulus_summary, derived (or NULL), convergence (or NULL).
summarize_mcmc <- function(chains, draws, setup, history = NULL, convergence = NULL) {
  nm <- colnames(chains[[1]]$theta)
  stats <- c(nm, "log_post")

  # Step 1: per chain.
  per_chain <- data.frame(
    chain = vapply(chains, function(ch) ch$chain, numeric(1)),
    seed = vapply(chains, function(ch) ch$settings$seed, numeric(1)),
    acceptance = vapply(chains, function(ch) ch$acceptance, numeric(1)),
    scale = vapply(chains, function(ch) ch$scale, numeric(1)),
    log_post_start = vapply(chains, function(ch) {
      if (is.null(ch$start)) NA_real_ else log_post_theta(ch$start, setup)
    }, numeric(1)),
    start_clipped = vapply(chains, function(ch) length(ch$start_clipped), numeric(1)),
    seconds = vapply(chains, function(ch) ch$seconds, numeric(1)))

  # Step 2: per hyperparameter (and the log posterior), over all chains.
  # The log posterior is a column of the chain matrix, so it gets the same
  # statistics; the multivariate R-hat uses the hyperparameters only.
  with_lp <- function(ch) cbind(ch$theta, log_post = ch$log_target)
  mc <- coda::mcmc.list(lapply(chains, function(ch) coda::mcmc(with_lp(ch))))
  theta_all <- do.call(rbind, lapply(chains, with_lp))
  h_all <- cbind(exp(theta_all[, nm, drop = FALSE]), log_post = theta_all[, "log_post"])
  ess <- coda::effectiveSize(mc)[stats]
  n_iter <- length(chains) * (chains[[1]]$settings$n_draw - chains[[1]]$settings$n_burn)
  rhat_mv <- NA_real_
  if (length(chains) >= 2) {
    gd <- coda::gelman.diag(mc, autoburnin = FALSE, multivariate = FALSE)$psrf
    rhat <- gd[stats, 1]; rhat_upper <- gd[stats, 2]
    mc_h <- coda::mcmc.list(lapply(chains, function(ch) coda::mcmc(ch$theta)))
    rhat_mv <- tryCatch(coda::gelman.diag(mc_h, autoburnin = FALSE, multivariate = TRUE)$mpsrf,
                        error = function(e) NA_real_)
  } else {
    rhat <- rhat_upper <- rep(NA_real_, length(stats))
  }
  geweke <- sapply(chains, function(ch) coda::geweke.diag(coda::mcmc(with_lp(ch)))$z[stats])
  geweke <- matrix(geweke, nrow = length(stats))
  near <- vapply(nm, function(v) {
    mean(h_all[, v] <= setup$lower[[v]] * 1.01 | h_all[, v] >= setup$upper[[v]] / 1.01)
  }, numeric(1))
  hyper <- data.frame(hyperparameter = stats,
                      q05 = apply(h_all, 2, quantile, 0.05), median = apply(h_all, 2, median),
                      q95 = apply(h_all, 2, quantile, 0.95),
                      ess = unname(ess), ess_per_1000 = unname(ess) / n_iter * 1000,
                      rhat = unname(rhat), rhat_upper = unname(rhat_upper),
                      geweke_max_abs_z = apply(abs(geweke), 1, max),
                      share_near_bound = c(unname(near), NA_real_), row.names = NULL)

  # Step 3: kept coefficient draws.
  n_kept <- dim(draws$beta)[1]
  sigma_pd <- vapply(seq_len(n_kept), function(i) {
    !inherits(try(chol(draws$sigma[i, , ]), silent = TRUE), "try-error")
  }, logical(1))
  max_mod <- vapply(seq_len(n_kept), function(i) companion_max_modulus(draws$beta[i, , ], draws$lags),
                    numeric(1))

  # Step 4: derived quantities and the convergence check (optional).
  derived <- if (is.null(history)) NULL else derived_quantities(draws, history, max_mod)
  out <- list(chains = per_chain, hyper = hyper, rhat_multivariate = rhat_mv,
              sigma_pd = mean(sigma_pd), max_modulus = max_mod,
              max_modulus_summary = c(median = median(max_mod), q95 = unname(quantile(max_mod, 0.95)),
                                      max = max(max_mod), share_above_1 = mean(max_mod > 1)),
              derived = derived)
  out$convergence <- if (is.null(convergence)) NULL else check_convergence(out, convergence)
  out
}

#' Quantities the engine uses, computed per kept draw and compared across chains.
#'
#' The hyperparameter chain is what the diagnostics test, but results come
#' from the coefficient draws. Two derived quantities are tracked:
#'   - the largest companion-matrix eigenvalue modulus (stability, M3 form);
#'   - rgdp growth over the next 4 quarters on the no-shock forecast path
#'     (the forecast mean given the draw), from the last observed quarter of
#'     `history`, in log points (100 ln, about percent).
#' For each, the median and 5/95 percent quantiles by chain, and an R-hat
#' over the chains. The kept draws are evenly spaced over the pooled chains,
#' so chains contribute nearly equal numbers; each chain's draws are cut to
#' the shortest chain's count for the R-hat.
#'
#' @param draws Output of pool_and_draw().
#' @param history Numeric matrix (quarters x n), named columns, model units.
#' @param max_mod Largest companion modulus per kept draw (summarize_mcmc()).
#' @param variable Variable whose growth is tracked.
#' @param horizon Quarters ahead.
#' @return List: values (data.frame per kept draw: chain, max_modulus,
#'   growth), by_chain (data.frame: quantity, chain, q05, median, q95), rhat
#'   (named vector, one per quantity; NA with one chain), variable, horizon.
derived_quantities <- function(draws, history, max_mod, variable = "rgdp", horizon = 4) {
  # Step 1: no-shock forecast path of every kept draw (R/forecast.R;
  # docs/methodology.md (M18) without shocks), growth from the last observed quarter.
  path <- forecast_unconditional(draws$beta, draws$sigma, history, draws$lags, horizon,
                                 seed = 1, shocks = FALSE)
  growth <- unname(path[, horizon, variable] - history[nrow(history), variable])
  values <- data.frame(chain = draws$chain, max_modulus = max_mod, growth = growth)

  # Step 2: quantiles by chain and R-hat over chains.
  qty <- c("max_modulus", "growth")
  chains <- split(values, values$chain)
  rows <- list()
  for (q in qty) {
    for (k in names(chains)) {
      x <- chains[[k]][[q]]
      rows[[length(rows) + 1]] <- data.frame(quantity = q, chain = as.numeric(k),
                                             q05 = unname(quantile(x, 0.05)), median = median(x),
                                             q95 = unname(quantile(x, 0.95)))
    }
  }
  by_chain <- do.call(rbind, rows)
  rhat <- setNames(rep(NA_real_, length(qty)), qty)
  if (length(chains) >= 2) {
    n_min <- min(vapply(chains, nrow, integer(1)))
    for (q in qty) {
      mc <- coda::mcmc.list(lapply(chains, function(d) coda::mcmc(d[[q]][seq_len(n_min)])))
      rhat[[q]] <- coda::gelman.diag(mc, autoburnin = FALSE)$psrf[1, 1]
    }
  }
  list(values = values, by_chain = by_chain, rhat = rhat, variable = variable, horizon = horizon)
}

#' The run's diagnostics against the convergence targets.
#'
#' Targets (config mcmc$convergence): for every hyperparameter and for the
#' log posterior, ESS summed over chains at least ess_min, and R-hat (point
#' estimate) below rhat_max. A flag is reported, not acted on: trace plots
#' and Geweke statistics are read alongside, and a run that fails is
#' reported as failing.
#'
#' Implements docs/methodology.md section 8.2.
#'
#' @param summary List with `hyper` (summarize_mcmc()).
#' @param convergence List: ess_min, rhat_max.
#' @return data.frame, one row per check: check, worst (the worst value and
#'   where it occurs), threshold, pass.
check_convergence <- function(summary, convergence) {
  h <- summary$hyper
  i_ess <- which.min(h$ess)
  rows <- list(data.frame(check = "ESS (sum over chains), minimum",
                          worst = sprintf("%.0f (%s)", h$ess[i_ess], h$hyperparameter[i_ess]),
                          threshold = sprintf(">= %s", convergence$ess_min),
                          pass = h$ess[i_ess] >= convergence$ess_min))
  if (all(is.finite(h$rhat))) {
    i_r <- which.max(h$rhat)
    rows[[2]] <- data.frame(check = "R-hat (point), maximum",
                            worst = sprintf("%.3f (%s)", h$rhat[i_r], h$hyperparameter[i_r]),
                            threshold = sprintf("< %s", convergence$rhat_max),
                            pass = h$rhat[i_r] < convergence$rhat_max)
  } else {
    rows[[2]] <- data.frame(check = "R-hat (point), maximum", worst = "not computed (one chain)",
                            threshold = sprintf("< %s", convergence$rhat_max), pass = FALSE)
  }
  do.call(rbind, rows)
}


# ---- Run identifier and files ----

#' Run identifier: <QD vintage>_<COVID method>_<MCMC profile>_<hash>.
#'
#' The hash (first 8 characters of an MD5) covers everything that determines
#' the results apart from the code: the data, sample, model, prior and COVID
#' settings, the series configuration (it sets each variable's prior mean)
#' and the MCMC settings including the seeds. The number of processor cores
#' is not included (it does not change the results). The code version is
#' recorded separately in meta.json (git commit).
#'
#' @param cfg Run configuration.
#' @param method COVID method of this run.
#' @param settings Output of mcmc_settings().
#' @param series Series configuration.
#' @return List: run_id, hash, vintage, hashed (the list that was hashed).
make_run_id <- function(cfg, method, settings, series) {
  covid <- cfg$covid
  covid$methods <- NULL
  covid$method <- method
  hashed <- list(data = cfg$data, sample = cfg$sample, model = cfg$model, prior = cfg$prior,
                 covid = covid, series = as.list(series), mcmc = settings,
                 mode_search = cfg$mcmc$mode_search)
  tmp <- tempfile(fileext = ".yml")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(yaml::as.yaml(hashed), tmp)
  hash <- substr(unname(tools::md5sum(tmp)), 1, 8)
  vintage <- regmatches(cfg$data$qd_vintage_file,
                        regexpr("[0-9]{4}m[0-9]{2}", cfg$data$qd_vintage_file))
  list(run_id = paste(vintage, method, settings$profile, hash, sep = "_"),
       hash = hash, vintage = vintage, hashed = hashed)
}

#' Git commit of the code, and whether tracked files had uncommitted changes.
git_state <- function() {
  commit <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
                     error = function(e) NA_character_, warning = function(w) NA_character_)
  status <- tryCatch(system2("git", c("status", "--porcelain", "--untracked-files=no"),
                             stdout = TRUE, stderr = FALSE),
                     error = function(e) NA_character_, warning = function(w) NA_character_)
  list(commit = if (length(commit)) commit[1] else NA_character_,
       uncommitted_changes = if (length(status) == 1 && is.na(status)) NA else length(status) > 0)
}

#' Write a run's files to <output_root>/<run_id>/.
#'
#' Files:
#'   covid.rds         est_data, origin_data, notes (COVID treatment, R/covid.R)
#'   draws.rds         pool_and_draw() output: the posterior draws the engine uses
#'   estimation.rds    mode search, proposal, chains (all retained draws), timings
#'   mcmc_summary.rds  summarize_mcmc() output
#'   meta.json         run_id, settings, sample, seeds, code version, R and
#'                     package versions, random-number generator
#' A rerun with the same run_id replaces these files.
#'
#' @param id Output of make_run_id().
#' @param output_root Folder that holds the runs (cfg$output_root).
#' @param method COVID method.
#' @param covid Output of covid_treatment().
#' @param mode,proposal,chains,draws,summary Outputs of the estimation steps.
#' @param settings Output of mcmc_settings().
#' @param cfg Run configuration.
#' @return Character vector of the paths written.
write_run_files <- function(id, output_root, method, covid, mode, proposal, chains, draws,
                            summary, settings, cfg) {
  dir <- file.path(output_root, id$run_id)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  paths <- file.path(dir, c("covid.rds", "draws.rds", "estimation.rds", "mcmc_summary.rds",
                            "meta.json"))

  saveRDS(covid, paths[1])
  saveRDS(draws, paths[2])
  saveRDS(list(mode = mode, proposal = proposal, chains = chains), paths[3])
  saveRDS(summary, paths[4])

  fmt <- function(d) format_quarter(d)
  meta <- list(
    run_id = id$run_id, config_hash = id$hash, created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    covid_method = method, mcmc_profile = settings$profile,
    data = list(qd_vintage_file = cfg$data$qd_vintage_file,
                api_realtime_date = cfg$data$api_realtime_date,
                estimation_sample = c(fmt(min(covid$est_data$quarter)), fmt(max(covid$est_data$quarter))),
                estimation_quarters = nrow(covid$est_data),
                origin_last_quarter = fmt(max(covid$origin_data$quarter))),
    model = cfg$model, prior = cfg$prior, covid = id$hashed$covid,
    mcmc = settings, mode_search = cfg$mcmc$mode_search,
    timings_seconds = list(mode = mode$seconds, proposal = proposal$seconds,
                           chains = vapply(chains, function(ch) ch$seconds, numeric(1)),
                           coefficients = draws$seconds),
    code = git_state(),
    r_version = R.version.string, blas = extSoftVersion()[["BLAS"]], rng_kind = RNGkind(),
    packages = as.list(vapply(sort(loadedNamespaces()),
                              function(p) as.character(utils::packageVersion(p)), character(1))),
    files = basename(paths))
  jsonlite::write_json(meta, paths[5], auto_unbox = TRUE, pretty = TRUE, digits = NA)
  paths
}
