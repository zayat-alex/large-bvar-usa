# estimate_checks.R: checks of the prior setup and of the hyperparameter
# posterior computed by BVAR::bvar().
#
# Why it exists: before the full model is estimated, each ingredient of the
# hyperparameter estimation is checked on a small model where every number can
# be inspected: where the prior mean is placed, what the posterior mean is
# under a very tight prior, whether the package's marginal likelihood agrees
# with an independent implementation (R/estimate_ml.R), and how the
# hyperparameter mode depends on the psi hyperprior, the psi starting values
# and the centring of the sum-of-coefficients prior.
#
# Inputs:  a small estimation matrix, the prior configuration, fitted BVAR objects.
# Outputs: data.frames of results (one per check) and base-graphics plots;
#          scripts/check_priors.R runs them and reports/prior_checks.qmd shows them.


# ---- Helpers ----

#' Indices of the psi entries in a named hyperparameter vector.
psi_index <- function(hyper) grep("^psi[0-9]+$", names(hyper))

#' Flag values within 1% (on a log scale) of a lower or upper bound.
#'
#' @return Character vector: "lower", "upper" or "".
at_bound <- function(x, lower, upper) {
  out <- rep("", length(x))
  out[x <= lower * 1.01] <- "lower"
  out[x >= upper / 1.01] <- "upper"
  out
}

#' Fit BVAR only far enough to obtain the hyperparameter mode.
#'
#' BVAR::bvar() always runs the Metropolis step after the mode search; with 20
#' draws it is negligible. The mode (fit$optim) and the stored prior settings
#' are what these checks use.
fit_for_mode <- function(y, lags, priors, seed) {
  set.seed(seed)
  BVAR::bvar(y, lags = lags, n_draw = 20, n_burn = 10, priors = priors, verbose = FALSE)
}


# ---- Check 1: prior mean placement ----

#' Compare the prior mean matrix stored by BVAR with the intended one.
#'
#' @param fit Fitted BVAR object.
#' @param b_own Intended own-first-lag prior means (named, length n).
#' @return data.frame, one row per variable: intended value, value found on the
#'   diagonal of the lag-1 block, and the largest absolute value anywhere else
#'   in that variable's column (should be 0).
check_prior_mean <- function(fit, b_own) {
  b_found <- fit$priors$b                     # K x n
  n_var   <- length(b_own)
  rows <- lapply(seq_len(n_var), function(j) {
    on_diag   <- b_found[1 + j, j]            # row 1 + j = lag 1 of variable j
    elsewhere <- b_found[-(1 + j), j]
    data.frame(variable = names(b_own)[j], intended = b_own[[j]], found = on_diag,
               max_abs_elsewhere = max(abs(elsewhere)),
               ok = on_diag == b_own[[j]] && all(elsewhere == 0))
  })
  do.call(rbind, rows)
}


# ---- Check 2: very tight prior ----

#' Posterior mean of the coefficients when lambda is tiny.
#'
#' As lambda -> 0 the Minnesota prior dominates and the posterior mean of the
#' lag coefficients must approach the prior mean b. Uses BVAR's own posterior
#' mean (beta_hat from BVAR:::bv_ml) at lambda = 1e-4, mu at its prior mode and
#' psi at its starting values.
#'
#' @return List: table (one row per variable: posterior mean of the own first
#'   lag and the largest absolute deviation from b over all its lag
#'   coefficients) and max_abs_dev over everything.
check_tight_prior <- function(fit, b_own, mu, psi, lambda = 1e-4) {
  hyper <- c(lambda = lambda, soc = mu, setNames(psi, paste0("psi", seq_along(psi))))
  res   <- bvar_log_posterior(fit, hyper, opt = FALSE)
  lags  <- fit$meta$lags
  b     <- prior_mean_matrix(b_own, lags)
  dev   <- abs(res$beta_hat[-1, ] - b[-1, ])   # drop the constant (diffuse prior)
  tab <- data.frame(variable = names(b_own),
                    prior_own_lag1 = b_own,
                    posterior_own_lag1 = diag(res$beta_hat[1 + seq_along(b_own), , drop = FALSE]),
                    max_abs_dev_all_lags = apply(dev, 2, max),
                    row.names = NULL)
  list(lambda = lambda, table = tab, max_abs_dev = max(dev))
}


#' Deviation from the prior mean along a decreasing sequence of lambda.
#'
#' The deviation should shrink in proportion to lambda^2 (the prior variance of
#' the coefficients), with any fixed mu.
#'
#' @return data.frame: lambda, max_abs_dev, ratio to the previous row's value.
check_tight_sequence <- function(fit, b_own, mu, psi, lambdas = 10^(-3:-6)) {
  dev <- vapply(lambdas, function(l) check_tight_prior(fit, b_own, mu, psi, l)$max_abs_dev,
                numeric(1))
  data.frame(lambda = lambdas, max_abs_dev = dev,
             ratio_to_previous = c(NA, dev[-1] / dev[-length(dev)]),
             lambda_sq_ratio = c(NA, (lambdas[-1] / lambdas[-length(lambdas)])^2))
}


# ---- Check 3: marginal likelihood cross-check ----

#' Compare BVAR's log posterior with the independent implementation.
#'
#' At each point: BVAR's value; the direct value computed the way BVAR treats
#' the dummy observations (at the prior mean b); and the direct value with the
#' textbook formula (at their posterior mean). Agreement with the first shows
#' the package computes what its source says; the gap to the second is the
#' effect of evaluating the dummy observations at b.
#'
#' @param fit Fitted BVAR object (sum-of-coefficients centred as in BVAR).
#' @param y Estimation matrix.
#' @param b_own Own-first-lag prior means.
#' @param prior_cfg Prior configuration used for the fit.
#' @param points data.frame with columns label, lambda, mu, and a list column
#'   psi (one numeric vector per row).
#' @return data.frame, one row per point.
check_ml_crosscheck <- function(fit, y, b_own, prior_cfg, points) {
  lags <- fit$meta$lags
  ybar <- soc_centring_means(y, lags, prior_cfg$soc_ybar)     # same centring as the fit
  rows <- lapply(seq_len(nrow(points)), function(i) {
    psi   <- points$psi[[i]]
    hyper <- c(lambda = points$lambda[i], soc = points$mu[i],
               setNames(psi, paste0("psi", seq_along(psi))))
    bv    <- bvar_log_posterior(fit, hyper)
    own_b <- log_posterior_direct(y, lags, points$lambda[i], points$mu[i], psi, b_own,
                                  prior_cfg, ybar, dummy_beta_at = "prior_mean")
    own_p <- log_posterior_direct(y, lags, points$lambda[i], points$mu[i], psi, b_own,
                                  prior_cfg, ybar, dummy_beta_at = "posterior")
    data.frame(point = points$label[i], lambda = points$lambda[i], mu = points$mu[i],
               bvar = bv, direct_dummy_at_b = own_b$log_posterior,
               direct_textbook = own_p$log_posterior,
               diff_bvar_vs_dummy_at_b = bv - own_b$log_posterior,
               diff_bvar_vs_textbook = bv - own_p$log_posterior)
  })
  do.call(rbind, rows)
}

#' Log posterior along mu, other hyperparameters fixed: BVAR vs textbook.
#'
#' @param mu_grid Values of mu.
#' @param lambda,psi Fixed values of the other hyperparameters.
#' @return data.frame: mu, bvar, direct_textbook.
profile_mu <- function(fit, y, b_own, prior_cfg, mu_grid, lambda, psi) {
  lags <- fit$meta$lags
  ybar <- soc_centring_means(y, lags, prior_cfg$soc_ybar)
  rows <- lapply(mu_grid, function(mu) {
    hyper <- c(lambda = lambda, soc = mu, setNames(psi, paste0("psi", seq_along(psi))))
    data.frame(mu = mu, bvar = bvar_log_posterior(fit, hyper),
               direct_textbook = log_posterior_direct(y, lags, lambda, mu, psi, b_own,
                                                      prior_cfg, ybar, "posterior")$log_posterior)
  })
  do.call(rbind, rows)
}


# ---- Mode of the textbook log posterior ----

#' Find the hyperparameter mode of the textbook log posterior.
#'
#' Same optimizer and bounds as BVAR::bvar() (L-BFGS-B, but with a higher
#' iteration limit), applied to
#' log_posterior_direct() with the dummy observations treated by the textbook
#' formula.
#'
#' @param start Named starting vector (lambda, soc, psi1..psin).
#' @param lower,upper Bounds, same order.
#' @param ybar Centring of the sum-of-coefficients dummies.
#' @return List: par (named), value, convergence, message.
find_mode_direct <- function(y, lags, b_own, prior_cfg, ybar, start, lower, upper) {
  objective <- function(h) {
    psi <- h[psi_index(start)]
    log_posterior_direct(y, lags, h[["lambda"]], h[["soc"]], psi, b_own, prior_cfg, ybar,
                         "posterior")$log_posterior
  }
  # maxit raised from optim's default of 100: with 2 + n hyperparameters the
  # search can need more iterations; convergence is reported, not assumed.
  opt <- optim(par = start, fn = objective, method = "L-BFGS-B", lower = lower,
               upper = upper, control = list(fnscale = -1, maxit = 2000))
  names(opt$par) <- names(start)
  opt[c("par", "value", "convergence", "message")]
}


# ---- Check 4: psi hyperprior x starting values ----

#' Hyperparameter modes for each psi hyperprior and starting-value option.
#'
#' For each configuration, two modes are found from the same starting values
#' and bounds: BVAR's (fit$optim) and the textbook log posterior's
#' (find_mode_direct()).
#'
#' @param y Estimation matrix; @param series Series configuration (subset).
#' @param prior_cfg Base prior configuration; @param lags Number of lags.
#' @param hyperpriors data.frame: option, shape, scale.
#' @param starts Character vector of psi start options.
#' @param seed Seed for the (negligible) Metropolis step of BVAR::bvar().
#' @return List: summary (one row per configuration and source) and psi (one
#'   row per configuration, source and variable).
check_psi_grid <- function(y, series, prior_cfg, lags, hyperpriors, starts, seed) {
  b_own  <- prior_mean_vector(series)
  ybar   <- soc_centring_means(y, lags, "estimation_sample")   # grid uses the package's centring
  summary_rows <- list()
  psi_rows     <- list()
  for (h in seq_len(nrow(hyperpriors))) {
    for (st in starts) {
      # Step 1: the prior for this configuration.
      pc <- prior_cfg
      pc$psi$shape <- hyperpriors$shape[h]
      pc$psi$scale <- hyperpriors$scale[h]
      pc$psi$start <- st
      pc$soc_ybar  <- "estimation_sample"
      built <- build_priors(pc, series, y, lags)
      tab   <- built$psi_table

      # Step 2: BVAR's mode.
      fit <- fit_for_mode(y, lags, built$priors, seed)
      start <- c(lambda = pc$lambda$mode, soc = pc$soc$mode,
                 setNames(tab$start, paste0("psi", seq_len(nrow(tab)))))
      lower <- c(lambda = pc$lambda$min, soc = pc$soc$min, setNames(tab$min, names(start)[-(1:2)]))
      upper <- c(lambda = pc$lambda$max, soc = pc$soc$max, setNames(tab$max, names(start)[-(1:2)]))

      # Step 3: the textbook mode from the same start and bounds.
      direct <- find_mode_direct(y, lags, b_own, pc, ybar, start, lower, upper)

      # Step 4: record both.
      for (src in c("bvar", "textbook")) {
        if (src == "bvar") {
          par <- fit$optim$par; value <- fit$optim$value
          conv <- fit$optim$convergence; msg <- fit$optim$message
          textbook_value <- log_posterior_direct(y, lags, par[["lambda"]], par[["soc"]],
                                                 par[psi_index(par)], b_own, pc, ybar,
                                                 "posterior")$log_posterior
        } else {
          par <- direct$par; value <- direct$value
          conv <- direct$convergence; msg <- direct$message
          textbook_value <- value
        }
        bound <- at_bound(par, lower, upper)
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          hyperprior = hyperpriors$option[h], psi_start = st, source = src,
          lambda = par[["lambda"]], mu = par[["soc"]],
          log_post_own_objective = value, log_post_textbook = textbook_value,
          convergence = conv, message = msg,
          at_bound = paste(names(par)[bound != ""], bound[bound != ""], collapse = ", "),
          stringsAsFactors = FALSE)
        psi <- par[psi_index(par)]
        psi_rows[[length(psi_rows) + 1]] <- data.frame(
          hyperprior = hyperpriors$option[h], psi_start = st, source = src,
          variable = colnames(y), psi = unname(psi), ar_variance = tab$ar_variance,
          ratio_to_ar_variance = unname(psi) / tab$ar_variance,
          start = tab$start, min = tab$min, max = tab$max,
          at_bound = bound[psi_index(par)], stringsAsFactors = FALSE)
      }
    }
  }
  list(summary = do.call(rbind, summary_rows), psi = do.call(rbind, psi_rows))
}


# ---- Check 5: centring of the sum-of-coefficients dummies ----

#' Textbook mode under the two centring conventions.
#'
#' "estimation_sample": ybar = mean of rows p+1..2p (BVAR's convention).
#' "presample": ybar = mean of rows 1..p (Giannone, Lenza and Primiceri 2015).
#'
#' @return data.frame, one row per convention: ybar summary, mode and value.
check_soc_centring <- function(y, series, prior_cfg, lags) {
  b_own <- prior_mean_vector(series)
  built <- build_priors(prior_cfg, series, y, lags)
  tab   <- built$psi_table
  start <- c(lambda = prior_cfg$lambda$mode, soc = prior_cfg$soc$mode,
             setNames(tab$start, paste0("psi", seq_len(nrow(tab)))))
  lower <- c(lambda = prior_cfg$lambda$min, soc = prior_cfg$soc$min, setNames(tab$min, names(start)[-(1:2)]))
  upper <- c(lambda = prior_cfg$lambda$max, soc = prior_cfg$soc$max, setNames(tab$max, names(start)[-(1:2)]))
  conventions <- list(estimation_sample = lags + seq_len(lags), presample = seq_len(lags))
  rows <- lapply(names(conventions), function(cv) {
    ybar <- soc_centring_means(y, lags, cv)
    m <- find_mode_direct(y, lags, b_own, prior_cfg, ybar, start, lower, upper)
    data.frame(convention = cv, rows_used = paste(rownames(y)[range(conventions[[cv]])], collapse = " to "),
               lambda = m$par[["lambda"]], mu = m$par[["soc"]], log_post = m$value,
               convergence = m$convergence,
               t(setNames(m$par[psi_index(m$par)], paste0("psi_", colnames(y)))),
               check.names = FALSE)
  })
  do.call(rbind, rows)
}


# ---- Plots (base graphics) ----

# Colours: categorical slots 1-3 of the reference palette (validated as a set),
# ink colours for text, a light grey for the grid.
check_colours <- c(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a")
ink <- c(primary = "#0b0b0b", secondary = "#52514e", grid = "#e4e3df")

#' Log density of the three psi hyperpriors over a range of psi.
#'
#' @param hyperpriors data.frame: option, shape, scale.
#' @param ar_var AR(p) residual variances of the variables (marked as ticks).
plot_psi_hyperpriors <- function(hyperpriors, ar_var) {
  psi <- 10^seq(-5, 4, length.out = 400)
  dens <- sapply(seq_len(nrow(hyperpriors)), function(h)
    log_dinvgamma(psi, hyperpriors$shape[h], hyperpriors$scale[h]))
  # The y range is cut at 30 log points below the highest value: further down,
  # one curve falls by thousands of log points and would flatten the others.
  ylim <- c(max(dens) - 30, max(dens) + 1)
  par(mar = c(5.5, 4.5, 2.5, 1), col.axis = ink[["secondary"]], col.lab = ink[["primary"]])
  plot(NA, xlim = range(psi), ylim = ylim, log = "x",
       xlab = "psi (log scale); ticks: AR(p) residual variances of the variables",
       ylab = "log density", main = "Inverse-gamma hyperpriors on psi", font.main = 1,
       col.main = ink[["primary"]], bty = "n")
  abline(v = 10^(-5:4), col = ink[["grid"]], lwd = 1)
  for (h in seq_len(nrow(hyperpriors))) {
    lines(psi, dens[, h], col = check_colours[h], lwd = 2)
  }
  rug(ar_var, col = ink[["secondary"]], lwd = 2)
  legend("topright", legend = paste0(hyperpriors$option, " (shape ", hyperpriors$shape,
                                       ", scale ", hyperpriors$scale, ")"),
         col = check_colours[seq_len(nrow(hyperpriors))], lwd = 2, bty = "n",
         text.col = ink[["primary"]])
}

#' Log posterior along mu: BVAR vs textbook formula.
plot_mu_profile <- function(prof, mu_mode_bvar) {
  par(mar = c(4.5, 4.5, 2.5, 1), col.axis = ink[["secondary"]], col.lab = ink[["primary"]])
  rng <- range(c(prof$bvar, prof$direct_textbook))
  plot(NA, xlim = range(prof$mu), ylim = rng, log = "x", xlab = "mu (log scale)",
       ylab = "log posterior (up to a constant)",
       main = "Log posterior along mu, other hyperparameters at BVAR's mode",
       font.main = 1, col.main = ink[["primary"]], bty = "n")
  abline(v = 10^(-4:2), col = ink[["grid"]])
  lines(prof$mu, prof$bvar, col = check_colours[["blue"]], lwd = 2)
  lines(prof$mu, prof$direct_textbook, col = check_colours[["orange"]], lwd = 2)
  legend("bottomleft", legend = c("BVAR (dummies evaluated at b)", "textbook formula"),
         col = check_colours[c("blue", "orange")], lwd = 2, bty = "n", text.col = ink[["primary"]])
}

#' psi at the mode relative to the AR(p) residual variance, per configuration.
#'
#' One panel per source (BVAR, textbook); x = variable, y = psi / AR variance
#' (log scale); colour = hyperprior option; point shape = starting-value option.
plot_psi_modes <- function(psi_tab, hyperpriors) {
  vars <- unique(psi_tab$variable)
  starts <- unique(psi_tab$psi_start)
  op <- par(mfrow = c(1, 2), mar = c(6, 4.5, 2.5, 1), oma = c(0, 0, 0, 11),
            col.axis = ink[["secondary"]], col.lab = ink[["primary"]])
  on.exit(par(op))
  for (src in c("bvar", "textbook")) {
    d <- psi_tab[psi_tab$source == src, ]
    plot(NA, xlim = c(0.5, length(vars) + 0.5), ylim = range(psi_tab$ratio_to_ar_variance),
         log = "y", xaxt = "n", xlab = "", ylab = "psi / AR(p) residual variance",
         main = paste0("Mode: ", src), font.main = 1, col.main = ink[["primary"]], bty = "n")
    abline(h = 10^(-3:3), col = ink[["grid"]])
    abline(h = 1, col = ink[["secondary"]], lty = 2)
    axis(1, at = seq_along(vars), labels = vars, las = 2)
    for (h in seq_len(nrow(hyperpriors))) {
      for (s in seq_along(starts)) {
        dd <- d[d$hyperprior == hyperpriors$option[h] & d$psi_start == starts[s], ]
        x  <- match(dd$variable, vars) + (h - 2) * 0.18 + (s - 1.5) * 0.06
        points(x, dd$ratio_to_ar_variance, col = check_colours[h], pch = c(16, 1)[s], cex = 1.4)
      }
    }
  }
  par(xpd = NA)    # allow the legend in the outer margin, outside both panels
  legend(x = grconvertX(1, from = "ndc"), y = grconvertY(0.8, from = "ndc"), xjust = 1,
         legend = c(hyperpriors$option, paste("start:", starts)),
         col = c(check_colours[seq_len(nrow(hyperpriors))], ink[["secondary"]], ink[["secondary"]]),
         pch = c(rep(16, nrow(hyperpriors)), 16, 1), bty = "n", text.col = ink[["primary"]], cex = 0.85)
}


# ---- The package's value, for comparison ----

#' BVAR's log posterior kernel of the hyperparameters at given values.
#'
#' Calls the package's internal function BVAR:::bv_ml() with the prior settings
#' stored in a fitted object, exactly as BVAR::bvar() does during the
#' Metropolis step.
#'
#' @param fit A fitted BVAR object (its prior settings and data are used).
#' @param hyper Named numeric vector: lambda, soc, psi1, ..., psin.
#' @param opt TRUE to return only the value; FALSE for the full list
#'   (includes beta_hat, the posterior mean of the coefficients).
bvar_log_posterior <- function(fit, hyper, opt = TRUE) {
  priors <- fit$priors
  # The fixed parameters (all hyperparameters at their modes), named as BVAR
  # names them; the values in `hyper` overwrite the estimated ones.
  pars_names <- names(priors)[!grepl("^hyper$|^var$|^b$|^psi[0-9]+$|^dummy$", names(priors))]
  pars_full  <- do.call(c, lapply(pars_names, function(x) priors[[x]][["mode"]]))
  names(pars_full) <- BVAR:::name_pars(pars_names, fit$meta$M)
  X <- fit$meta$X
  BVAR:::bv_ml(hyper = hyper, hyper_min = -Inf, hyper_max = Inf, pars = pars_full,
               priors = priors, Y = fit$meta$Y, X = X, XX = crossprod(X),
               K = fit$meta$K, M = fit$meta$M, N = fit$meta$N, lags = fit$meta$lags,
               opt = opt)
}
