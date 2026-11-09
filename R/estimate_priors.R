# estimate_priors.R: build the prior of the Bayesian VAR from the run
# configuration (estimation stage).
#
# Why it exists: the prior is where the model's assumptions live (Crump et al.
# 2021, Section 2.2; Giannone, Lenza and Primiceri 2015, Section 3). Building it
# in one place from config/run.yml and config/series.csv keeps every setting
# visible and makes the prior object that goes into BVAR::bvar() reproducible.
#
# The prior (normal-inverse-Wishart):
#   Sigma            ~ IW(Psi, d),  d = n + 2,  Psi = diag(psi_1, ..., psi_n)
#   vec(B) | Sigma   ~ N(vec(b), Sigma (x) Omega)
#   Minnesota prior: prior mean b = 1 on each variable's own first lag for
#     random-walk (RW) variables and 0 for white-noise (WN) variables (only
#     realized volatility); 0 on everything else. Prior variance factor of the
#     coefficient on variable j at lag s: lambda^2 / (s^alpha * psi_j).
#   Sum-of-coefficients prior: n dummy observations, tightness mu.
#   Hyperpriors: lambda ~ Gamma, mu ~ Gamma (given by mode and sd),
#     psi_j ~ inverse-gamma(shape, scale).
#
# The formulas and settings are set out in docs/methodology.md (section 4);
# code comments cite its equation labels (M1), (M2), ...
#
# Inputs:  the `prior` and `model` sections of the run configuration, the
#          series configuration, and the estimation data.
# Outputs: a list with the BVAR prior object and the intermediate values used
#          to build it (for inspection and reporting).


# ---- Data as a matrix ----

#' Turn the estimation data.frame into the numeric matrix BVAR expects.
#'
#' @param est_data data.frame: `quarter` plus the model variables.
#' @param series Series configuration (gives the expected column order).
#' @return Numeric matrix, one row per quarter, one column per variable, in the
#'   order of config/series.csv; row names are the quarter labels.
est_data_to_matrix <- function(est_data, series) {
  # Step 1: the columns must be exactly the configured variables, in order,
  # because the column position defines each variable's place in the VAR.
  vars <- setdiff(names(est_data), "quarter")
  if (!identical(vars, series$short_name)) {
    stop("est_data columns do not match config/series.csv (names or order)")
  }

  # Step 2: drop the quarter column and convert to a numeric matrix.
  y <- as.matrix(est_data[, vars, drop = FALSE])
  storage.mode(y) <- "double"
  rownames(y) <- format_quarter(est_data$quarter)
  y
}


# ---- Minnesota prior mean ----

#' Prior mean of each variable's own first-lag coefficient.
#'
#' Implements docs/methodology.md (M6).
#'
#' @param series Series configuration (column `prior`: "RW" or "WN").
#' @return Named numeric vector, 1 for RW and 0 for WN, in series order.
prior_mean_vector <- function(series) {
  b <- ifelse(series$prior == "RW", 1, 0)
  if (any(!series$prior %in% c("RW", "WN"))) stop("prior must be RW or WN")
  names(b) <- series$short_name
  b
}


# ---- Starting values and bounds of psi ----

#' Residual variance of an AR(p) regression with a constant, by OLS.
#'
#' Used to set the scale of each variable: psi_j is the prior mean of the
#' innovation variance Sigma_jj (with d = n + 2), so a natural starting value
#' is the residual variance of a univariate autoregression on the same data.
#'
#' @param x Numeric vector (one variable, in model units).
#' @param lags Number of lags p.
#' @param keep NULL, or a logical vector of length T: the observations t whose
#'   equation is used (as in var_data_matrices(); excluded quarters still
#'   serve as lags).
#' @return List: variance (sum of squared residuals divided by the degrees of
#'   freedom, n_obs - (p + 1), where n_obs = T - p observations are used and
#'   p + 1 coefficients are estimated; fewer when quarters are excluded), sd
#'   (its square root), n_obs.
ar_residual_variance <- function(x, lags, keep = NULL) {
  # Step 1: build the regression y_t = c + a_1 y_{t-1} + ... + a_p y_{t-p} + e_t
  # on the observations t = p+1, ..., T.
  n_total <- length(x)
  y_dep   <- x[(lags + 1):n_total]
  x_reg   <- matrix(1, nrow = n_total - lags, ncol = lags + 1)   # column 1: constant
  for (s in seq_len(lags)) {
    x_reg[, s + 1] <- x[(lags + 1 - s):(n_total - s)]
  }
  if (!is.null(keep)) {
    use   <- keep[(lags + 1):n_total]
    y_dep <- y_dep[use]
    x_reg <- x_reg[use, , drop = FALSE]
  }

  # Step 2: OLS coefficients and residuals.
  coef  <- solve(crossprod(x_reg), crossprod(x_reg, y_dep))
  resid <- y_dep - x_reg %*% coef

  # Step 3: residual variance with a degrees-of-freedom correction.
  n_obs <- length(y_dep)
  dof   <- n_obs - ncol(x_reg)
  variance <- sum(resid^2) / dof
  list(variance = variance, sd = sqrt(variance), n_obs = n_obs)
}

#' Table of psi starting values and search bounds, one row per variable.
#'
#' Two options:
#'   "ar_variance": start at the OLS AR(p) residual variance of each variable
#'     (psi is a variance scale, so this matches its units); bounds are
#'     start / bound_factor and start * bound_factor.
#'   "ar_sd_auto": the package default (BVAR's auto_psi): start at the residual
#'     standard deviation of an AR(p) fitted by maximum likelihood (arima) on
#'     the estimation rows after the initial lags; bounds start / 100 and
#'     start * 100. Recomputed here with the package's own function so the
#'     values can be shown; the prior then uses mode = "auto".
#'
#' Implements docs/methodology.md section 4.4.
#'
#' @param y Numeric matrix of estimation data (est_data_to_matrix()).
#' @param lags Number of lags p.
#' @param start "ar_variance" or "ar_sd_auto".
#' @param bound_factor Ratio between the start and each bound.
#' @param keep NULL, or the observations whose equation is used
#'   (ar_residual_variance()); only with start = "ar_variance".
#' @return data.frame: variable, ar_variance, ar_sd (OLS), start, min, max.
psi_start_table <- function(y, lags, start, bound_factor, keep = NULL) {
  if (!is.null(keep) && start != "ar_variance") stop("excluded quarters need prior$psi$start = ar_variance")
  # Step 1: OLS AR(p) residual variance and sd of every variable (shown under
  # both options, for comparison).
  ar <- lapply(seq_len(ncol(y)), function(j) ar_residual_variance(y[, j], lags, keep))
  tab <- data.frame(
    variable    = colnames(y),
    ar_variance = vapply(ar, function(a) a$variance, numeric(1)),
    ar_sd       = vapply(ar, function(a) a$sd, numeric(1)),
    stringsAsFactors = FALSE
  )

  # Step 2: starting values and bounds under the chosen option.
  if (start == "ar_variance") {
    tab$start <- tab$ar_variance
    tab$min   <- tab$start / bound_factor
    tab$max   <- tab$start * bound_factor
  } else if (start == "ar_sd_auto") {
    # BVAR's auto_psi() works on the rows after the first `lags` (the rows used
    # as left-hand side in the VAR), as inside BVAR::bvar().
    y_lhs <- y[(lags + 1):nrow(y), , drop = FALSE]
    auto  <- suppressMessages(BVAR:::auto_psi(y_lhs, lags))
    tab$start <- auto$mode
    tab$min   <- auto$min
    tab$max   <- auto$max
  } else {
    stop("prior$psi$start must be ar_variance or ar_sd_auto")
  }
  tab
}


# ---- Sum-of-coefficients dummies ----

#' Sum-of-coefficients dummy observations centred on the first p observations.
#'
#' Giannone, Lenza and Primiceri (2015) centre the dummies on ybar0, the mean
#' of the first p observations of the sample. BVAR's built-in version uses the
#' first p rows after the initial lags instead. This function returns a dummy
#' generator in the format BVAR::bv_dummy() expects, with ybar0 fixed from the
#' sample start:
#'   Y+ = diag(ybar0) / mu                 (n x n)
#'   X+ = [0, Y+, Y+, ..., Y+]             (n x (1 + n p), Y+ repeated p times)
#'
#' @param ybar0 Numeric vector: mean of the first p rows of the estimation data.
#' @return function(Y, lags, par) returning list(Y = Y+, X = X+); `par` is mu.
make_soc_presample_dummy <- function(ybar0) {
  force(ybar0)
  function(Y, lags, par) {
    y_plus <- diag(ybar0, nrow = length(ybar0)) / par
    x_plus <- cbind(0, do.call(cbind, rep(list(y_plus), lags)))
    list(Y = y_plus, X = x_plus)
  }
}


#' Centring means of the sum-of-coefficients dummies.
#'
#' "presample": mean of the first p rows of the estimation data (Giannone,
#'   Lenza and Primiceri 2015).
#' "estimation_sample": mean of rows p+1 to 2p, the first p rows used as
#'   left-hand side (BVAR's built-in convention).
#'
#' Implements docs/methodology.md (M9).
#'
#' @param y Numeric estimation matrix T x n.
#' @param lags Number of lags p.
#' @param option "presample" or "estimation_sample".
#' @return Numeric vector of length n.
soc_centring_means <- function(y, lags, option) {
  if (option == "presample") {
    rows <- seq_len(lags)
  } else if (option == "estimation_sample") {
    rows <- lags + seq_len(lags)
  } else {
    stop("prior$soc_ybar must be estimation_sample or presample")
  }
  colMeans(y[rows, , drop = FALSE])
}


# ---- The prior object ----

#' Build the BVAR prior from the configuration.
#'
#' @param prior_cfg The `prior` section of the run configuration.
#' @param series Series configuration.
#' @param y Numeric matrix of estimation data (est_data_to_matrix()).
#' @param lags Number of lags p.
#' @return List:
#'   priors    : the BVAR::bv_priors() object passed to BVAR::bvar();
#'   b         : prior mean of each variable's own first lag;
#'   psi_table : psi starting values and bounds (psi_start_table());
#'   soc_ybar  : the soc_ybar option used, and ybar0 if "presample".
build_priors <- function(prior_cfg, series, y, lags) {
  # Step 1: prior mean of the own first lag (1 = random walk, 0 = white noise).
  b <- prior_mean_vector(series)

  # Step 2: psi starting values and bounds.
  psi_table <- psi_start_table(y, lags, prior_cfg$psi$start, prior_cfg$psi$bound_factor)
  if (prior_cfg$psi$start == "ar_variance") {
    psi <- BVAR::bv_psi(scale = prior_cfg$psi$scale, shape = prior_cfg$psi$shape,
                        mode = psi_table$start, min = psi_table$min, max = psi_table$max)
  } else {
    psi <- BVAR::bv_psi(scale = prior_cfg$psi$scale, shape = prior_cfg$psi$shape,
                        mode = "auto")
  }

  # Step 3: lambda (estimated) and alpha (fixed: listed with sd 0.25 and
  # bounds [1, 3], the package defaults, but not estimated, so only its mode
  # is used).
  lam <- prior_cfg$lambda
  lambda <- BVAR::bv_lambda(mode = lam$mode, sd = lam$sd, min = lam$min, max = lam$max)
  alpha  <- BVAR::bv_alpha(mode = prior_cfg$alpha)

  # Step 4: the Minnesota prior.
  mn <- BVAR::bv_minnesota(lambda = lambda, alpha = alpha, psi = psi,
                           var = prior_cfg$constant_var, b = unname(b))

  # Step 5: the sum-of-coefficients prior (mu), centred on either convention.
  s <- prior_cfg$soc
  if (prior_cfg$soc_ybar == "estimation_sample") {
    soc <- BVAR::bv_soc(mode = s$mode, sd = s$sd, min = s$min, max = s$max)
    ybar0 <- NULL
  } else if (prior_cfg$soc_ybar == "presample") {
    ybar0 <- soc_centring_means(y, lags, "presample")
    soc <- BVAR::bv_dummy(mode = s$mode, sd = s$sd, min = s$min, max = s$max,
                          fun = make_soc_presample_dummy(ybar0))
  } else {
    stop("prior$soc_ybar must be estimation_sample or presample")
  }

  # Step 6: put them together; lambda, soc (mu) and psi are estimated.
  priors <- BVAR::bv_priors(hyper = c("lambda", "soc", "psi"), mn = mn, soc = soc)

  list(priors = priors, b = b, psi_table = psi_table,
       soc_ybar = list(option = prior_cfg$soc_ybar, ybar0 = ybar0))
}
