# prior_sensitivity_checks.R: how much do forecasts and impulse responses
# move when the hyperparameters are moved away from their mode?
#
# Why it exists: the hyperparameters (lambda, overall tightness; mu, the
# sum-of-coefficients tightness; psi, the scale of each variable's residual
# variance in the prior) are estimated (Crump et al. 2021, Section 2.2), so
# the results carry their uncertainty. This check measures how sensitive the
# results are to them: one at a time, a hyperparameter is set to a stressed
# value (half or twice its mode) or to the edge of its posterior (mode plus
# or minus two posterior standard deviations), the coefficients are drawn
# given that value, and forecasts and impulse responses are compared with
# those at the mode. Read together with the stability check across origins
# (R/hyper_stability_checks.R): if results respond little to hyperparameter
# shifts of the size seen between origins, the hyperparameters need not be
# re-estimated at every origin.
#
# The comparison machinery (draws given fixed hyperparameters, band
# comparison, figures) is shared with R/hyper_stability_checks.R.
#
# Inputs:  the estimation setup, the current mode and posterior draws.
# Outputs: the list of variants; figures used by reports/prior_sensitivity.qmd.


#' Hyperparameter vectors to compare with the mode.
#'
#' All on theta = log(h). Stress variants multiply h by 1/2 or 2 (theta
#' -/+ log 2); for psi all 31 values are scaled together. Posterior-range
#' variants move one of lambda, mu by -/+ 2 posterior s.d. of theta.
#'
#' @param theta_mode Named mode on the log scale (lambda, soc, psi1..psin).
#' @param sd_theta Named posterior s.d. of theta (from the MCMC draws).
#' @return Named list of theta vectors, with attribute "info" (data.frame:
#'   variant, hyperparameter, kind, change in theta).
sensitivity_variants <- function(theta_mode, sd_theta) {
  psi <- grepl("^psi[0-9]+$", names(theta_mode))
  shift <- function(which, by) { th <- theta_mode; th[which] <- th[which] + by; th }
  v <- list(
    `lambda x 0.5` = shift("lambda", -log(2)), `lambda x 2` = shift("lambda", log(2)),
    `mu x 0.5`     = shift("soc", -log(2)),    `mu x 2`     = shift("soc", log(2)),
    `psi x 0.5`    = shift(psi, -log(2)),      `psi x 2`    = shift(psi, log(2)),
    `lambda -2 sd` = shift("lambda", -2 * sd_theta[["lambda"]]),
    `lambda +2 sd` = shift("lambda",  2 * sd_theta[["lambda"]]),
    `mu -2 sd`     = shift("soc", -2 * sd_theta[["soc"]]),
    `mu +2 sd`     = shift("soc",  2 * sd_theta[["soc"]]))
  attr(v, "info") <- data.frame(
    variant = names(v),
    hyperparameter = rep(c("lambda", "mu", "psi (all)", "lambda", "mu"), each = 2),
    kind = rep(c("stress", "posterior range"), c(6, 4)),
    theta_change = c(-log(2), log(2), -log(2), log(2), -log(2), log(2),
                     -2 * sd_theta[["lambda"]], 2 * sd_theta[["lambda"]],
                     -2 * sd_theta[["soc"]], 2 * sd_theta[["soc"]]))
  v
}

#' Difference of medians by horizon for every variant, one variable per panel.
#'
#' Colour = hyperparameter (lambda blue, mu orange, psi aqua); solid = moved
#' up, dashed = moved down; thick = stress (x 2, x 0.5), thin = posterior
#' range (+/- 2 s.d.). The grey band is the largest seed-noise difference.
#'
#' @param cmp Named list of compare_bands() results, one per variant.
#' @param variable Variable name.
#' @param title Panel title.
plot_sensitivity_panel <- function(cmp, variable, title) {
  d <- lapply(cmp, function(x) x[x$variable == variable & !is.na(x$median_diff_share), ])
  noise <- max(sapply(d, function(x) max(abs(x$median_diff_share_seed))))
  ylim <- range(c(-noise, noise, unlist(lapply(d, `[[`, "median_diff_share"))))
  h <- d[[1]]$horizon
  plot(NA, xlim = range(h), ylim = ylim, bty = "n", xlab = "quarters", ylab = "",
       main = title, font.main = 1, cex.main = 0.9, col.main = ink[["primary"]],
       col.axis = ink[["secondary"]])
  polygon(c(h, rev(h)), c(rep(-noise, length(h)), rep(noise, length(h))),
          col = "#e4e3df99", border = NA)
  abline(h = 0, col = ink[["secondary"]])
  col_of <- c(lambda = check_colours[["blue"]], mu = check_colours[["orange"]],
              psi = check_colours[["aqua"]])
  for (nm in names(d)) {
    par_name <- sub(" .*", "", nm)
    up <- grepl("x 2|\\+2", nm)
    stress <- grepl(" x ", nm)
    lines(d[[nm]]$horizon, d[[nm]]$median_diff_share, col = col_of[[par_name]],
          lty = if (up) 1 else 2, lwd = if (stress) 2 else 1)
  }
}
