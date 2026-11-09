# mcmc_diagnostics.R: figures and report of an estimation run's MCMC
# diagnostics (estimation stage, after the chains).
#
# Why it exists: the hyperparameter chain is the only Markov-chain part of
# the estimator (Giannone, Lenza and Primiceri 2015, Appendix B; the
# coefficients are drawn exactly given each hyperparameter draw), so whether
# the run's posterior draws can be used depends on that chain having
# converged and mixed. summarize_mcmc() (R/estimate_run.R) computes the
# numbers; this file draws the pictures a reader needs to judge them:
#   - trace plots of every hyperparameter and of the log posterior, chains
#     overlaid, burn-in shaded (chains that started apart and now overlap
#     have forgotten their start);
#   - the tuning path of the proposal scale c during the burn-in;
#   - the posterior of each hyperparameter within its search bounds, and of
#     lambda and mu against their hyperpriors;
#   - the largest companion-matrix eigenvalue modulus across kept draws;
#   - a derived quantity (rgdp growth over 4 quarters, no-shock path) by chain.
# The report reports/mcmc_diagnostics.qmd shows these with the tables of
# summarize_mcmc(); render_mcmc_report() renders it into the run's folder.
#
# Notation as in R/estimate_bvar.R: h = (lambda, soc, psi1..psin), theta = log(h).
#
# The formulas and settings are set out in docs/methodology.md (section 8.2);
# code comments cite its equation labels (M1), (M2), ...
#
# Inputs:  chains (run_chains()), the summary (summarize_mcmc()), the setup
#          (estimation_setup()).
# Outputs: PNG files in <output_root>/<run_id>/plots/, and
#          <output_root>/<run_id>/mcmc_diagnostics.html.


# Colours of the chains: categorical slots 1-8 of the reference palette, in
# fixed order (slots 1-4 validated as a set for the four-chain runs). The
# palette extends check_colours (R/estimate_checks.R).
chain_colours <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
                   "#6250d6", "#e34948")


#' Colours for k chains, in fixed order.
#'
#' @param k Number of chains (at most 8; more chains would need small
#'   multiples rather than more hues).
#' @return Character vector of k colours.
colours_for_chains <- function(k) {
  if (k > length(chain_colours)) stop("at most ", length(chain_colours), " chains can be coloured")
  chain_colours[seq_len(k)]
}

#' The whole recorded path of one chain: burn-in, then retained draws.
#'
#' @param ch One element of run_chains() output.
#' @param column Column of theta, or "log_post" for the log posterior.
#' @return data.frame: iteration (1..n_draw, every n_thin-th), value, burn
#'   (TRUE during the burn-in).
chain_path <- function(ch, column) {
  thin <- ch$settings$n_thin
  if (column == "log_post") {
    burn <- ch$burn_log_target; kept <- ch$log_target
  } else {
    burn <- ch$burn_theta[, column]; kept <- ch$theta[, column]
  }
  n_b <- length(burn)
  data.frame(iteration = c(seq_len(n_b), n_b + seq_along(kept)) * thin,
             value = c(burn, kept), burn = c(rep(TRUE, n_b), rep(FALSE, length(kept))))
}

#' Trace plots: one panel per quantity, chains overlaid, burn-in shaded.
#'
#' Values are on the sampler's scale: log(h) for the hyperparameters and log
#' p(theta | Y) for the log posterior.
#'
#' @param chains Output of run_chains().
#' @param columns Columns to plot (names of theta, or "log_post").
#' @param labels Display names, one per column.
#' @param n_col Panels per row.
#' @param max_points Upper limit on points drawn per chain and panel
#'   (evenly spaced), to keep the files small; the statistics use every draw.
plot_traces <- function(chains, columns, labels, n_col = 4, max_points = 3000) {
  cols <- colours_for_chains(length(chains))
  n_row <- ceiling(length(columns) / n_col)
  op <- par(mfrow = c(n_row, n_col), mar = c(2, 3, 1.8, 0.5), oma = c(2, 0, 0, 0),
            col.axis = ink[["secondary"]], cex.axis = 0.7)
  on.exit(par(op))
  n_burn <- chains[[1]]$settings$n_burn
  for (j in seq_along(columns)) {
    paths <- lapply(chains, chain_path, column = columns[j])
    yr <- range(unlist(lapply(paths, function(d) d$value)), finite = TRUE)
    xr <- range(unlist(lapply(paths, function(d) d$iteration)))
    plot(NA, xlim = xr, ylim = yr, bty = "n", xlab = "", ylab = "", main = labels[j],
         font.main = 1, cex.main = 0.8, col.main = ink[["primary"]])
    rect(0, yr[1], n_burn, yr[2], col = ink[["grid"]], border = NA)
    for (k in seq_along(paths)) {
      d <- paths[[k]]
      keep <- unique(round(seq(1, nrow(d), length.out = min(nrow(d), max_points))))
      lines(d$iteration[keep], d$value[keep], col = cols[k], lwd = 0.6)
    }
  }
  mtext(sprintf("iteration (grey: burn-in, %d iterations); chains: %s", n_burn,
                paste(sprintf("%d", seq_along(chains)), collapse = ", ")),
        side = 1, outer = TRUE, line = 0.5, cex = 0.7, col = ink[["secondary"]])
}

#' Legend strip naming the chain colours (for the report, above trace plots).
#'
#' @param k Number of chains.
plot_chain_legend <- function(k) {
  cols <- colours_for_chains(k)
  op <- par(mar = c(0, 0, 0, 0))
  on.exit(par(op))
  plot.new()
  legend("center", legend = paste("chain", seq_len(k)), col = cols, lwd = 2, horiz = TRUE,
         bty = "n", text.col = ink[["primary"]], cex = 0.9)
}

#' Tuning path of the proposal scale c during the burn-in, by chain.
#'
#' @param chains Output of run_chains().
plot_tuning <- function(chains) {
  cols <- colours_for_chains(length(chains))
  paths <- lapply(chains, function(ch) ch$scale_path)
  if (all(vapply(paths, is.null, logical(1)))) return(invisible(NULL))
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 1, 1), col.axis = ink[["secondary"]],
            col.lab = ink[["primary"]])
  on.exit(par(op))
  xr <- range(unlist(lapply(paths, function(p) p[, "iteration"])))
  plot(NA, xlim = xr, ylim = range(unlist(lapply(paths, function(p) p[, "scale"]))), log = "y",
       bty = "n", xlab = "iteration (end of tuning window)", ylab = "proposal scale c (log axis)")
  for (k in seq_along(paths)) lines(paths[[k]][, "iteration"], paths[[k]][, "scale"], col = cols[k], lwd = 2)
  plot(NA, xlim = xr, ylim = c(0, 1), bty = "n", xlab = "iteration (end of tuning window)",
       ylab = "acceptance in the window")
  target <- chains[[1]]$settings$target
  rect(xr[1], target[1], xr[2], target[2], col = ink[["grid"]], border = NA)
  for (k in seq_along(paths)) lines(paths[[k]][, "iteration"], paths[[k]][, "acceptance"], col = cols[k], lwd = 2)
  legend("topright", legend = paste("chain", seq_along(chains)), col = cols, lwd = 2, bty = "n",
         cex = 0.8, text.col = ink[["primary"]])
}

#' Posterior of every hyperparameter within its search bounds.
#'
#' One row per hyperparameter, on the scale log10(h / start), where start is
#' the configured starting value (lambda and mu: hyperprior modes; psi: AR(p)
#' residual variances). The grey bar spans the search bounds, the blue bar the
#' 5-95 percent posterior interval, the dot the median. A posterior pressed
#' against the end of its grey bar would be shaped by the bound, not the data.
#'
#' @param summary Output of summarize_mcmc().
#' @param setup Output of estimation_setup().
#' @param labels Display names, one per hyperparameter (in setup order).
plot_hyper_position <- function(summary, setup, labels) {
  nm <- names(setup$start)
  h <- summary$hyper[match(nm, summary$hyper$hyperparameter), ]
  s <- setup$start[nm]
  pos <- function(x) log10(x / s)
  lo <- pos(setup$lower[nm]); hi <- pos(setup$upper[nm])
  k <- length(nm); yy <- rev(seq_len(k))
  op <- par(mar = c(4, 8, 1, 1), col.axis = ink[["secondary"]], col.lab = ink[["primary"]])
  on.exit(par(op))
  plot(NA, xlim = range(c(lo, hi)), ylim = c(0.5, k + 0.5), yaxt = "n", bty = "n",
       xlab = "log10(hyperparameter / starting value)", ylab = "")
  segments(lo, yy, hi, yy, col = ink[["grid"]], lwd = 6)
  abline(v = 0, col = ink[["secondary"]], lty = 3)
  segments(pos(h$q05), yy, pos(h$q95), yy, col = check_colours[["blue"]], lwd = 3)
  points(pos(h$median), yy, pch = 16, cex = 0.8, col = check_colours[["blue"]])
  axis(2, at = yy, labels = labels, las = 1, cex.axis = 0.7, tick = FALSE)
}

#' Posterior of lambda and mu against their Gamma hyperpriors.
#'
#' Histograms of the retained draws (all chains) on the original scale, with
#' the hyperprior density (Giannone, Lenza and Primiceri 2015, Section 3;
#' config prior$lambda, prior$soc) drawn over the same range.
#'
#' @param chains Output of run_chains().
#' @param prior_cfg The `prior` section of the run configuration.
plot_lambda_mu_prior <- function(chains, prior_cfg) {
  theta <- do.call(rbind, lapply(chains, function(ch) ch$theta))
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 1.8, 1), col.axis = ink[["secondary"]],
            col.lab = ink[["primary"]])
  on.exit(par(op))
  for (v in c("lambda", "soc")) {
    x <- exp(theta[, v])
    g <- gamma_shape_scale(prior_cfg[[v]]$mode, prior_cfg[[v]]$sd)
    hist(x, breaks = 40, freq = FALSE, col = check_colours[["blue"]], border = "white",
         main = if (v == "soc") "mu" else "lambda", font.main = 1, cex.main = 0.9,
         xlab = "posterior draws (all chains)", ylab = "density")
    grid_x <- seq(min(x), max(x), length.out = 200)
    lines(grid_x, dgamma(grid_x, shape = g$shape, scale = g$scale), col = check_colours[["orange"]], lwd = 2)
    legend("topright", legend = c("posterior", "hyperprior"), fill = c(check_colours[["blue"]], NA),
           border = c("white", NA), col = c(NA, check_colours[["orange"]]), lwd = c(NA, 2),
           bty = "n", cex = 0.8, text.col = ink[["primary"]])
  }
}

#' A derived quantity by chain: rgdp growth over the forecast horizon.
#'
#' Density of each chain's kept draws (derived_quantities()); chains that
#' have converged give overlapping curves.
#'
#' @param derived Output of derived_quantities().
plot_derived <- function(derived) {
  v <- derived$values
  chains <- split(v$growth, v$chain)
  cols <- colours_for_chains(length(chains))
  dens <- lapply(chains, density)
  op <- par(mar = c(4, 4, 1, 1), col.axis = ink[["secondary"]], col.lab = ink[["primary"]])
  on.exit(par(op))
  plot(NA, xlim = range(unlist(lapply(dens, function(d) d$x))),
       ylim = c(0, max(unlist(lapply(dens, function(d) d$y)))), bty = "n",
       xlab = sprintf("%s growth over %d quarters, no-shock path (log points)", derived$variable,
                      derived$horizon),
       ylab = "density")
  for (k in seq_along(dens)) lines(dens[[k]], col = cols[k], lwd = 2)
  legend("topright", legend = paste("chain", names(chains)), col = cols, lwd = 2, bty = "n",
         cex = 0.8, text.col = ink[["primary"]])
}

#' Save the run's diagnostic figures as PNG files.
#'
#' @param chains Output of run_chains().
#' @param summary Output of summarize_mcmc() (with derived quantities).
#' @param setup Output of estimation_setup().
#' @param dir Folder to write to (created if needed).
#' @return Named character vector of the paths written.
save_mcmc_plots <- function(chains, summary, setup, dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  nm <- colnames(chains[[1]]$theta)
  labels <- hyper_labels(nm, colnames(setup$y))
  psi <- grep("^psi[0-9]+$", nm)
  png_file <- function(name, width, height, f) {
    path <- file.path(dir, paste0(name, ".png"))
    png(path, width = width, height = height, res = 110)
    on.exit(dev.off())
    f()
    path
  }
  c(
    chain_legend = png_file("chain_legend", 900, 60, function() plot_chain_legend(length(chains))),
    trace_main = png_file("trace_main", 1100, 330, function() {
      plot_traces(chains, c("lambda", "soc", "log_post"), c("log lambda", "log mu", "log p(theta | Y)"),
                  n_col = 3)
    }),
    trace_psi = png_file("trace_psi", 1100, 1500, function() {
      plot_traces(chains, nm[psi], paste("log", labels[psi]), n_col = 4)
    }),
    tuning = png_file("tuning", 1100, 380, function() plot_tuning(chains)),
    hyper_position = png_file("hyper_position", 900, 900, function() {
      plot_hyper_position(summary, setup, labels)
    }),
    lambda_mu_prior = png_file("lambda_mu_prior", 1100, 380, function() {
      plot_lambda_mu_prior(chains, setup$prior_cfg)
    }),
    companion = png_file("companion", 900, 380, function() plot_companion(summary$max_modulus)),
    derived = png_file("derived", 900, 380, function() plot_derived(summary$derived))
  )
}

#' Render the MCMC diagnostics report into the run's folder.
#'
#' @param qmd Path of reports/mcmc_diagnostics.qmd.
#' @param run_dir The run's folder (<output_root>/<run_id>/), holding the
#'   run files and plots/.
#' @param depends_on Paths the report reads (passed so that the pipeline
#'   re-renders when they change; not used otherwise).
#' @return Path of the rendered HTML file.
render_mcmc_report <- function(qmd, run_dir, depends_on = NULL) {
  dir_abs <- normalizePath(run_dir, winslash = "/")
  quarto::quarto_render(qmd, execute_params = list(run_dir = dir_abs), quiet = TRUE)
  rendered <- sub("\\.qmd$", ".html", qmd)
  target   <- file.path(run_dir, "mcmc_diagnostics.html")
  if (file.exists(target)) file.remove(target)
  if (!file.rename(rendered, target)) stop("could not move the rendered report to ", target)
  target
}

#' Histogram of the largest companion eigenvalue modulus across draws.
#'
#' @param max_mod Vector from companion_max_modulus(), one value per draw.
plot_companion <- function(max_mod) {
  op <- par(mar = c(4, 4, 1, 1), col.axis = ink[["secondary"]], col.lab = ink[["primary"]])
  on.exit(par(op))
  hist(max_mod, breaks = 30, col = check_colours[["blue"]], border = "white", main = "",
       xlab = "largest eigenvalue modulus of the companion matrix", ylab = "draws")
  abline(v = 1, col = check_colours[["orange"]], lwd = 2)
}

#' Display names: lambda, mu, and psi with the variable's short name.
hyper_labels <- function(names_h, variables) {
  out <- names_h
  out[out == "soc"] <- "mu"
  psi <- grepl("^psi[0-9]+$", out)
  out[psi] <- paste0("psi ", variables[as.integer(sub("psi", "", out[psi]))])
  out
}
