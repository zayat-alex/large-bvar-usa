# check_priors.R: run the prior and hyperparameter checks on a small model and
# write the results and the report.
#
# Run from the project root (so .Rprofile activates renv):
#   Rscript scripts/check_priors.R
#
# What it does (R/estimate_checks.R has the checks, R/estimate_ml.R the
# independent marginal likelihood):
#   1. prior mean placement in BVAR's prior object;
#   2. posterior mean under a very tight prior;
#   3. BVAR's log posterior of the hyperparameters against the independent
#      implementation, at several points, for the model as specified (realized
#      volatility with a white-noise prior mean) and with every variable given
#      a random-walk prior mean;
#   4. hyperparameter modes for three psi hyperpriors x two psi starting-value
#      options, from BVAR and from the textbook formula;
#   5. centring of the sum-of-coefficients dummies;
#   6. names of the estimated hyperparameters; psi starting values.
#
# Small model: 6 of the 31 variables (one real activity, one price, the
# unemployment rate, the 2-year yield, the stock index, realized volatility),
# 5 lags, estimation sample from COVID method none. Takes a few
# minutes. Reads the data build (outputs/data_<vintage>/data.rds); writes
# outputs/checks_priors/ (prior_checks.rds, PNG plots, prior_checks.html).

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

# ---- Setup ----
cfg     <- yaml::read_yaml("config/run.yml")
cfg$covid$method <- "none"                 # these checks use the full sample
series  <- read_series_config("config/series.csv")
lags    <- cfg$model$lags
seed    <- cfg$mcmc$seed
out_dir <- "outputs/checks_priors"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

small_vars <- c("rgdp", "corepce", "urate", "gs2", "sp500", "rvol")
ser_small  <- series[match(small_vars, series$short_name), ]
model_data <- readRDS(file.path(cfg$data$output_dir, "data.rds"))$model_data
covid      <- covid_treatment(model_data, cfg)
y          <- est_data_to_matrix(covid$est_data[, c("quarter", small_vars)], ser_small)

# The three psi hyperprior options (shape, scale of the inverse gamma):
#   glp          : Giannone, Lenza and Primiceri (2015) as written, (0.02)^2,
#                  calibrated by them for data in 4 * ln units;
#   glp_rescaled : the same belief expressed for 100 * ln units (scale x 625);
#   bvar_default : the BVAR package default.
hyperpriors <- data.frame(option = c("glp", "glp_rescaled", "bvar_default"),
                          shape  = c(0.0004, 0.0004, 0.004),
                          scale  = c(0.0004, 0.25, 0.004))
psi_starts <- c("ar_variance", "ar_sd_auto")

results <- list(
  run = list(time = Sys.time(), variables = small_vars, lags = lags,
             sample = paste(rownames(y)[1], "to", rownames(y)[nrow(y)], "(", nrow(y), "quarters )"),
             covid_method = cfg$covid$method, prior_cfg = cfg$prior, hyperpriors = hyperpriors,
             bvar_version = as.character(packageVersion("BVAR")), seed = seed)
)

# ---- Base fit: the configured prior ----
built <- build_priors(cfg$prior, ser_small, y, lags)
fit   <- fit_for_mode(y, lags, built$priors, seed)
results$base <- list(b = built$b, psi_table = built$psi_table, optim = fit$optim,
                     hyper_names = colnames(fit$hyper))

# ---- Check 1: prior mean placement ----
results$check1_prior_mean <- check_prior_mean(fit, built$b)

# ---- Check 2: very tight prior ----
results$check2_tight <- check_tight_prior(fit, built$b, mu = cfg$prior$soc$mode,
                                          psi = built$psi_table$start)
results$check2_sequence <- check_tight_sequence(fit, built$b, mu = cfg$prior$soc$mode,
                                                psi = built$psi_table$start)

# ---- Check 3: marginal likelihood cross-check ----
psi0 <- built$psi_table$start
points <- data.frame(label = c("BVAR mode", "prior modes", "mu = 0.01", "mu = 0.1", "mu = 10",
                               "lambda = 1"),
                     lambda = c(fit$optim$par[["lambda"]], 0.2, 0.2, 0.2, 0.2, 1),
                     mu = c(fit$optim$par[["soc"]], 1, 0.01, 0.1, 10, 1))
points$psi <- c(list(unname(fit$optim$par[psi_index(fit$optim$par)])), rep(list(psi0), 5))
results$check3_as_specified <- check_ml_crosscheck(fit, y, built$b, cfg$prior, points)

# Same, with every variable given a random-walk prior mean (rvol as RW).
ser_rw <- ser_small
ser_rw$prior[ser_rw$short_name == "rvol"] <- "RW"
built_rw <- build_priors(cfg$prior, ser_rw, y, lags)
fit_rw   <- fit_for_mode(y, lags, built_rw$priors, seed)
points_rw <- points
points_rw$lambda[1] <- fit_rw$optim$par[["lambda"]]
points_rw$mu[1]     <- fit_rw$optim$par[["soc"]]
points_rw$psi[[1]]  <- unname(fit_rw$optim$par[psi_index(fit_rw$optim$par)])
results$check3_all_rw <- check_ml_crosscheck(fit_rw, y, built_rw$b, cfg$prior, points_rw)
results$check3_all_rw_optim <- fit_rw$optim

# Profile along mu at BVAR's mode (model as specified).
mu_grid <- 10^seq(log10(cfg$prior$soc$min), log10(cfg$prior$soc$max), length.out = 60)
results$check3_mu_profile <- profile_mu(fit, y, built$b, cfg$prior, mu_grid,
                                        lambda = fit$optim$par[["lambda"]],
                                        psi = unname(fit$optim$par[psi_index(fit$optim$par)]))

# ---- Check 4: psi hyperprior x starting values ----
results$check4_grid <- check_psi_grid(y, ser_small, cfg$prior, lags, hyperpriors, psi_starts, seed)

# ---- Check 5: centring of the sum-of-coefficients dummies ----
results$check5_soc_centring <- check_soc_centring(y, ser_small, cfg$prior, lags)

# ---- Check 6: psi starting values under both options ----
results$check6_psi_starts <- merge(
  psi_start_table(y, lags, "ar_variance", cfg$prior$psi$bound_factor)[, c("variable", "ar_variance", "ar_sd")],
  setNames(psi_start_table(y, lags, "ar_sd_auto", 100)[, c("variable", "start", "min", "max")],
           c("variable", "auto_start", "auto_min", "auto_max")),
  by = "variable", sort = FALSE)

# ---- Save results and plots ----
saveRDS(results, file.path(out_dir, "prior_checks.rds"))

png(file.path(out_dir, "psi_hyperpriors.png"), width = 1600, height = 1000, res = 200)
plot_psi_hyperpriors(hyperpriors, built$psi_table$ar_variance)
dev.off()
png(file.path(out_dir, "mu_profile.png"), width = 1600, height = 1000, res = 200)
plot_mu_profile(results$check3_mu_profile, fit$optim$par[["soc"]])
dev.off()
png(file.path(out_dir, "psi_modes.png"), width = 2200, height = 1100, res = 200)
plot_psi_modes(results$check4_grid$psi, hyperpriors)
dev.off()

# ---- Render the report ----
out_abs <- normalizePath(out_dir, winslash = "/")
quarto::quarto_render("reports/prior_checks.qmd", execute_params = list(output_dir = out_abs),
                      quiet = TRUE)
target <- file.path(out_dir, "prior_checks.html")
if (file.exists(target)) file.remove(target)
if (!file.rename("reports/prior_checks.html", target)) stop("could not move the report to ", target)
cat("Wrote", target, "\n")
