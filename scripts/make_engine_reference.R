# make_engine_reference.R: write the stored reference outputs of the engine
# (tests/testthat/fixtures/engine_reference.rds).
#
# Run from the project root (so .Rprofile activates renv):
#   Rscript scripts/make_engine_reference.R
#
# Why it exists: the engine's results depend on many pieces (the
# conditioning algebra, the random-number stream, the quantiles, the IRF
# formulas). The reference stores the inputs and outputs of a small known
# VAR with a fixed seed; tests/testthat/test-engine-run.R recomputes them
# and fails if any number moves by more than 1e-8, so a later change to the
# code cannot change results without notice. Rerun this script only when a
# change of results is intended, and say so in the commit message.
#
# The model: a VAR(2) with three variables (a, b, c), three posterior
# "draws" (the same coefficients scaled down, the covariance scaled up),
# horizon 6, two paths per draw, a scenario with a single fixed cell, two
# cells in one quarter and a skipped quarter, one generalized IRF and one
# recursive IRF; and the same scenario with two soft targets added (a mean
# of a level and a median of a 4-quarter growth rate), tilted on 50 paths
# per draw.

for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)

B1 <- matrix(c(0.6, 0.1, 0.0, 0.2, 0.5, 0.1, 0.0, 0.1, 0.7), 3)   # (B_1)_ij
B2 <- matrix(c(0.2, 0.0, 0.05, -0.1, 0.2, 0.0, 0.0, 0.0, 0.1), 3)
B  <- rbind(c(0.5, -0.2, 0.1), t(B1), t(B2))                       # K x n
Sigma <- matrix(c(1, 0.4, 0.2, 0.4, 0.5, 0.1, 0.2, 0.1, 0.8), 3)
variables <- c("a", "b", "c")
n_draws <- 3
beta  <- array(NA_real_, c(n_draws, dim(B)))
sigma <- array(NA_real_, c(n_draws, dim(Sigma)))
for (i in seq_len(n_draws)) {
  Bi <- B
  Bi[-1, ] <- Bi[-1, ] * (1 - 0.05 * (i - 1))
  beta[i, , ] <- Bi
  sigma[i, , ] <- Sigma * (1 + 0.1 * (i - 1))
}
history <- matrix(c(1.0, 0.2, -0.5, 0.8, 1.1, 0.3, 0.0, -0.4, 0.6, 1.2, 0.9, 0.1, 0.4, 0.7, -0.2),
                  5, 3, dimnames = list(c("2025Q2", "2025Q3", "2025Q4", "2026Q1", "2026Q2"), variables))
horizon <- 6

inputs <- list(beta = beta, sigma = sigma, history = history, lags = 2L, variables = variables,
               last_observed = "2026Q2", quarters = forecast_quarters("2026Q2", horizon))
settings <- list(horizon = horizon, method = "gaussian", seed = 20261016L, seed_offset = 0L,
                 paths_per_draw = 2L, quantiles = c(0.05, 0.16, 0.25, 0.5, 0.75, 0.84, 0.95),
                 tilting = list(paths_per_draw = 50L, tol = 1e-10, maxit = 100L, ess_warn_share = 0.1),
                 irf = list(horizon = 8, girf_shocks = "a",
                            recursive = list(slow = "b", policy = "a", fast = "c")))
conditions <- matrix(NA_real_, horizon, 3, dimnames = list(inputs$quarters, variables))
conditions[1, "a"] <- 2
conditions[3, c("b", "c")] <- c(0, 1)

drop_time <- function(x) { x$seconds <- NULL; x }
baseline <- drop_time(engine_baseline(inputs, settings))
scenario <- drop_time(engine_scenario(conditions, inputs, settings, baseline$mean))
irf <- drop_time(engine_irfs(inputs, settings))

# Tilted scenario: the hard cells above plus soft targets on cells that are
# not fixed: the mean of a in quarter 4 and the median of b's 4-quarter
# growth in quarter 6 (base quarter 2, a forecast quarter). The targets are
# set away from the untilted values so that the weights must move.
targets <- data.frame(variable = c("a", "b"), quarter = inputs$quarters[c(4, 6)], h = c(4L, 6L),
                      measure = c("level", "yoy"), statistic = c("mean", "median"),
                      target = c(1.2, 0.3), stringsAsFactors = FALSE)
tilted <- drop_time(engine_scenario(conditions, inputs, settings, baseline$mean, targets))

dir.create("tests/testthat/fixtures", showWarnings = FALSE, recursive = TRUE)
path <- "tests/testthat/fixtures/engine_reference.rds"
saveRDS(list(inputs = inputs, settings = settings, conditions = conditions,
             baseline = baseline, scenario = scenario, irf = irf,
             targets = targets, tilted = tilted,
             created = format(Sys.time(), "%Y-%m-%d"), r_version = R.version.string),
        path)
cat("Wrote", path, "\n")
