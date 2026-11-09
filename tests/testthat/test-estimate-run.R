# test-estimate-run.R: the estimation stage as the pipeline runs it
# (R/estimate_run.R): settings, chains, pooled draws, summary, run id, files.
#
# A small simulated VAR (3 variables, 2 lags), so the tests are fast and need
# no data build. The full-size run is scripts/estimate.R.

set.seed(321)
run_y <- cbind(
  a = 100 + cumsum(rnorm(80, 0.5, 1)),
  b = 50 + cumsum(rnorm(80, 0.2, 0.5)),
  c = 20 + rnorm(80, 0, 3)
)
run_series <- data.frame(short_name = c("a", "b", "c"), prior = c("RW", "RW", "WN"),
                         stringsAsFactors = FALSE)
run_cfg <- list(
  data = list(qd_vintage_file = "fred_qd_vintages/FRED-QD_2026m09.csv", api_realtime_date = "2026-09-30"),
  sample = list(start = "1973Q1", end = "latest"),
  model = list(lags = 2),
  prior = list(
    lambda = list(mode = 0.2, sd = 0.4, min = 0.0001, max = 5), alpha = 2,
    soc = list(mode = 1, sd = 1, min = 0.0001, max = 50), soc_ybar = "presample",
    constant_var = 1e7,
    psi = list(shape = 0.0004, scale = 0.0004, start = "ar_variance", bound_factor = 100)),
  covid = list(methods = c("none", "precovid_params"), window = c("2020Q1", "2021Q4"),
               cutoff = "2019Q4"),
  mcmc = list(profile = "smoke", seed = 100, n_cores = 1,
              mode_search = list(start_factors = 1, maxit = 2000),
              chain_start = list(type = "mode"),
              convergence = list(ess_min = 100, rhat_max = 1.05),
              tuning = list(target = c(0.20, 0.35), adapt_share = 0.5, window = 50),
              profiles = list(
                smoke = list(n_draw = 400, n_burn = 200, n_thin = 1, n_chains = 2, n_engine_draws = 20),
                full = list(n_draw = NULL, n_burn = NULL, n_thin = NULL, n_chains = NULL,
                            n_engine_draws = NULL))),
  output_root = "outputs")

run_setup    <- estimation_setup(run_y, run_series, run_cfg$prior, lags = 2)
run_settings <- mcmc_settings(run_cfg)
run_mode     <- find_hyper_mode(run_setup, start_factors = 1)
run_prop     <- hyper_proposal(run_mode$theta_mode, run_setup)
run_ch       <- run_chains(run_setup, run_mode$theta_mode, run_prop$W, run_settings, n_cores = 1)

test_that("mcmc_settings applies the seed rule and rejects missing or inconsistent counts", {
  expect_equal(run_settings$chain_seeds, c(101, 102))
  expect_equal(run_settings$coefficient_seed, 1100)
  expect_equal(run_settings$n_keep_per_chain, 200)

  cfg <- run_cfg; cfg$mcmc$profile <- "full"
  expect_error(mcmc_settings(cfg), "has no value for: n_draw")
  cfg <- run_cfg; cfg$mcmc$profile <- "long"
  expect_error(mcmc_settings(cfg), "must be one of")
  cfg <- run_cfg; cfg$mcmc$profiles$smoke$n_burn <- 400
  expect_error(mcmc_settings(cfg), "n_burn must be smaller")
  cfg <- run_cfg; cfg$mcmc$profiles$smoke$n_engine_draws <- 401
  expect_error(mcmc_settings(cfg), "exceeds the retained draws")
  cfg <- run_cfg; cfg$mcmc$profiles$smoke$n_thin <- 1.5
  expect_error(mcmc_settings(cfg), "positive whole number")
})

test_that("chains are reproducible, differ by seed, and pass the tuning settings on", {
  expect_length(run_ch, 2)
  again <- run_chains(run_setup, run_mode$theta_mode, run_prop$W, run_settings, n_cores = 1)
  expect_identical(again[[1]]$theta, run_ch[[1]]$theta)
  expect_false(identical(run_ch[[1]]$theta, run_ch[[2]]$theta))
  expect_equal(run_ch[[2]]$settings$seed, 102)
  expect_equal(run_ch[[1]]$settings$window, 50)
  expect_equal(nrow(run_ch[[1]]$theta), 200)
  # The chain equals a direct call of run_metropolis() with the same settings.
  direct <- run_metropolis(function(t) log_post_theta(t, run_setup), run_mode$theta_mode, run_prop$W,
                           400, 200, 1, 101, target = c(0.20, 0.35), adapt_share = 0.5, window = 50)
  expect_identical(direct$theta, run_ch[[1]]$theta)
})

test_that("parallel chains equal sequential chains", {
  par_ch <- run_chains(run_setup, run_mode$theta_mode, run_prop$W, run_settings, n_cores = 2,
                       r_dir = file.path(project_root(), "R"))
  for (k in 1:2) {
    expect_identical(par_ch[[k]]$theta, run_ch[[k]]$theta)
    expect_identical(par_ch[[k]]$acceptance, run_ch[[k]]$acceptance)
  }
})

test_that("mcmc_settings records the start seeds and checks the chain-start settings", {
  expect_equal(run_settings$start_seeds, c(2101, 2102))
  expect_equal(run_settings$chain_start, list(type = "mode", scale = 0))
  cfg <- run_cfg; cfg$mcmc$chain_start <- list(type = "random")
  expect_error(mcmc_settings(cfg), "must be 'mode' or 'overdispersed'")
  cfg <- run_cfg; cfg$mcmc$chain_start <- list(type = "overdispersed")
  expect_error(mcmc_settings(cfg), "scale must be a positive number")
  cfg$mcmc$chain_start$scale <- 2
  expect_equal(mcmc_settings(cfg)$chain_start, list(type = "overdispersed", scale = 2))
})

test_that("chain_start: the mode, or mode + scale L z, reproducible and inside the bounds", {
  th <- run_mode$theta_mode; W <- run_prop$W
  expect_identical(chain_start(th, W, run_setup, "mode", 0, 1)$theta, th)
  a <- chain_start(th, W, run_setup, "overdispersed", 2, 2101)
  expect_identical(chain_start(th, W, run_setup, "overdispersed", 2, 2101), a)
  expect_false(identical(chain_start(th, W, run_setup, "overdispersed", 2, 2102)$theta, a$theta))
  expect_length(a$clipped, 0)
  # Unclipped, the start is exactly mode + 2 L z with z the seed's normals.
  set.seed(2101); z <- rnorm(length(th))
  expect_equal(unname(a$theta), unname(th + 2 * as.vector(t(chol(W)) %*% z)))
  expect_equal(names(a$theta), names(th))
  # A huge scale pushes coordinates outside: they are clipped 1% inside the bounds.
  b <- chain_start(th, W, run_setup, "overdispersed", 1e4, 2101)
  expect_true(length(b$clipped) > 0)
  expect_true(all(b$theta >= log(run_setup$lower * 1.01) - 1e-12 &
                  b$theta <= log(run_setup$upper / 1.01) + 1e-12))
})

test_that("overdispersed chains start where chain_start says, and parallel equals sequential", {
  cfg <- run_cfg; cfg$mcmc$chain_start <- list(type = "overdispersed", scale = 2)
  st <- mcmc_settings(cfg)
  seq_ch <- run_chains(run_setup, run_mode$theta_mode, run_prop$W, st, n_cores = 1)
  for (k in 1:2) {
    expect_identical(seq_ch[[k]]$start,
                     chain_start(run_mode$theta_mode, run_prop$W, run_setup, "overdispersed", 2,
                                 st$start_seeds[k])$theta)
  }
  expect_false(identical(seq_ch[[1]]$start, seq_ch[[2]]$start))
  par_ch <- run_chains(run_setup, run_mode$theta_mode, run_prop$W, st, n_cores = 2,
                       r_dir = file.path(project_root(), "R"))
  for (k in 1:2) {
    expect_identical(par_ch[[k]]$theta, seq_ch[[k]]$theta)
    expect_identical(par_ch[[k]]$burn_theta, seq_ch[[k]]$burn_theta)
  }
})

test_that("run_metropolis keeps the burn-in path, thinned like the retained draws", {
  ch <- run_ch[[1]]
  expect_equal(dim(ch$burn_theta), c(200, ncol(ch$theta)))
  expect_length(ch$burn_log_target, 200)
  expect_true(all(is.finite(ch$burn_log_target)))
  thin <- run_metropolis(function(t) log_post_theta(t, run_setup), run_mode$theta_mode, run_prop$W,
                         400, 200, 3, 101, window = 50)
  expect_equal(nrow(thin$burn_theta), floor(200 / 3))
  # The burn-in record does not change the retained draws (same seed, same sampler).
  plain <- run_metropolis(function(t) log_post_theta(t, run_setup), run_mode$theta_mode, run_prop$W,
                          400, 200, 1, 101, window = 50)
  expect_identical(plain$theta, ch$theta)
  # The last burn-in point and the retained draws are consecutive states of one chain:
  # each log posterior value matches its theta.
  expect_equal(plain$burn_log_target[200], log_post_theta(plain$burn_theta[200, ], run_setup))
})

test_that("check_convergence flags low ESS and chains that disagree", {
  # Two chains of independent normals: high ESS, R-hat about 1.
  fake_chain <- function(k, shift) {
    set.seed(k)
    th <- matrix(rnorm(2000 * 2, mean = shift), ncol = 2, dimnames = list(NULL, c("lambda", "soc")))
    list(theta = th, log_target = rnorm(2000), chain = k, acceptance = 0.3, scale = 1, seconds = 1,
         settings = list(n_draw = 2500, n_burn = 500, seed = k))
  }
  setup2 <- run_setup
  setup2$lower <- c(lambda = 1e-10, soc = 1e-10); setup2$upper <- c(lambda = 1e10, soc = 1e10)
  dr <- pool_and_draw(run_ch, run_setup, run_settings)
  conv <- list(ess_min = 1000, rhat_max = 1.05)
  good <- summarize_mcmc(list(fake_chain(1, 0), fake_chain(2, 0)), dr, setup2, convergence = conv)
  expect_true(all(good$convergence$pass))
  expect_equal(good$hyper$ess_per_1000[1], good$hyper$ess[1] / 4000 * 1000)
  # Chain 2 shifted by 3 standard deviations: R-hat far above 1.05.
  bad <- summarize_mcmc(list(fake_chain(1, 0), fake_chain(2, 3)), dr, setup2, convergence = conv)
  expect_false(bad$convergence$pass[bad$convergence$check == "R-hat (point), maximum"])
  expect_true(bad$rhat_multivariate > 1.05)
  # Too few draws for the ESS target.
  strict <- check_convergence(good, list(ess_min = 1e6, rhat_max = 1.05))
  expect_false(strict$pass[1])
})

test_that("derived_quantities forecasts each kept draw's no-shock path from the history", {
  dr <- pool_and_draw(run_ch, run_setup, run_settings)
  mm <- vapply(seq_len(20), function(i) companion_max_modulus(dr$beta[i, , ], 2), numeric(1))
  d <- derived_quantities(dr, run_y, mm, variable = "a", horizon = 4)
  # By hand for draw 1: iterate y_t = B' x_t four times from the last 2 rows.
  B <- dr$beta[1, , ]; lagged <- run_y[80:79, ]
  for (h in 1:4) {
    y <- as.vector(crossprod(B, c(1, as.vector(t(lagged)))))
    lagged <- rbind(y, lagged[1, ])
  }
  expect_equal(d$values$growth[1], unname(y[1] - run_y[80, "a"]))
  expect_equal(d$values$max_modulus, mm)
  expect_equal(nrow(d$by_chain), 4)
  expect_true(all(is.finite(d$rhat)))
})

test_that("draw_posterior_coefficients reproduces the per-draw loop it replaced", {
  hyper <- t(apply(exp(run_ch[[1]]$theta[1:5, ]), 1, complete_hyper, setup = run_setup))
  got <- draw_posterior_coefficients(hyper, run_setup, seed = 7)
  # The loop as it was written in estimate_bvar().
  set.seed(7)
  for (i in 1:5) {
    d <- draw_coefficients(niw_posterior(hyper[i, ], run_setup), 1)
    expect_identical(got$beta[i, , ], d$beta[1, , ])
    expect_identical(got$sigma[i, , ], d$sigma[1, , ])
  }
})

test_that("pool_and_draw keeps evenly spaced rows of the stacked chains", {
  dr <- pool_and_draw(run_ch, run_setup, run_settings)
  pick <- round(seq(1, 400, length.out = 20))
  expect_equal(dr$chain, ifelse(pick <= 200, 1, 2))
  expect_equal(dr$draw, ifelse(pick <= 200, pick, pick - 200))
  stacked <- rbind(run_ch[[1]]$theta, run_ch[[2]]$theta)
  expect_equal(unname(dr$hyper), unname(exp(stacked[pick, ])))
  expect_equal(dr$log_post, c(run_ch[[1]]$log_target, run_ch[[2]]$log_target)[pick] -
                 rowSums(stacked[pick, ]))
  expect_equal(dim(dr$beta), c(20, 1 + 3 * 2, 3))
  expect_equal(dim(dr$sigma), c(20, 3, 3))
  # Same seed, same draws.
  expect_identical(pool_and_draw(run_ch, run_setup, run_settings)$beta, dr$beta)
})

test_that("summarize_mcmc reports every hyperparameter and the kept draws", {
  dr <- pool_and_draw(run_ch, run_setup, run_settings)
  s  <- summarize_mcmc(run_ch, dr, run_setup)
  expect_equal(s$hyper$hyperparameter, c(names(run_setup$start), "log_post"))
  expect_true(all(is.finite(s$hyper$rhat)))
  expect_true(all(s$hyper$ess > 0))
  expect_equal(nrow(s$chains), 2)
  expect_equal(s$sigma_pd, 1)
  expect_length(s$max_modulus, 20)
  # One chain: no R-hat.
  s1 <- summarize_mcmc(run_ch[1], dr, run_setup)
  expect_true(all(is.na(s1$hyper$rhat)))
})

test_that("make_run_id is stable and changes only with settings that change results", {
  a <- make_run_id(run_cfg, "none", run_settings, run_series)
  expect_match(a$run_id, "^2026m09_none_smoke_[0-9a-f]{8}$")
  expect_identical(make_run_id(run_cfg, "none", run_settings, run_series)$run_id, a$run_id)
  # A different method, seed or prior mean gives a different id.
  expect_false(make_run_id(run_cfg, "precovid_params", run_settings, run_series)$hash == a$hash)
  cfg <- run_cfg; cfg$mcmc$seed <- 101
  expect_false(make_run_id(cfg, "none", mcmc_settings(cfg), run_series)$hash == a$hash)
  ser <- run_series; ser$prior[3] <- "RW"
  expect_false(make_run_id(run_cfg, "none", run_settings, ser)$hash == a$hash)
  # The number of cores does not.
  cfg <- run_cfg; cfg$mcmc$n_cores <- 4
  expect_identical(make_run_id(cfg, "none", mcmc_settings(cfg), run_series)$hash, a$hash)
})

test_that("write_run_files writes the run folder and a readable meta.json", {
  q <- seq(as.Date("1973-01-01"), by = "quarter", length.out = 80)
  covid <- list(est_data = data.frame(quarter = q, run_y), origin_data = data.frame(quarter = q, run_y),
                est_weights = NULL, notes = data.frame())
  dr <- pool_and_draw(run_ch, run_setup, run_settings)
  s  <- summarize_mcmc(run_ch, dr, run_setup)
  id <- make_run_id(run_cfg, "none", run_settings, run_series)
  root <- file.path(tempdir(), "runs")
  paths <- write_run_files(id, root, "none", covid, timed_call(function() run_mode),
                           timed_call(function() run_prop), run_ch, dr, s, run_settings, run_cfg)
  expect_true(all(file.exists(paths)))
  expect_equal(normalizePath(dirname(paths[1]), winslash = "/"),
               normalizePath(file.path(root, id$run_id), winslash = "/"))
  meta <- jsonlite::read_json(paths[5])
  expect_equal(meta$run_id, id$run_id)
  expect_equal(meta$data$estimation_quarters, 80)
  expect_equal(unlist(meta$data$estimation_sample), c("1973Q1", "1992Q4"))
  expect_equal(unlist(meta$mcmc$chain_seeds), c(101, 102))
  expect_identical(readRDS(paths[2])$beta, dr$beta)
  unlink(root, recursive = TRUE)
})
