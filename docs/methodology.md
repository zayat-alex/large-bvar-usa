# Methodology

This document sets out how the model is specified, estimated and used for forecasting: the equations, the priors and hyperparameters, the sampling algorithm, its settings and diagnostics, and how the posterior draws are used: baseline forecasts, scenarios (conditional forecasts, with hard conditions and with soft conditions imposed by entropic tilting) and impulse responses. It follows Crump, Eusepi, Giannone, Qian and Sbordone, "A Large Bayesian VAR of the United States Economy" (Federal Reserve Bank of New York Staff Report 976, 2021; "the paper"), which in turn follows Giannone, Lenza and Primiceri (2015; "GLP"). Where the paper is silent, the choice made here and its source are stated. Section 12 lists every difference from the paper.

The data are defined in `docs/data.md`. The code is in `R/` and runs as stages 2 to 5 of the `targets` pipeline (`_targets.R`: COVID treatment, estimation, forecasting engine, export to the website); all settings are in `config/run.yml`.

**Data and the numbers quoted.** The model was developed, and every check and comparison in this document was run, on data through 2026:Q2 (FRED-QD vintage 2026m09). The published forecasts have since been updated to data through 2026:Q3 (FRED-QD vintage 2026m10, other sources as of 9 November 2026): the published model (`exclude_window`) was re-estimated and its forecasts, scenarios and impulse responses rebuilt, with the same settings and seeds. The checks, the comparisons of COVID treatments and the numbers quoted in this document (hyperparameters, roots, drift, scenario effects, check results) were not rerun for the update; they describe the 2026:Q2 runs.

**Equation labels.** Equations are labelled (M1), (M2), .... Code comments that implement an equation cite its label (for example `docs/methodology.md (M14)`), so the code and the formulas can be read side by side.

---

## 1. Notation

| Symbol | Meaning | Dimension / value |
|---|---|---|
| n | number of variables | 31 |
| p | number of lags | 5 |
| T | quarters in the estimation sample | 214 (1973:Q1–2026:Q2) or 188 (1973:Q1–2019:Q4), section 11 |
| N = T − p | quarters used as left-hand side | 209 or 183 |
| K = 1 + n·p | regressors per equation (constant plus p lags of n variables) | 156 |
| y_t | the 31 variables in quarter t, in model units (`docs/data.md`, section 1) | n × 1 |
| c | intercepts | n × 1 |
| B_s | coefficients on lag s; (B_s)_ij is the effect of variable j at lag s on variable i | n × n |
| B | all coefficients stacked, B = [c, B_1, ..., B_p]' | K × n |
| ε_t | one-step-ahead forecast error | n × 1 |
| Σ | covariance matrix of ε_t | n × n |
| λ | overall tightness of the Minnesota prior | scalar |
| α | lag decay of the Minnesota prior variance | 2 |
| ψ_j, Ψ | scale of variable j; Ψ = diag(ψ_1, ..., ψ_n) | n values |
| µ | tightness of the sum-of-coefficients prior | scalar |
| d | degrees of freedom of the inverse-Wishart prior | n + 2 = 33 |
| b | prior mean of B | K × n |
| Ω | prior covariance factor of B (diagonal) | K × K |
| h | hyperparameters treated as unknown, h = (λ, µ, ψ_1, ..., ψ_n) | 33 |
| θ = log h | the hyperparameters on the log scale, where the sampler moves | 33 |
| ȳ_0 | mean of the first p observations of each variable | n × 1 |
| H | forecast horizon in quarters; T + 1 is the first forecast quarter. As a subscript, h indexes quarters ahead (unrelated to the hyperparameters h) | 20 (forecasts), 40 (impulse responses) |
| L | lower Cholesky factor of Σ, L L' = Σ | n × n |
| Ψ_h | moving-average coefficients: response of y_{t+h} to the forecast error ε_t (Ψ_0 = I); not to be confused with Ψ = diag(ψ) | n × n |
| J | number of simulated forecast paths | draws × paths per draw |
| w_j | weight of path j after entropic tilting | Σ_j w_j = 1 |

---

## 2. Data entering the model

The model is estimated on the 31 series of `docs/data.md`, in levels: `100·ln(x)` for the 22 log series and the raw value for the 9 rates and ratios. Growth rates are computed from forecasts afterwards, not modelled (paper, Section 2.1). The sample starts in 1973:Q1; the first p = 5 quarters serve only as initial lags, so the first quarter on the left-hand side is 1974:Q2.

Units matter in two places: the inverse-gamma hyperprior on ψ (section 4.4) is not scale-free, and scenario paths and shock sizes are stated in model units (1 unit = 1 log point, about 1 percent, for log series; 1 percentage point for rates).

---

## 3. The VAR

One multivariate VAR(5) (paper, eq. 1):

    y_t = c + B_1 y_{t-1} + ... + B_p y_{t-p} + ε_t,    ε_t ~ N(0, Σ) i.i.d.        (M1)

Every equation contains the lags of all 31 variables, and Σ is a full matrix. Unrestricted, this is 31 × 156 = 4,836 coefficients and 496 distinct covariance terms, estimated from about 210 quarters; the prior (section 4) provides the shrinkage that makes this feasible.

**Stacked form.** With x_t = (1, y_{t-1}', ..., y_{t-p}')' and rows t = p+1, ..., T stacked:

    Y = X B + E,     Y: N × n,  X: N × K,  B: K × n        (M2)

Column 1 of X is the constant; columns 2 to n+1 are lag 1 of variables 1..n; then lag 2; and so on. Hence B[1 + (s−1)n + j, i] = (B_s)_ij.

**Companion form**, used for forecasting and stability checks. With z_t = (y_t', ..., y_{t-p+1}')':

    z_t = c̃ + F z_{t-1} + ẽ_t,   F = [ B_1 B_2 ... B_p ; I_{n(p-1)} 0 ]        (M3)

The VAR is stable if every eigenvalue of F lies inside the unit circle. No stability restriction is imposed (the paper imposes none); section 16 reports what the draws show.

---

## 4. Prior

### 4.1 Family

Normal–inverse-Wishart (paper, eqs. 2–3; GLP, eqs. 2.3–2.4):

    Σ | h      ~ IW(Ψ, d)                                  (M4)
    vec(B) | Σ, h ~ N(vec(b), Σ ⊗ Ω)                       (M5)

with d = n + 2, the smallest value for which the prior mean of Σ exists; that mean is Ψ/(d − n − 1) = Ψ. Conjugacy gives the marginal likelihood in closed form (section 5).

### 4.2 Minnesota prior

Prior mean (paper, eq. 4 and Table 1):

    E[(B_s)_ij] = δ_i if i = j and s = 1;  0 otherwise        (M6)
    δ_i = 1 (random walk) for 30 variables;  δ_i = 0 (white noise) for realized volatility

Prior covariance (paper, eq. 5):

    cov[(B_s)_ij, (B_r)_hm] = λ² (1/s^α) Σ_ih / (ψ_j/(d − n − 1))   if m = j and r = s;  0 otherwise        (M7)

In the Kronecker form (M5), Ω is diagonal:

    Ω[1, 1] = v_c;     Ω[1 + (s−1)n + j, 1 + (s−1)n + j] = λ² / (s^α ψ_j),   s = 1..p, j = 1..n        (M8)

- λ scales all prior variances: as λ → 0 the coefficients are held at the prior mean; as λ → ∞ the prior stops mattering.
- 1/s^α with α = 2 ("1/s²" in the paper) shrinks longer lags harder.
- Σ_ih/ψ_j adjusts for units, so a yield in percent and a log level in log points are treated comparably.
- The intercept has a diffuse prior, v_c = 10⁷ (GLP use a diffuse prior; the paper does not state one).

### 4.3 Sum-of-coefficients prior

Written in error-correction form, Δy_t = c + Π y_{t−1} + Σ_s A_s Δy_{t−s} + ε_t with Π = B_1 + ... + B_p − I, the prior shrinks Π towards zero, a unit root in each equation (paper, Section 2.2; GLP, Section 3). It is implemented with n dummy observations added to (Y, X):

    Y⁺ = diag(ȳ_0) / µ,     X⁺ = [ 0_{n×1}, Y⁺, Y⁺, ..., Y⁺ ]   (Y⁺ repeated p times)        (M9)

where ȳ_0 is the mean of the first p observations of each variable (1973:Q1–1974:Q1), as in GLP. µ → ∞ switches the prior off; µ → 0 imposes exact differencing.

No single-unit-root (dummy-initial-observation) prior is used: the paper lists only the Minnesota and sum-of-coefficients priors, and treats λ, Ψ and µ as the unknown hyperparameters.

The sum-of-coefficients prior pulls the sum of each variable's own-lag coefficients towards 1, including realized volatility, whose Minnesota mean is 0. The paper applies both priors, and so does this model.

### 4.4 Hyperpriors and bounds

The paper treats λ, Ψ and µ as unknown but does not state their hyperpriors. They are taken from GLP (Section 3):

| Hyperparameter | Hyperprior | Search bounds | Starting value of the mode search |
|---|---|---|---|
| λ | Gamma, mode 0.2, s.d. 0.4 | [0.0001, 5] | 0.2 |
| µ | Gamma, mode 1, s.d. 1 | [0.0001, 50] | 1 |
| ψ_j | inverse-gamma, shape 0.0004, scale 0.0004 | start/100 to start × 100 | residual variance of an OLS AR(p) of variable j |
| α | fixed at 2 | — | — |

A Gamma with mode m and standard deviation s has shape k = (2 + r + √(r(r + 4)))/2, r = m²/s², and scale m/(k − 1).

The hyperpriors are truncated to the bounds. GLP state no bounds; the bounds follow the `BVAR` package's defaults for λ and µ and its ×/÷ 100 rule for ψ. In all estimated posteriors no draw lies within 1 percent of a bound (section 8.3).

**Units of the ψ hyperprior.** GLP's inverse-gamma (shape = scale = 0.02²) is calibrated for data in 4·ln(x) units, whereas this model, like the paper, uses 100·ln(x). The values are used as GLP write them. The hyperprior is extremely diffuse either way: the posterior mode of ψ was practically unchanged when the hyperprior was rescaled to 100·ln units or replaced by the `BVAR` default (`reports/prior_checks.qmd`).

---

## 5. Marginal likelihood and the posterior of the hyperparameters

The hyperparameters are estimated from their posterior (GLP, Section 2):

    p(h | Y) ∝ p(Y | h) · p(h)                                        (M10)
    p(Y | h) = ∫ p(Y | B, Σ) p(B, Σ | h) dB dΣ                         (M11)

For the prior (M4)–(M5), with

    B̂ = (X'X + Ω⁻¹)⁻¹ (X'Y + Ω⁻¹ b)                                  (M12)
    Ŝ = (Y − X B̂)'(Y − X B̂) + (B̂ − b)' Ω⁻¹ (B̂ − b)                    (M13)

the marginal likelihood is (GLP, Appendix A)

    log p(Y | h) = −(nN/2) log π + log Γ_n((N + d)/2) − log Γ_n(d/2)
                   − (n/2) log|I + Ω^{1/2} X'X Ω^{1/2}|
                   − (N/2) log|Ψ| − ((N + d)/2) log|I + Ψ^{−1/2} Ŝ Ψ^{−1/2}|        (M14)

where Γ_n is the multivariate gamma function. The determinant |I + Ω^{1/2} X'X Ω^{1/2}| is computed from a symmetric eigendecomposition: its eigenvalues are all at least 1, so the computation stays stable even with the diffuse intercept and a tight sum-of-coefficients prior.

**Dummy observations.** The sum-of-coefficients prior enters as the dummy observations (M9). The marginal likelihood of the actual data is the marginal likelihood of the augmented data divided by that of the dummy observations alone (GLP; Bańbura, Giannone and Lenza 2015, Appendix B):

    p(Y | h) = p([Y⁺; Y] | h) / p(Y⁺ | h)                              (M15)

Both terms are computed with (M12)–(M14), the second with B̂ and Ŝ evaluated for the dummy observations themselves. Section 13 explains why this detail matters.

The log posterior is (M14)–(M15) plus the log hyperprior densities of section 4.4, computed in `log_posterior_direct()` (`R/estimate_ml.R`).

**Interpretation** (GLP, eq. 2.5). The marginal likelihood trades off in-sample fit against a complexity penalty. It amounts to rewarding good one-step-ahead forecasts within the sample, which is why it selects tighter priors when there are many parameters relative to observations.

---

## 6. Posterior of the coefficients given the hyperparameters

Given h, the posterior is again normal–inverse-Wishart (GLP, Appendix A):

    Σ | h, Y      ~ IW(Ψ + Ŝ, N + d)                                  (M16)
    vec(B) | Σ, h, Y ~ N( vec(B̂), Σ ⊗ (X'X + Ω⁻¹)⁻¹ )                  (M17)

where X, Y and N include the n sum-of-coefficients dummy rows. These draws are exact: given h, no Markov chain is involved.

---

## 7. Sampling algorithm

The algorithm is GLP's (Appendix B), which the paper follows. The hyperparameters are drawn from their marginal posterior (the coefficients are integrated out), and Σ and B are drawn exactly given each retained hyperparameter draw. Only the hyperparameter step is a Markov chain.

**Sampling scale.** The sampler moves on θ = log h. The change of variables adds the log Jacobian:

    log p(θ | Y) = log p(h | Y) + Σ_k θ_k,   h = exp(θ)        (M19)

and the density is zero outside the bounds. The posterior of h is unchanged by this choice.

**Steps** (`R/estimate_bvar.R`, `R/estimate_run.R`):

1. **Mode.** Maximize (M19) with L-BFGS-B within the log bounds, from three starts: the starting values of section 4.4 multiplied by 1, 0.5 and 2. All three must reach the same mode. The objective clips its argument to the bounds, because the optimizer's finite-difference gradient can step one rounding unit outside a bound.
2. **Proposal.** W = (−H)⁻¹, with H the numerical Hessian of (M19) at the mode: the covariance of the normal approximation to the posterior.
3. **Starting points.** Each of the chains starts at an overdispersed point:

       θ_0^(k) = θ̂ + s L z^(k),   L L' = W,  z^(k) ~ N(0, I),  s = 2        (M20)

   about two approximate posterior standard deviations from the mode θ̂ in a random direction, moved 1 percent inside a bound if it would fall outside. Chains that start apart and come to agree must each have forgotten their start, which is what makes R-hat (Gelman and Rubin 1992) informative.
4. **Metropolis.** Propose θ* = θ + L_c u, L_c L_c' = c W, u ~ N(0, I); accept with probability min{1, exp(log p(θ* | Y) − log p(θ | Y))}. The scale c starts at 2.38²/33 and is tuned during the first half of the burn-in in windows of 100 iterations (× 0.8 if the window's acceptance is below 20 percent, × 1.25 if above 35 percent), then fixed, so that the retained draws come from a fixed Metropolis kernel.
5. **Retained draws.** After the burn-in, every n_thin-th state is kept. The chains run in parallel; each has its own seeds, so the results do not depend on the number of processor cores.
6. **Draws for the forecasts.** The retained draws of all chains are pooled and 1,000 of them, evenly spaced, are kept. For each, one Σ is drawn from (M16) and one B from (M17). Each triple (h, Σ, B) is a draw from the joint posterior.

The algorithm was checked on a 6-variable version of the model against cases with known answers: the coefficient draws against their exact moments, the sampler against the posterior computed on a grid, and the whole estimator against the `BVAR` package in a configuration where the package is exact (a check not included in this repository; see the README).

---

## 8. Settings, diagnostics and results

### 8.1 Settings

The table lists the settings that determine the estimates, as they are in `config/run.yml`. A test (`tests/testthat/test-docs-methodology.R`) checks that the two agree.

| Setting (`config/run.yml`) | Value | Meaning |
|---|---|---|
| `model.lags` | 5 | lags p |
| `prior.lambda.mode` | 0.2 | Gamma hyperprior on λ: mode |
| `prior.lambda.sd` | 0.4 | Gamma hyperprior on λ: standard deviation |
| `prior.lambda.min` | 0.0001 | lower search bound of λ |
| `prior.lambda.max` | 5 | upper search bound of λ |
| `prior.alpha` | 2 | lag decay α, fixed |
| `prior.soc.mode` | 1 | Gamma hyperprior on µ: mode |
| `prior.soc.sd` | 1 | Gamma hyperprior on µ: standard deviation |
| `prior.soc.min` | 0.0001 | lower search bound of µ |
| `prior.soc.max` | 50 | upper search bound of µ |
| `prior.soc_ybar` | presample | ȳ_0 is the mean of the first p observations (GLP) |
| `prior.constant_var` | 1e7 | prior variance factor of the intercept v_c |
| `prior.psi.shape` | 0.0004 | inverse-gamma hyperprior on ψ_j: shape |
| `prior.psi.scale` | 0.0004 | inverse-gamma hyperprior on ψ_j: scale |
| `prior.psi.start` | ar_variance | starting value: residual variance of an OLS AR(p) |
| `prior.psi.bound_factor` | 100 | search bounds of ψ_j: start/100 to start × 100 |
| `covid.methods` | none, precovid_params, exclude_window | the models estimated (section 11) |
| `covid.window` | 2020Q1, 2021Q4 | first and last quarter whose equations `exclude_window` leaves out |
| `covid.cutoff` | 2019Q4 | last estimation quarter of `precovid_params` |
| `mcmc.seed` | 20261016 | master seed |
| `mcmc.mode_search.start_factors` | 1, 0.5, 2 | starts of the mode search, as multiples of the starting values |
| `mcmc.mode_search.maxit` | 2000 | iteration limit of the mode search |
| `mcmc.chain_start.type` | overdispersed | chain starts (M20) |
| `mcmc.chain_start.scale` | 2 | s in (M20) |
| `mcmc.tuning.target` | 0.2, 0.35 | acceptance band used to tune c |
| `mcmc.tuning.adapt_share` | 0.5 | share of the burn-in during which c is tuned |
| `mcmc.tuning.window` | 100 | iterations per tuning window |
| `mcmc.convergence.ess_min` | 1000 | target: effective sample size of every hyperparameter |
| `mcmc.convergence.rhat_max` | 1.05 | target: R-hat of every hyperparameter |
| `mcmc.profiles.full.n_chains` | 4 | chains |
| `mcmc.profiles.full.n_draw` | 42000 | iterations per chain, including burn-in |
| `mcmc.profiles.full.n_burn` | 2000 | burn-in iterations per chain (discarded) |
| `mcmc.profiles.full.n_thin` | 10 | every 10th iteration after the burn-in is retained |
| `mcmc.profiles.full.n_engine_draws` | 1000 | draws of (h, Σ, B) passed to the forecasts |

**Seeds.** Chain k uses seed + k for its proposals and seed + 2000 + k for its starting point; the draws of Σ and B use seed + 1000. Every run records its settings, seeds, code version, R and package versions in `meta.json` in its output folder.

**How the chain length was chosen.** The paper does not report its draw counts (for reference, Bańbura, Giannone and Lenza 2015 use 25,000 draws with 5,000 discarded). The counts here were set from a pilot run of 4 chains × 12,000 iterations: the chains reached the posterior's level within about 500 iterations of their overdispersed starts, and the slowest-mixing quantity yielded about 8.4 effective draws per 1,000 iterations per chain. An effective sample size of 1,000 with a 30 percent margin needs about 155,000 retained iterations, hence 4 chains × 40,000 after a burn-in of 2,000. Thinning by 10 is well below the autocorrelation time (about 110 iterations per independent draw), so it loses little information.

### 8.2 Diagnostics

Computed for every run and reported in `mcmc_diagnostics.html` in its output folder (`R/estimate_run.R`, `R/mcmc_diagnostics.R`, `reports/mcmc_diagnostics.qmd`). They are computed on θ, for every hyperparameter and for the log posterior (M19), with the `coda` package:

- **Effective sample size (ESS):** the number of independent draws the autocorrelated chains are worth, summed over chains. Target ≥ 1,000.
- **R-hat** (Gelman and Rubin 1992): the variance between chains relative to the variance within them; near 1 when the chains agree. Target < 1.05. The multivariate version (Brooks and Gelman 1998) is reported as well.
- **Geweke z** (Geweke 1992): the mean of the first 10 percent of a chain against the mean of its last 50 percent; the largest |z| over chains is reported.
- **Trace plots** of every hyperparameter and of the log posterior, burn-in included; acceptance rate and tuning path per chain.
- **Bounds and hyperpriors:** each posterior within its search bounds; λ and µ against their hyperpriors; share of draws within 1 percent of a bound.
- **Kept coefficient draws:** whether each Σ is positive definite; the largest companion-matrix eigenvalue modulus (M3) of each draw; and, compared across chains, real GDP growth over the 4 quarters after the last observed quarter on each draw's no-shock forecast path. These check that the conclusions also hold for quantities the forecasts use, not only for the hyperparameters.

The targets are used to judge a run, not to tune it: a run that misses them is reported as such.

### 8.3 Results

Estimation runs on the 2026m09 data vintage (QD vintage of September 2026; series through 2026:Q2), 5 lags, 4 chains × 42,000 iterations (burn-in 2,000, thinning 10: 16,000 retained draws), 1,000 draws passed to the forecasts. Intervals are 90 percent posterior (credible) intervals.

| | `none` (1973:Q1–2026:Q2) | `precovid_params` (1973:Q1–2019:Q4) | `exclude_window` (1973:Q1–2026:Q2 without 2020:Q1–2021:Q4) |
|---|---|---|---|
| Acceptance rate (by chain) | 0.296–0.300 | 0.243–0.248 | 0.197–0.297 |
| Smallest ESS (sum over chains) | 1,440 (ψ of baa) | 1,326 (log posterior) | 1,391 (log posterior) |
| Largest R-hat | 1.005 | 1.005 | 1.004 |
| Multivariate R-hat | 1.012 | 1.021 | 1.015 |
| Quantities with \|Geweke z\| > 2, of 34 (largest) | 6 (3.15) | 7 (3.66) | 7 (3.29) |
| Draws within 1% of a bound | 0 | 0 | 0 |
| λ: median (90% interval) | 0.394 (0.356–0.433) | 0.409 (0.368–0.453) | 0.391 (0.353–0.431) |
| µ: median (90% interval) | 0.278 (0.228–0.338) | 0.352 (0.289–0.430) | 0.297 (0.244–0.360) |
| Σ positive definite | every draw | every draw | every draw |
| Largest companion modulus: median, max | 1.011, 1.048 | 1.013, 1.048 | 1.012, 1.042 |
| R-hat of the derived quantities (largest modulus; GDP growth) | 1.0006; 1.0006 | 1.002; 1.002 | 1.004; 1.000 |

All three runs meet both targets (one `exclude_window` chain accepted 19.7 percent of proposals, just below the tuning band of 20–35 percent; its draws enter the ESS and R-hat above). The trace plots show flat, overlapping chains for all 33 hyperparameters. With 34 quantities and the largest |z| of four chains, about six values of |Geweke z| above 2 are expected by chance. For comparison, the posterior mode is λ = 0.401, µ = 0.279 (`none`), λ = 0.416, µ = 0.350 (`precovid_params`) and λ = 0.397, µ = 0.296 (`exclude_window`). In the 6-variable version of the model the posterior favours looser priors (λ ≈ 0.6), as the marginal likelihood selects tighter priors when there are more parameters per observation.

To reproduce: `Rscript scripts/estimate.R` with `mcmc.profile: full`. The run takes about 75 minutes on a laptop (both models, four chains each in parallel), plus about 25 minutes for the two mode searches and Hessians when they are not cached.

---

## 9. Forecasts

### 9.1 Baseline (unconditional) forecast

For each of the 1,000 draws m of (B, Σ), starting from the last p observed quarters:

    y_{T+h}^(m) = c^(m) + Σ_{s=1}^{p} B_s^(m) y_{T+h−s}^(m) + L^(m) u_{T+h},   L^(m) L^(m)' = Σ^(m),  u ~ N(0, I_n)        (M18)

for h = 1, ..., H. Each simulated shock feeds into later quarters through the lags, so the bands widen with the horizon. The collection of paths is a sample from the predictive density (paper, Section 2.3), which combines shock, coefficient and hyperparameter uncertainty. Setting u = 0 gives each draw's forecast mean (the no-shock path). Implemented in `forecast_unconditional()` (`R/forecast.R`); checked against the exact forecast moments of a known VAR(2) (20,000 paths: largest |z| of the means 1.31, variance ratios 0.989–1.016; `reports/forecast_checks.qmd`).

The forecasts start from the latest observed data under both COVID treatments (section 11).

### 9.2 Scenarios (conditional forecasts)

A scenario fixes the values of some variables in some future quarters (hard conditions) and asks for the predictive density of everything else given them: p(y_{T+1:T+H} | y_{1:T}, C), where C is the set of fixed cells (paper, Section 2.3, following Bańbura, Giannone and Lenza 2015). The model fills in every free cell, including the conditioned variables in the quarters that are not fixed.

A scenario is a conditional forecast, not a causal effect. The model explains the fixed values by the combination of forecast errors that is most likely given Σ and the dynamics; it does not say which economic shock produced them. Two scenarios with the same fixed values but different stories behind them give the same answer.

**Stacked form.** For one draw of (B, Σ), stack the future quarters, y = (y_{T+1}', ..., y_{T+H}')' (nH × 1). It is linear in the future shocks:

    y = m + A u,   u ~ N(0, I_{nH}),   A = Θ (I_H ⊗ L)        (M21)

where m is the no-shock path ((M18) with u = 0) and Θ is block lower triangular with block (h, k) equal to Ψ_{h−k}, the moving-average coefficients

    Ψ_0 = I,   Ψ_j = Σ_{s=1}^{min(j, p)} B_s Ψ_{j−s}.

**Gaussian conditioning** (the Waggoner and Zha 1999 form, written on observables). The fixed cells are a selection R y = r. Then G u = d with G = R A and d = r − R m, and

    E[y | R y = r] = m + A G'(GG')⁻¹ d,   Var[y | R y = r] = A (I − P) A',   P = G'(GG')⁻¹ G        (M22)
    u = G'(GG')⁻¹ d + (I − P) z,   z ~ N(0, I),   y = m + A u                                    (M23)

(M23) is an exact draw from the conditional distribution. A is block lower triangular with L on its diagonal, so it is invertible and G has full row rank for any set of distinct cells: every scenario that fixes each cell at most once is well posed. Conditions on closely related variables (for example CPI and core CPI in the same quarter) can still make GG' badly conditioned numerically.

**Kalman simulation smoother** (the paper's algorithm). The VAR in state-space form (paper, eqs. 6–8), with the constant in the state:

    y_t = C_t S_t,   S_{t+1} = G̃ S_t + w_t,   w_t ~ N(0, Q),   S_t = (y_t', ..., y_{t−p+1}', c')'        (M24)

(186 states; G̃ is the paper's transition matrix, written with a tilde to distinguish it from G in (M22); Q has Σ in its top-left block and zeros elsewhere; C_t selects the fixed cells of quarter t, so the free cells are treated as missing). The state at T is known, so the filter starts from a_1 = G̃ S_T, P_1 = Q. The Kalman filter and the state smoother of Durbin and Koopman (2012, Sections 4.3–4.4) give the conditional means and variances:

    v_t = y_t − C_t a_t,  F_t = C_t P_t C_t',  K_t = G̃ P_t C_t' F_t⁻¹,  a_{t+1} = G̃ a_t + K_t v_t,  P_{t+1} = G̃ P_t (G̃ − K_t C_t)' + Q
    E[S_t | C] = a_t + P_t r_{t−1},   Var[S_t | C] = P_t − P_t N_{t−1} P_t        (M25)

with r_t and N_t from the backward recursion. Draws use the mean-correction simulation smoother of Durbin and Koopman (2002): simulate an unconditional path S⁺, smooth both the scenario's values and S⁺'s own values at the fixed cells, and set

    S̃ = Ŝ − Ŝ⁺ + S⁺        (M26)

The filter covariances do not depend on the scenario's values and are computed once per draw.

**The two give the same paths.** With y⁺ = m + A z, the simulation smoother's draw is E[y | C] + A z − A P z = E[y | C] + A (I − P) z, the Gaussian draw (M23) for the same z. Both take z from the same random stream as (M18) (one vector of n standard normals per quarter), so with the same seed they give the same paths up to rounding, and with no fixed cell both reproduce the baseline exactly. The method is a setting (`forecast.conditional_method`, section 9.6); the default is the Gaussian form, which is 5 to 20 times faster here.

**Scenario effect.** For each draw θ (paper, Section 2.3):

    D(θ) = E[y | C, θ] − E[y | θ]

the conditional mean (M22) minus the no-shock path. Both terms are exact, so the spread of D over the draws shows parameter uncertainty only, as for the impulse responses (section 10). The scenario's own fan chart comes from the simulated paths (M23), which include future shocks. The baseline is computed by the same code with no fixed cell, so an empty scenario has an effect of exactly zero.

**Checks** (`reports/conditional_checks.qmd`; tests in `tests/testthat/test-forecast-conditional.R`):
- *Against a textbook computation.* On a known VAR(2) with three variables, both methods equal the conditional normal computed from the joint covariance in companion form, to about 10⁻¹⁵ (means and covariances). On the full model (1,000 draws per COVID treatment, five test scenarios from one cell to 200 cells), the two methods' conditional means agree to 1.2 × 10⁻¹⁰ (levels of about 10³) and their variances to 5 × 10⁻¹⁴ in relative terms; with the same seed their paths differ by at most 1.3 × 10⁻⁹.
- *Simulated against exact.* 5,000 to 20,000 simulated paths against the exact conditional mean and covariance of one draw (593 free cells): joint χ² p-values 0.05–0.32; variance ratios 0.95–1.07. The two lowest p-values (0.054 and 0.065) were re-run with 80,000 paths, with the original seed and two new ones: p = 0.14–0.79, no sign of a bias.
- *Exact properties, per draw, both methods.* Fixed cells are met in every path (to 10⁻¹² on the production draws); fixed cells have zero conditional variance and no conditional variance exceeds the unconditional one; a scenario that fixes cells at their baseline means returns the baseline mean everywhere; scaling a scenario's deviation from the baseline by α scales its effect by α, and effects of deviations on the same cells add up; reordering the variables or the columns of a scenario file changes nothing. Largest relative difference 4.7 × 10⁻¹³ (Gaussian) and 9.3 × 10⁻¹¹ (Kalman).
- *Against the `BVAR` package* (an independent Waggoner–Zha implementation), on the same draws: with two or more fixed cells, its draws match the exact conditional distribution (χ² p 0.10–0.91), and at production size its medians differ from the project's by 1.1 percent of the 90 percent band width (median over cells; at most 4.9 percent), the same as the difference between two seeds of the project's own code (0.9 and 4.6 percent). With a single fixed cell the package is wrong (section 13).
- *Cost.* About 10–70 ms per draw for the Gaussian method and 130–570 ms for the Kalman smoother, which keeps the paper's dense 186 × 186 matrices; about 35–50 s per scenario for 1,000 draws with the default.

### 9.3 Scenario files

A scenario is a CSV file `scenarios/<id>.csv` with an entry in `scenarios/scenarios.yml` (`id`, `label`, `description`, and optionally `short_label`, the name shown on the website, and `soft`, section 9.4). Implemented in `R/scenarios.R`. Every scenario listed in `config/run.yml` is run by the engine and exported to the website (section 15).

| Rule | Detail |
|---|---|
| Layout | One row per quarter, a column `quarter` (for example `2027Q1`), then any subset of the 31 short names of `docs/data.md`, in any order |
| Units | Model units: 100·ln(x) for log series, natural units for rates. Values are levels (absolute paths), not changes, because the VAR is in levels |
| Blank cell | Free (blank or `NA`); the model fills it in |
| Quarters | Calendar quarters after the last observed quarter and within the horizon; rows may skip quarters. A file whose quarters are no longer in the future after a data update is rejected, not shifted |
| Errors (all reported together) | missing `quarter` column; unknown or repeated column; malformed or repeated quarter; quarter out of range; non-numeric or infinite value; all 31 variables fixed in one quarter |

**Helpers.** `to_model_units()` converts natural units (stops on non-positive values of a log series). `path_from_growth()` builds a level path from published growth rates and the last observed values, using exact percent changes: quarterly at an annual rate, g = 100[(x_t/x_{t−1})⁴ − 1], inverted as y_t = y_{t−1} + 100 ln(1 + g/100)/4; four-quarter, g = 100(x_t/x_{t−4} − 1), inverted as y_t = y_{t−4} + 100 ln(1 + g/100). A growth rate whose base quarter is itself free cannot be turned into a level without fixing that base quarter too, so it is an error rather than being filled in from the baseline. `make_scenario_template()` writes the baseline median path of all variables, to edit and trim.

**Soft columns hold values in their measure.** In a column marked soft (section 9.4) every non-blank cell is a target in that column's measure: a four-quarter growth rate in percent for `yoy`, not a level. A level of about 490 in a `yoy` column would be read as a target of 490 percent, and the engine stops because no path comes near it.

**Not supported: conditions on linear combinations.** Each condition fixes one cell. Conditions on combinations of cells, such as a growth rate between two free quarters, or the paper's tax-reform example and its Figure 16a (four-quarter inflation fixed 12 to 20 quarters ahead), are not available; the Gaussian form (M22) could take them as a general matrix R.

**Example files.** `anchor_pce_2pct` is the paper's inflation anchor (section 9.4). `jpm_midyear_2026` and `gs_inflation_oil_2026` translate a few headline numbers from two published 2026 bank outlooks into conditions (`scripts/build_example_scenarios.R` records each number, its source and the mapping); they are rough illustrations of the file format, not those banks' forecasts. `spf_2026q3_2028` tilts the baseline so that the medians of headline PCE inflation (fourth quarter over fourth quarter, `yoy` in 2028:Q4) and real GDP growth (annual average over annual average, `annual_avg` in 2028:Q4) equal the 2028 medians of the Survey of Professional Forecasters' 2026:Q3 survey (2.12 and 2.06 percent); its file is written by hand from the survey's published medians, and its effective sample share is 87.5 percent. `oil_rises` and `oil_normalizes` are built on the SPF-tilted forecast (section 9.4, "Scenarios built on a tilted forecast") and fix only the WTI oil price over 2026:Q3–2027:Q2, on symmetric paths from its 2026:Q2 average ($95.75 a barrel): halfway in logs in 2026:Q4, then the end level in 2027:Q1 and 2027:Q2, either $110 (halfway $102.63) or the 2025 average, $65.46 (halfway $79.17); both are built by `scripts/build_example_scenarios.R`. All six were written for data through 2026:Q2. With data through 2026:Q3 they are used as written, except that the two oil files no longer fix 2026:Q3, now observed (its average was $87.22): they fix 2026:Q4–2027:Q2 only.

### 9.4 Soft conditions: entropic tilting

A soft condition is a target for a feature of the forecast distribution, for example "the median of four-quarter PCE inflation in 2029:Q2 is 2 percent", instead of a fixed value. It is imposed by entropic tilting (paper, Section 4.2; Robertson, Tallman and Whiteman 2005; Krüger, Clark and Ravazzolo 2017): the forecast paths are kept as they are and reweighted, so the model's spread around the target survives instead of being collapsed to a point.

**Problem and solution.** Start from J paths y^(j), each with weight 1/J, and m targets written as moment conditions E[g(y)] = ḡ. Among all weightings that meet them, take the one closest to equal weights in Kullback–Leibler divergence (paper, eqs. 14–16):

    min_w Σ_j w_j log(J w_j)   subject to   Σ_j w_j g(y^(j)) = ḡ,   Σ_j w_j = 1
    w_j = exp(γ' g_j) / Σ_k exp(γ' g_k),   γ = argmin_γ Σ_j exp(γ'(g_j − ḡ)),   g_j = g(y^(j))        (M30)

The objective for γ is convex, so it has a single minimum when the targets can be met (ḡ strictly inside the convex hull of the g_j) and none otherwise.

**Targets.** A target is a statistic of a quantity s(y) computed on each path:
- *median* of s at x: g = 1[s(y) ≤ x], ḡ = 0.5 (the paper's anchor; the default);
- *mean* of s at x: g = s(y), ḡ = x.

The quantity is one variable in one forecast quarter, as a level in model units (`level`) or as an exact percent change of x = exp(y/100): `yoy` = 100[exp((y_t − y_{t−4})/100) − 1], `qoq_ann` = 100[exp(4(y_t − y_{t−1})/100) − 1] or `annual_avg` = 100[Σ_{i=0}^{3} x_{t−i} / Σ_{i=4}^{7} x_{t−i} − 1] (log series only). `annual_avg` is the growth of a four-quarter average over the previous four; in a Q4 quarter it is the growth of the calendar-year average over the previous year's, the form in which annual real GDP growth is usually quoted (for example by the Survey of Professional Forecasters). A base quarter after T is taken from the same path, so a four-quarter rate whose base quarter is also in the forecast is a single soft target; base quarters up to T come from the observed history. In `scenarios.yml` a soft condition is written per variable, `soft: {pce: yoy}`, or in long form `pce: {measure: yoy, statistic: mean}`; unlisted columns stay hard.

**Solver.** Newton's method from γ = 0 on the log of the objective, which has the same minimizer and is better scaled:

    f(γ) = log Σ_j exp(γ'(g_j − ḡ)),   ∇f = Σ_j w_j (g_j − ḡ),   ∇²f = Σ_j w_j (g_j − ḡ)(g_j − ḡ)' − ∇f ∇f'
    γ ← γ − t (∇²f)⁻¹ ∇f        (M31)

with w_j the weights (M30) at the current γ, so the gradient is the current miss of the targets and the Hessian the weighted covariance of the g_j. The step length t is halved until f decreases (backtracking), all sums are computed in log-sum-exp form, and the iterations stop when every miss is below `forecast.tilting.tol`. Before solving, each target on its own must lie strictly inside the range of its g_j (for a median: paths on both sides of the value), otherwise the engine stops with an error; targets that are jointly impossible show up as non-convergence or a singular Hessian, also errors. A target that the paths already meet gives γ = 0 and equal weights.

**Weighted quantiles.** Tilted fan charts use the weighted version of R's default quantile (Hyndman–Fan type 7) of Akinshin (2023, Section 6). With the values sorted, x_(1) ≤ ... ≤ x_(J), and their normalized weights w_(i):

    Q(p) = Σ_i W_i x_(i),   W_i = F(t_i) − F(t_{i−1}),   t_i = Σ_{k≤i} w_(k),   n* = 1 / Σ_i w_i²
    F = CDF of the uniform distribution on [(η − 1)/n*, η/n*],   η = (n* − 1) p + 1        (M32)

With equal weights this is exactly `quantile(x, p)`; values with zero weight drop out.

**Effective sample size.**

    ESS = 1 / Σ_j w_j²        (M33)

the number of equally weighted paths the tilted sample is worth (J with equal weights, 1 if one path carries all the weight). It is reported with every tilted result, and a warning is given below 10 percent of J (`forecast.tilting.ess_warn_share`); a low ESS means the tilted fan rests on few paths. For a single median target the weights take only two values, 0.5/(qJ) on the share q of paths at or below the target and 0.5/((1 − q)J) on the rest, so ESS/J = 4q(1 − q): it depends on where the target lies in the untilted distribution, not on J.

**Hard and soft together.** Hard cells are imposed first, by conditioning (section 9.2); the paths that result are then tilted. Tilting only changes weights, so the hard cells stay exactly met. A tilted baseline is a scenario with soft conditions only.

**Paths and effect.** A scenario with soft targets uses 10 paths per posterior draw (J = 10,000; `forecast.tilting.paths_per_draw`), and so do the baseline and scenarios with hard conditions (`forecast.paths_per_draw`; section 14.4); the paper uses one path per draw. With more paths the weight is spread over more posterior draws (in the anchor example under `precovid_params`, the share of draws carrying half the weight is 17 percent with one path per draw and 33 percent with ten). The effect of a tilted scenario is one path, the weighted mean of its paths minus the baseline mean, without bands, because the weights belong to paths and not to posterior draws. Impulse responses are not tilted.

**Checks** (`tests/testthat/test-tilting.R`, `reports/tilting_checks.qmd`): targets met to the tolerance; weights positive and summing to 1; equal weights when there is no target or the target is already met; on a Gaussian sample, a mean target shifts the mean and leaves the variance unchanged (the closed form); the solution has a smaller divergence than perturbed weights that also meet the targets; hard cells unchanged; impossible targets rejected; growth quantities against hand calculations; weighted quantiles against Akinshin's Example 13 and against `quantile()` with equal weights.

**The paper's inflation anchor** (paper, Figure 16b): the median of four-quarter headline PCE inflation 12 quarters ahead (2029:Q2) set to 2 percent, applied to the baseline of each model (`scenarios/anchor_pce_2pct`).

| | `none` | `precovid_params` |
|---|---|---|
| Untilted median of the target quantity | 1.89% | 0.41% |
| Effective sample size (of 10,000 paths) | 9,970 (99.7%) | 5,947 (59.5%) |
| Newton steps | 2 | 4 |
| Largest miss of the target | 4 × 10⁻¹³ | 4 × 10⁻¹⁶ |

As in the paper, the tilted bands do not pinch at the target quarter: the interquartile range of the inflation quantity keeps growing with the horizon (`precovid_params`: 1.99, 2.67 and 2.70 percentage points at 6, 11 and 12 quarters ahead, against 1.90, 2.26 and 2.29 untilted; under `none` the anchor barely changes the distribution). Hard conditioning at 2 percent for 12–20 quarters ahead (the paper's Figure 16a) is a condition on a linear combination of cells and is not reproduced (section 9.3). A mixed example with oil prices hard and core inflation soft meets all targets, keeps the hard cells to 10⁻¹², and has an ESS of 64 percent (`none`) and 31 percent (`precovid_params`).

#### Scenarios built on a tilted forecast

A scenario is built on the unconditional forecast unless its entry names a `base` (`scenarios/scenarios.yml`): a tilted scenario with soft conditions only. Its paths are then drawn given its own fixed cells, as in section 9.2, and its quantiles are weighted with the base's tilting weights,

  Q(scenario built on b) = weighted quantiles of the paths y_j given the scenario's conditions, with weights w_j of the tilt of b, j = 1, ..., J.

This is valid because every scenario of a run uses the same seed and the same number of paths per draw. Path j of the scenario and path j of the base, which fixes no cell and so is an unconditional path, come from the same standard normal draws z (M23); the scenario's are adjusted to meet its conditions. The tilt reweights these shock histories, and the scenario runs the reweighted histories through its own conditions. Its median minus the base's median is the effect of the scenario's own conditions under the base's distribution. The base's targets are not imposed again: they move with the conditions (a higher oil path, for example, raises inflation in the target quarter above the survey value that the base meets). Rules: the base has soft conditions only and no base of its own, and the scenario has no soft conditions of its own; its effective sample size is the base's. This construction is not in the paper. Tests check it: with no condition of its own, a scenario built on a base reproduces the base exactly; with conditions, its paths are the hard scenario's paths, weighted with the base's weights.

### 9.5 Simulation recovery

Every part of the estimator and the engine is checked against an exact answer above. This check runs them together where the truth is known (`scripts/check_recovery.R`, `reports/recovery_checks.qmd`). A 6-variable VAR(5) (real GDP, core PCE prices, the unemployment rate, the 2-year yield, the S&P 500 and realized volatility) is set to the posterior mean of the model of section 3 estimated on those six series of the real data; its largest root is 0.99991, so it is highly persistent but stable. One hundred histories of 214 quarters are simulated from it, the model is estimated on each with the prior and sampler of sections 4–7 (1,000 draws), and the bands are checked against the truth.

| Quantity | Coverage of 68% bands | Coverage of 90% bands |
|---|---|---|
| Forecasts 1–8 quarters ahead, pooled | 66.8% | 88.7% |
| Recursive 2-year yield shock, responses at h = 1–12 | 73.9% | 93.8% |
| Coefficients on own first lags | 29.5% | 59.5% |
| Coefficients on other variables' lags | 75% | 93% |
| Σ: variances, covariances | 71%, 70% | 91%, 90% |

One coverage rate over 100 replications has a standard error of about 4.7 percentage points at 68 percent and 3.0 at 90 percent. The forecast bands have close to nominal coverage, so the estimation and the forecasts are consistent with each other. Coverage is not expected to be exactly nominal: these are credible bands under a shrinkage prior, not frequentist confidence intervals. Two departures are kept in view:
- Own first-lag coefficients are biased down (by 0.045 on average, against a truth of about 0.99) and the longer own lags up, so total own persistence is nearly unbiased (about −0.01) but split differently between the lags. The impact-plus-one-quarter responses (h = 1), which depend on the first lags only, are undercovered (2-year yield 44% and 78%; real GDP 64% and 85%); from h = 2 coverage is at or above nominal.
- The largest root of the estimated posterior mean has a median of 1.0000 across replications (range 0.9996–1.0036) although the truth is stable, so estimated roots at or slightly above 1 (section 16) are what a stable, very persistent economy would produce.

### 9.6 Engine settings and outputs

The engine is stage 4 of the pipeline (`R/engine_run.R`, `scripts/run_engine.R`). For each estimated model it reads the saved draws and the forecast history (the latest data, section 11), and writes the baseline, each listed scenario and the impulse responses (section 10) to `outputs/<run_id>/engine/`: the simulated paths, each draw's mean path, the scenario effects, the response draws, quantile tables (5, 16, 25, 50, 75, 84 and 95 percent) and, for tilted scenarios, the weights and the effective sample size. Its `meta.json` records the settings, the seed and the MD5 checksums of the draws and of each scenario file.

**Generalized responses of the other variables.** For the published model only, a further step (stage 5, `run_engine_irf_extra()`) computes the generalized response (M28) to a shock in every variable not listed in `irf.girf_shocks`. Each is computed draw by draw and summarized straight away into the export's quantile tables, with the four-quarter changes and the spreads taken per draw first (section 15), and its draws are not kept: about 10 MB each, against a table of about 0.1 MB. They are written to `engine/irf_extra.rds` with the checksum of the draws they come from, which the export checks.

**Common random numbers.** The engine's seed is the master seed plus 3,000, the same for the baseline and every scenario, so their simulated paths are built from the same random numbers and the differences between their fan charts are not sampling noise.

**Reference outputs.** A small known VAR with a fixed seed is run through the engine (both methods, and a tilted case) and its results are stored with the tests (`tests/testthat/fixtures/engine_reference.rds`); a test fails if a code change moves any of them by more than 10⁻⁸.

| Setting (`config/run.yml`) | Value | Meaning |
|---|---|---|
| `forecast.horizon` | 20 | forecast horizon H, quarters |
| `forecast.conditional_method` | gaussian | scenario algorithm: gaussian (M21)–(M23) or kalman (M24)–(M26) |
| `forecast.seed_offset` | 3000 | engine seed = `mcmc.seed` + 3000 |
| `forecast.paths_per_draw` | 10 | simulated paths per posterior draw, baseline and hard scenarios |
| `forecast.tilting.paths_per_draw` | 10 | paths per posterior draw for scenarios with soft targets |
| `forecast.tilting.tol` | 1e-10 | largest acceptable miss of a target (M31) |
| `forecast.tilting.maxit` | 100 | Newton steps (M31) |
| `forecast.tilting.ess_warn_share` | 0.1 | warning when the ESS (M33) is below this share of J |

---

## 10. Impulse responses

The paper defines an impulse response as a scenario minus the baseline (Section 2.3): for each posterior draw, the forecast mean given a condition on the first forecast quarter minus the forecast mean without it, with no simulated shocks. For the two conditions used here the difference has a closed form, a Cholesky impulse response, which is what the code computes (`R/irf.R`).

A shock is summarised by its impact vector a, the forecast errors it causes in T + 1. The response h quarters after impact (h = 0 is T + 1) is

    r_h = Ψ_h a,   h = 0, 1, ..., 40        (M27)

For a Cholesky identification with ordering π and a shock to variable j, a is the column of the lower Cholesky factor of Σ (variables in the order π) that belongs to j, rescaled so that a_j = 1; variables ordered before j do not move on impact. Because the VAR is linear, r_h does not depend on the forecast origin.

**Generalized impulse response of variable j** (Koop, Pesaran and Potter 1996): condition y_{j,T+1} on its baseline plus 1 and leave everything else free. Every variable moves on impact by its regression coefficient on j's forecast error:

    a = Σ e_j / Σ_jj        (M28)

which equals the Cholesky response with j ordered first. It has no structural interpretation: it is the mix of shocks that usually moves j (paper, Section 3.1). Computed, with every draw kept, for real GDP and the 2-year yield (the paper's Figures 6 and 8), oil, the S&P 500 and the Baa yield; for the published model it is also computed for every other variable, kept as quantile tables only (section 9.6), so the website can offer a generalized response to a shock in any variable.

**Recursive shock to the 2-year yield**: the slow (economic) variables held at their baseline in T + 1, the 2-year yield at baseline plus 1, the fast (financial) variables free (paper, Section 2.3 and Figure 9). This equals the Cholesky response with ordering (slow, 2-year yield, fast) and the shock on the yield:

    a = Cov(ε_{T+1}, ε̃) / Var(ε̃),   ε̃ = ε_{gs2,T+1} − E[ε_{gs2,T+1} | ε_{slow,T+1}]        (M29)

where ε̃ is the part of the yield's forecast error not explained by the slow variables' errors. Hence a_slow = 0 and a_gs2 = 1; a is the column of the Cholesky factor for the ordering (slow, gs2, fast) that belongs to gs2, rescaled to a unit shock. The response does not depend on the order within each block. The paper's text does not list the blocks; its Table 1 orders the economic variables before the financial ones (appendix, p. 52), and that order is used:
- slow (21): rgdp, rpce, rresinv, rnresinv, rexp, rimp, rgov, rfedgov, gdpdefl, pce, corepce, cpi, corecpi, comphr, payems, urate, indpro, caputil, houst, rdpi, umich;
- policy: gs2;
- fast (9): gs5, gs10, aaa, baa, usd, sp500, wti, gsci_nonenergy, rvol.

**Shock size and bands.** Unit shocks: 1 log point (about 1 percent) for log series, 1 percentage point for rates; the paper's 25 basis point yield exercises are 0.25 times the unit responses. One response per posterior draw, so the bands (median, 68 and 90 percent) show parameter uncertainty only.

**Reading.** An impulse response is a statement about the model's forecasts, conditional on the identification: the generalized responses carry no economic label, and the recursive response depends on the stated split into slow and fast variables. Every draw has a largest root slightly above 1 (section 16), so responses of log levels keep widening at long horizons; across draws, the largest root and the size of the 40-quarter response are almost uncorrelated (rank correlation −0.01 to 0.18), so the explosive roots are not what drives the long-horizon spread.

**Checks** (`tests/testthat/test-irf.R`; the figures below are for the production draws, from a check not included in this repository): the closed forms (M27)–(M29) equal the scenario engine's exact conditional means (M22) with the corresponding conditions, to a relative 1.1 × 10⁻¹³, with the slow variables exactly 0 on impact; reversing or shuffling the order within blocks changes the responses by at most 1.8 × 10⁻¹⁵; the `BVAR` package's `irf()` (Cholesky in data order, one-standard-deviation shocks) agrees to 2.3 × 10⁻¹⁴ after reordering and rescaling. All six identifications for 1,000 draws take about 17 seconds.

| Setting (`config/run.yml`) | Value | Meaning |
|---|---|---|
| `irf.horizon` | 40 | responses at h = 0, ..., 40 |
| `irf.girf_shocks` | rgdp, gs2, wti, sp500, baa | generalized impulse responses computed (M28) |
| `irf.recursive.slow` | rgdp, rpce, rresinv, rnresinv, rexp, rimp, rgov, rfedgov, gdpdefl, pce, corepce, cpi, corecpi, comphr, payems, urate, indpro, caputil, houst, rdpi, umich | slow block (M29) |
| `irf.recursive.policy` | gs2 | policy variable (M29) |
| `irf.recursive.fast` | gs5, gs10, aaa, baa, usd, sp500, wti, gsci_nonenergy, rvol | fast block (M29) |

---

## 11. COVID treatment

The paper's sample ends before the pandemic. The quarters 2020:Q2–2020:Q3 contain changes far outside the range of the rest of the sample (real GDP −8.2 and +7.5 log points), which affect the estimated Σ and the hyperparameters. Three models are estimated side by side (`R/covid.R`):

| Method | Estimation sample | Forecasts start from |
|---|---|---|
| `none` | 1973:Q1 to the last quarter (2026:Q2) | the latest data |
| `precovid_params` | 1973:Q1 to 2019:Q4 | the latest data |
| `exclude_window` | 1973:Q1 to the last quarter, without the equations of 2020:Q1–2021:Q4 | the latest data |

Under `precovid_params` the coefficients reflect pre-2020 dynamics and are applied to the post-2020 history. Under `exclude_window` the equations (y_t on its p lags) of the quarters in the window are left out of the estimation and every other equation is used, including those after the window; the left-out quarters' data still enter later equations as lags (with p = 5, the first equation used after the window, 2022:Q1, has 2020:Q4–2021:Q4 as lags, so 2020:Q1–Q3 enter no equation). The marginal likelihood, the posterior and the ψ starting values are computed on the equations used. Excluding the pandemic observations from estimation follows Schorfheide and Song (2024); Lenza and Primiceri (2022) note that it suits the estimation of the parameters but understates uncertainty for forecasts made during the pandemic. Every dropped quarter or left-out equation is listed in the run's `covid.rds`. Other treatments (interpolating the pandemic quarters, or scaling their volatility as in Lenza and Primiceri 2022) are not implemented.

**Treatment used for the published forecasts: `exclude_window`.** It keeps the responses of the pre-pandemic estimates, which are close to the paper's, and uses the data from 2022 on. Under `none` the 2020 quarters halve the activity responses and leave the pre-2020 payroll and unemployment residuals autocorrelated (lag-1 autocorrelation 0.40 and 0.30); under `precovid_params` the forecasts made from 2022–2024 with the same estimates were too low by about 1.5 percentage points of GDP growth, twice the error of the other two (section 14). The other two models remain in the pipeline for comparison.

---

## 12. Differences from the paper

| Item | Paper | This model | Source or reason |
|---|---|---|---|
| Data | 31 series, 2012-based vintages | the same 31 concepts, current vintages, with documented substitutions | `docs/data.md`, section 5 |
| Sample end | before the pandemic | 2026:Q2, with three COVID treatments | section 11 |
| Hyperpriors of λ, µ, ψ | not stated | GLP's | section 4.4 |
| Units of the ψ hyperprior | not stated | GLP's numbers, calibrated by GLP for 4·ln data, applied to 100·ln data | section 4.4; the mode is insensitive to the choice |
| Hyperparameter bounds | not stated | λ [0.0001, 5], µ [0.0001, 50], ψ start ×/÷ 100 | `BVAR` package defaults; never binding (section 8.3) |
| Intercept prior | not stated | diffuse, v_c = 10⁷ | GLP |
| ȳ_0 in the sum-of-coefficients dummies | not stated | mean of the first p observations | GLP |
| Sampler details (log scale, starts, tuning, chain lengths) | not stated | section 7 and 8 | GLP, Appendix B; diagnostics |
| Software | not stated | estimator and engine written for this project; the `BVAR` package is used only for checks | section 13 |
| Scenario algorithm | Kalman simulation smoother (Section 2.3) | Gaussian conditioning on the stacked shocks by default; the paper's Kalman simulation smoother as an option | same conditional distribution and, for the same seed, the same paths (section 9.2); faster |
| Conditions on linear combinations of cells | used (tax-reform example; Figure 16a) | not supported: each condition fixes one cell | section 9.3 |
| Blocks of the recursive yield shock | economic before financial variables; not listed in the text | the order of Table 1 (section 10) | paper, Table 1 and appendix |
| Shock size | 25 basis points for the yield exercises | unit shocks | responses are linear in the shock size (section 10) |
| Entropic tilting | one path per posterior draw; unconditional forecasts | 10 paths per draw; also combined with hard conditions | section 9.4 |
| Quantiles of tilted forecasts | not stated | weighted Hyndman–Fan type 7 (Akinshin 2023) | section 9.4 |

---

## 13. Why the estimator and engine are written for this project

The `BVAR` R package (version 1.0.5; Kuschnig and Vashold 2021) implements the same model class and the GLP hierarchical prior. Checks against formulas written from GLP found five problems that matter for this model. They were found with the scripts named below; the reports they produce are regenerated by running the scripts.

1. **Marginal likelihood of the dummy observations.** `BVAR` evaluates p(Y⁺ | h) in (M15) at the prior mean b instead of at the dummies' posterior mean. This is exact only when the prior mean fits the dummy observations, which fails for a variable with a white-noise prior mean combined with the sum-of-coefficients prior: the paper's specification for realized volatility. On a 6-variable version of the model, the package's log marginal likelihood differs from the formula by 2.3, 26.5, 58.6 and 133 log points at µ = 1, 0.1, 0.01 and 0.0001, and its mode search drives µ to its lower bound. With every variable given a random-walk prior mean, the package and the formula agree to about 10⁻⁸ (`scripts/check_priors.R`).
2. **Proposal of the Metropolis step.** `BVAR`'s proposal is diagonal and scaled by a heuristic rather than by the Hessian (package GitHub issue #94). For hyperparameters with large values the heuristic gives proposal standard deviations of order 10⁻⁶⁰ and smaller, so those chains do not move (a check not included in this repository; see the README). With 31 ψ's, several in the hundreds, the full model would not sample them. The mode search also stops at 100 iterations, with no option to raise it.
3. **Unconditional forecast bands.** `predict()` adds to each quarter's no-shock forecast a noise term with covariance Σ² instead of shocks with covariance Σ propagated through the dynamics (package GitHub issue #86). Mean forecasts are correct; bands are not. With this model's draws, the package's one-quarter-ahead bands for real GDP are 5.6 times too wide (`scripts/check_forecast.R`).
4. **Selecting conditioned variables by name.** The package matches variable names partially (`grep()`), so "pce" also matches "rpce" and "corepce". Conditioned variables must be passed by position.
5. **A single fixed cell in a conditional forecast.** The package inverts the singular values of the restriction matrix with `diag(1/d)`; with one fixed cell d is a single number, and `diag()` of a number x returns an identity matrix of size floor(x) instead of the 1 × 1 matrix x. When 1 ≤ 1/d < 2 the fixed cell moves by d times the requested amount, without a warning (real GDP +1 in T + 1 moved by 0.81 under `none` and 0.54 under `precovid_params`); otherwise the call fails. A generalized impulse response computed as a conditional forecast has exactly this shape (`tests/testthat/test-forecast-conditional.R`).

The project therefore computes the hyperparameter posterior, the sampler, the coefficient draws, the forecasts, the scenarios and the impulse responses in its own code. `BVAR` remains a cross-check where it is exact: the log posterior with every variable a random walk, mean forecasts, conditional forecasts with two or more fixed cells, and impulse responses; there the implementations agree.

---

## 14. Validation

The checks of the code and the engine (sections 7–10) compare each computation with an exact result. The checks summarized here ask whether the estimated model behaves sensibly. Each has a script and a report (`scripts/check_*.R`, `reports/*.qmd`); the numbers are for the 2026m09 vintage.

### 14.1 Estimation

- **Residuals** (`reports/estimation_diagnostics.qmd`). One-quarter-ahead residuals at the posterior mean of B and Σ. Beyond ±5 standard deviations: 15 cells under `none` (ten in 2020:Q2) and 7 under `precovid_params` (2008:Q4 prices and commodities, the 1987 equity crash); beyond ±3, about four times the normal-distribution count, the fat tails usual in macroeconomic shocks. Residual autocorrelation (Ljung–Box to lag 8) is significant at 5 percent for 5 of 31 variables under `precovid_params`; under `none`, the pre-2020 payroll and unemployment residuals have lag-1 autocorrelations of 0.40 and 0.30. Residual correlations are as expected (CPI and core CPI 0.56, PCE and core PCE 0.70, the 5- and 10-year yields 0.97).
- **Hyperparameters across origins** (`reports/hyper_stability.qmd`). The mode was re-estimated on samples ending in 2024:Q2, 2016:Q2 and 2007:Q4 (method `none`). λ hardly moves; µ moves by +0.05, +1.8 and +4.5 posterior standard deviations. Fixing the hyperparameters at the current values changes forecasts from 2024:Q2 by less than simulation noise, and from 2007:Q4 by about a quarter of a percentage point of core PCE inflation. Exercises at past origins therefore re-estimate them.
- **Sensitivity to the hyperparameters** (`reports/prior_sensitivity.qmd`). With λ or µ at the edge of their posterior (±2 standard deviations) forecasts change by up to 0.2 percentage points (year-on-year growth or level); with µ halved, by up to 0.6 (unemployment). λ and ψ matter little; impulse responses in growth terms change by less than 0.1 percentage point.

### 14.2 Impulse responses

`reports/irf_review.qmd`. For 32 expected responses (signs within 12 quarters, for the recursive 2-year yield shock, oil +10%, S&P 500 −10%, real GDP +1% and a Baa spread widening), the 68 percent band excludes zero with the expected sign for 30 under `exclude_window` and `precovid_params` and 29 under `none`; none has the opposite sign. The Baa spread widening (Baa +1 percentage point with the 10-year yield held at its baseline on impact) is a conditional forecast, since a generalized Baa shock mostly raises all yields. Against the paper's figures (medians read from the charts), under `exclude_window`: after a 1 percent GDP shock (Figure 6) GDP peaks at 1.19 (paper about 1.1), payrolls at 0.73 (0.65), unemployment falls by 0.36 (0.35), core PCE prices are 0.27 higher after three years (0.26); after the recursive 2-year shock of 25 basis points (Figure 9) GDP is unchanged one quarter after impact and falls by 0.07 (0.09), Baa rises by 0.14 (0.14). GDP responses are more persistent than the paper's (0.92 three years after the GDP shock, paper 0.5).

### 14.3 Forecasts

`reports/forecast_review.qmd` and `reports/covid_comparison.qmd`; origin 2026:Q2.

- **Against the Survey of Professional Forecasters** (2026:Q3 survey, whose panel had data through 2026:Q2). The SPF medians lie between the 30th and 80th percentile of the `exclude_window` forecast distribution for the comparable variables, except CPI inflation in 2026:Q3: the model carries part of the 6.1 percent annualized CPI inflation of 2026:Q2 into the next quarter, the SPF (which also saw July's data) does not. The model expects weaker activity than the SPF (GDP growth about 1 percentage point lower, annualized, over 2026:Q4–2027:Q3) and unemployment rising to about 4.7 percent by 2027:Q3, against a flat 4.3.
- **Forecasts from 2022–2024** made with the same estimates (no-shock paths, against the outcomes): GDP growth too low by 0.6–0.9 percentage points (year-on-year) under `exclude_window` and `none`, which were estimated on data that include these years, and by 1.5–1.7 under `precovid_params`.
- **Long horizons.** Inflation and interest rates drift down at long horizons (section 16): under `exclude_window`, the median of core PCE inflation is 1.5 percent three years ahead and 0.6 percent at ten years; 12 percent of the paths have a negative 2-year yield three years ahead. Growth and unemployment settle at sensible values (GDP growth about 2 percent, unemployment about 5 percent at ten years).
- **Inflation anchor.** Tilting the median of PCE inflation to 2 percent three years ahead keeps 91 percent of the paths effective under `exclude_window`. For a median target the effective sample share is exactly 4q(1 − q), where q is the untilted share of paths below the target.

### 14.4 Simulation noise

With one simulated path per posterior draw, changing only the seed moved the published medians by up to 0.3 percentage points of GDP growth (year-on-year) and 0.5 (quarterly, annualized). The baseline and hard scenarios therefore use ten paths per draw (`forecast.paths_per_draw`), which reduces this noise by a factor of about three (√10).

---

## 15. Export and the website

The website comes in two versions that show the same views and the same numbers. Neither runs model code: both read the bundle, written by stage 5 of the pipeline (`R/export.R`, `scripts/export.R`) for one estimated model, the published one (`export.run`). The bundle holds quantile tables (5, 16, 25, 50, 75, 84 and 95 percent) of the baseline and of every scenario in `config/run.yml`, the observed history, the impulse responses and metadata (the run identifier, data vintage and settings). It is about 4 MB.

- **The static site** (`site/`), published with GitHub Pages at <https://zayat-alex.github.io/large-bvar-usa/>. `site/build_site.R` writes the bundle's tables to two data files, `site/data/core.json` (about 0.8 MB) and `site/data/irf.json` (about 3.9 MB, loaded when the impulse-response tab first opens), and copies the third-party libraries into `site/vendor/`: Bootstrap with the app's theme, plotly.js and selectize. The page (`site/index.html`, `site/site.js`) builds the charts, tables and downloads in the browser. The workflow `.github/workflows/pages.yml` publishes it when `site/` changes on `main`.
- **The Shiny app** (`app/`), run locally in R with `shiny::runApp("app")` after `scripts/publish_bundle.R` has copied the bundle into `app/bundle/`.

### 15.1 Display transforms

A log series is stored as y = 100 ln x (`docs/data.md`, section 2). The site shows it in one of four transforms, with x = exp(y/100):

| Transform | Formula | Used for |
|---|---|---|
| level | y | model units (not shown on the charts) |
| natural | x | levels in published units: dollars, index points |
| quarterly, annualized | 100[(x_t/x_{t−1})⁴ − 1] | the reader's growth setting |
| year over year | 100(x_t/x_{t−4} − 1) | the reader's growth setting (default) |

These are the published conventions (exact percent changes, as in `path_from_growth()`, section 9.3, and the tilting measures, section 9.4), not the log approximations 4Δy and Δ₄y. In the first forecast quarters the base quarter t − 1 or t − 4 is observed and comes from the history. Rates, ratios and indices kept as published (`raw` in `docs/data.md`) are shown in levels. Some log series are prices or quantities usually read as levels; they are shown in natural units whatever the growth setting (`export.shown_in_levels`): the oil price in dollars a barrel, the S&P 500, the dollar index, the non-energy commodity index and housing starts. Three spreads are formed path by path, in percentage points: the 10-year minus the 2-year Treasury yield, Baa minus Aaa, and Baa minus the 10-year yield.

**Transform first, summarize last.** A transform is applied to every simulated path (or impulse-response draw) and the quantiles are taken afterwards. The quantiles of a growth rate are not the growth rates of the quantiles of the level: the path at the 5th percentile of the level is not the path at the 5th percentile of growth. Quantiles of a tilted scenario are weighted by its tilting weights (M32).

| Setting (`config/run.yml`) | Value | Meaning |
|---|---|---|
| `export.run` | exclude_window | COVID treatment of the published model (section 11) |
| `export.horizon` | 12 | forecast quarters exported and shown |
| `export.shown_in_levels` | wti, sp500, usd, gsci_nonenergy, houst | log series always shown in natural units |

### 15.2 What the site shows

- **Horizon.** Twelve quarters. Longer horizons mostly extrapolate trends in the sample (section 16) and are computed for diagnostics only.
- **Forecasts.** For each listed scenario, the median path and, if selected, the 68 and 90 percent bands; the observed history joins the forecasts at the last observed quarter. The site opens with the forecast tilted to the Survey of Professional Forecasters, labelled "Baseline", and the two oil scenarios; the unconditional forecast and the other scenarios can be listed on the Settings tab. The reader chooses which scenarios are listed (all those in the bundle can be), renames them and changes their colours. The static site remembers these and the other display choices in the reader's browser (local storage, nothing is sent) and restores them on the next visit; the Shiny app keeps them for the visit.
- **Differences.** In the difference view each line is a scenario's median minus the median of a baseline the reader chooses (any listed scenario, the SPF-tilted "Baseline" by default), quarter by quarter, in the chart's transform. It is a difference between two conditional forecasts, not a causal effect, and a difference of medians rather than the paper's per-draw difference of means (section 9.2), so it has no band; each scenario's own uncertainty is on its fan chart. All scenarios share the engine's random numbers (section 9.6), so the simulation noise in these differences is small: a standard error of about 0.01–0.03 percentage points of four-quarter growth, measured by splitting the draws into halves, against differences of 0.06–0.3 in the example scenarios.
- **What a scenario assumes.** Cells a scenario fixes are marked on level charts (hollow circles); the median targets of a tilted scenario are marked where the chart shows their measure (hollow diamonds). A fixed level is not marked on a growth chart, since its growth also depends on the quarter before.
- **Impulse responses.** The recursive 2-year yield shock and the generalized response to a shock in each of the 31 variables (section 10), in levels or, for log series, as the change over four quarters of the response (the response of the four-quarter log change, in percentage points), up to 20 quarters, scaled by a shock size the reader sets. Responses are linear, so scaling is exact; a negative size swaps the lower and upper bands.
- **Tables and downloads.** Every chart has a table of its numbers and each tab a table of all its charts, with a CSV download, in the current settings. The observed history of series whose publishers restrict redistribution is shown on the charts with its source credited but left out of tables and downloads (`docs/data.md`, section 6); the model's forecasts of those series are included.

### 15.3 Checks

The tests (`tests/testthat/test-export.R`, `test-app.R`) check the transforms against a hand calculation and against the inverse of `path_from_growth()`, the bundle's level quantiles against the engine's, the chart and table numbers against the bundle, and every chart in every view. Website checks (not included in this repository; see the README) compare what the site shows, for the published model, with quantiles recomputed from the engine's simulated paths and impulse-response draws, and with the data tables: every comparison agrees to the last digit (the history to 5 × 10⁻¹², the rounding of the data file). They also confirm that no table offers restricted history, that the app starts from a copy of its folder alone, the colour contrast and the labels of the controls. The static site's data files are tested against the app's chart data (`tests/testthat/test-site-data.R`), and further checks, also not included, compare the two versions directly: the 335 forecast and impulse-response tables the static site builds, and its CSV downloads, are identical to the app's.

---

## 16. Known limitations

- **Explosive roots.** In every posterior draw the largest eigenvalue modulus of the companion matrix (M3) is slightly above 1 (median about 1.012). At the posterior mean given the mode, 12 eigenvalues lie within 10⁻⁴ of 1: a cluster of unit roots from the random-walk prior mean and the sum-of-coefficients prior, which parameter uncertainty scatters around 1. No restriction is imposed, as in the paper. The simulation recovery (section 9.5) finds estimated roots around 1 even when the true model is stable.
- **Long-run drift.** A VAR in levels with intercepts and roots near 1 extrapolates trends in the sample. Forecasts of Treasury yields keep falling at long horizons (still about 0.07 percentage points a year 15–20 years ahead, close to the 1973–2019 average fall of 0.10), and the pre-2020 model's core inflation forecast drifts towards zero over several years. The paper calls its sum-of-coefficients prior "rudimentary" (footnote 4) and anchors long-run inflation by entropic tilting (Section 4.2; section 9.4 here). Under `exclude_window`, core PCE inflation settles at about 0.6 percent and the 2-year yield at about 0.2 percent ten years ahead (section 14.3); the drift comes from the 1973–2019 trends, which every COVID treatment keeps. For this reason the website shows forecasts only up to 12 quarters ahead (section 15); longer horizons are computed for diagnostics only.
- **The COVID treatment matters.** The three treatments give different baselines: the median of four-quarter PCE inflation in 2029:Q2 is 1.89 percent under `none`, 1.40 under `exclude_window` and 0.41 under `precovid_params`, and scenario effects differ accordingly. The published forecasts use `exclude_window` (section 11).
- **Pandemic outliers** under `none` enter as very large innovations, which raise Σ and affect the hyperparameters (section 11).
- **Scenarios are conditional forecasts.** They show what the model expects given the conditions, not what causes what (section 9.2). Conditions on linear combinations of cells are not supported (section 9.3).
- **Impulse responses depend on the identification.** Generalized responses have no structural interpretation, and the recursive yield response depends on the split into slow and fast variables (section 10).
- **First-lag coefficients.** In the simulation recovery the model splits persistence between the first and later own lags differently from the truth, and responses one quarter after impact are undercovered (section 9.5).
- **Tilting with many or distant targets.** When targets lie far from the untilted distribution, or many are imposed at once, a few paths carry most of the weight and the tilted results become noisy; the effective sample size (M33) measures this. Median targets make the weights jump at the target value, and the tilted joint distribution can have several modes (as in the paper's Section 4.2).
- **Website speed.** On the development machine, with no network time, the overview's six charts took about 4–7 seconds to appear on the static site, which draws them in the browser, and about 8–18 seconds in the Shiny app, which builds each chart on the server. The numbers are not affected.

---

## 17. References

- Akinshin, A. (2023). "Weighted Quantile Estimators." arXiv:2304.07265.
- Bańbura, M., D. Giannone and M. Lenza (2015). "Conditional Forecasts and Scenario Analysis with Vector Autoregressions for Large Cross-Sections." *International Journal of Forecasting* 31(3), 739–756.
- Brooks, S. P. and A. Gelman (1998). "General Methods for Monitoring Convergence of Iterative Simulations." *Journal of Computational and Graphical Statistics* 7(4), 434–455.
- Crump, R. K., S. Eusepi, D. Giannone, E. Qian and A. Sbordone (2021). "A Large Bayesian VAR of the United States Economy." Federal Reserve Bank of New York Staff Report 976.
- Durbin, J. and S. J. Koopman (2002). "A Simple and Efficient Simulation Smoother for State Space Time Series Analysis." *Biometrika* 89(3), 603–616.
- Durbin, J. and S. J. Koopman (2012). *Time Series Analysis by State Space Methods*, 2nd edition. Oxford University Press.
- Gelman, A. and D. B. Rubin (1992). "Inference from Iterative Simulation Using Multiple Sequences." *Statistical Science* 7(4), 457–472.
- Geweke, J. (1992). "Evaluating the Accuracy of Sampling-Based Approaches to the Calculation of Posterior Moments." In *Bayesian Statistics 4*, Oxford University Press, 169–193.
- Giannone, D., M. Lenza and G. E. Primiceri (2015). "Prior Selection for Vector Autoregressions." *Review of Economics and Statistics* 97(2), 436–451.
- Giannone, D., M. Lenza and G. E. Primiceri (2019). "Priors for the Long Run." *Journal of the American Statistical Association* 114(526), 565–580.
- Koop, G., M. H. Pesaran and S. M. Potter (1996). "Impulse Response Analysis in Nonlinear Multivariate Models." *Journal of Econometrics* 74(1), 119–147.
- Krüger, F., T. E. Clark and F. Ravazzolo (2017). "Using Entropic Tilting to Combine BVAR Forecasts with External Nowcasts." *Journal of Business & Economic Statistics* 35(3), 470–485.
- Kuschnig, N. and L. Vashold (2021). "BVAR: Bayesian Vector Autoregressions with Hierarchical Prior Selection in R." *Journal of Statistical Software* 100(14).
- Lenza, M. and G. E. Primiceri (2022). "How to Estimate a Vector Autoregression after March 2020." *Journal of Applied Econometrics* 37(4), 688–699.
- Robertson, J. C., E. W. Tallman and C. H. Whiteman (2005). "Forecasting Using Relative Entropy." *Journal of Money, Credit and Banking* 37(3), 383–401.
- Schorfheide, F. and D. Song (2024). "Real-Time Forecasting with a (Standard) Mixed-Frequency VAR During a Pandemic." *International Journal of Central Banking* 20(4), 275–320.
- Waggoner, D. F. and T. Zha (1999). "Conditional Forecasts in Dynamic Multivariate Models." *Review of Economics and Statistics* 81(4), 639–651.
- `BVAR` source and issues: https://github.com/nk027/bvar
