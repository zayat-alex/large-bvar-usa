# US Economy – BVAR Forecasts

A large Bayesian vector autoregression (BVAR) of the US economy, with a scenario engine and an interactive website. It follows Crump, Eusepi, Giannone, Qian and Sbordone (2021), "A Large Bayesian VAR of the United States Economy" (Federal Reserve Bank of New York Staff Report 976): 31 quarterly series in levels, five lags, a Minnesota prior with a sum-of-coefficients prior whose tightness is estimated from the data (Giannone, Lenza and Primiceri 2015).

The project estimates the model, forecasts from the latest data, and builds three kinds of results:

- **Forecasts** with 68% and 90% Bayesian credible bands.
- **Scenarios**: conditional forecasts given fixed paths for some variables (for example an oil price path), or forecasts reweighted so that their medians meet targets, such as a survey's inflation forecast (entropic tilting, as in the paper's Section 4.2).
- **Impulse responses**: generalized responses to a shock in any variable, and a recursive shock to the 2-year Treasury yield.

Everything is written in R; the estimator and the scenario engine are written for this project and checked against closed forms, simulations and the `BVAR` package (`docs/methodology.md`, sections 13 and 14).

## The website

**<https://zayat-alex.github.io/large-bvar-usa/>** shows the forecasts, scenarios and impulse responses of the published model. Readers choose which scenarios to show and which to treat as the baseline, switch between levels and growth rates, view differences between scenarios, and download the numbers behind every chart.

The forecasts use data through 2026Q3 (FRED-QD vintage 2026m10, other sources as of 9 November 2026). The model was developed and checked on data through 2026Q2; the update to 2026Q3 re-estimated the published model and rebuilt its forecasts without rerunning the checks (`docs/methodology.md`, introduction).

The website is a static page (`site/`), published with GitHub Pages: it reads tables exported from the model run and draws the charts in the browser, and it runs no model code. The same views are available as a Shiny app (`app/`), which reads the same export and runs locally in R (`docs/methodology.md`, section 15).

How to read it:

- Scenarios are **conditional forecasts**: what the model expects given a scenario's assumptions. They are not predictions, and a difference between two scenarios is not a causal effect.
- The bands are **credible intervals under the model**: they cover parameter and shock uncertainty, not the uncertainty about the model itself.
- Forecasts are shown **12 quarters ahead**. A model in levels mostly extrapolates the trends of its sample beyond two to three years.
- The estimation leaves out the equations of 2020Q1–2021Q4 (the pandemic's extreme movements; `docs/methodology.md`, section 11).
- Non-energy commodity prices are a **proxy**: the World Bank's non-energy index stands in for the S&P GSCI non-energy index used in the paper.
- Impulse responses depend on the identification: generalized responses describe co-movement, not structural causes.

## Repository layout

| Path | Contents |
|---|---|
| `R/` | functions: data, COVID treatment, estimation, forecasting engine, impulse responses, export, checks |
| `_targets.R` | the pipeline (`targets`): data, COVID treatment, estimation, engine, export |
| `config/` | `run.yml` (every setting, with the seeds) and `series.csv` (the 31 series) |
| `scenarios/` | scenario files and their index |
| `scripts/` | entry points for each stage and for the check reports |
| `site/` | the website: the static page, its data (written from the export by `site/build_site.R`) and third-party libraries |
| `app/` | the same website as a Shiny app, run locally |
| `.github/workflows/` | publication of the website with GitHub Pages |
| `reports/` | Quarto sources of the validation and check reports; their results in `reports/results/` |
| `tests/` | automated tests (`testthat`) |
| `docs/` | `methodology.md` (the model, the engine, the website) and `data.md` (the series, sources and construction) |
| `data/raw/SOURCE.md` | the raw files used, their sources and retrieval dates |

## Reproducing the results

Requirements: R 4.6 or later, Quarto, and a free FRED API key. Package versions are pinned with `renv`.

1. Put the FRED API key in your user-level `.Renviron` as `FRED_API_KEY=<key>`, open the project in R from its root folder, and run `renv::restore()`.
2. Get the raw data listed in `data/raw/SOURCE.md`. The raw files are not included in the repository (several sources do not permit redistribution). `Rscript scripts/pull_raw.R` downloads the FRED/ALFRED and Yahoo Finance series; the FRED-QD vintages, the Gürkaynak–Sack–Wright yields, the World Bank commodity file and the survey files are downloaded separately, from the addresses given in `data/raw/SOURCE.md`.
3. Run the stages from the project root:

```
Rscript scripts/build_data.R       # data and validation report
Rscript scripts/estimate.R         # estimation, one model per COVID treatment (a few hours)
Rscript scripts/run_engine.R       # baseline, scenarios, impulse responses
Rscript scripts/export.R           # the website's bundle
Rscript scripts/publish_bundle.R   # copy it into app/
Rscript site/build_site.R          # the static website's data and libraries
```

Then view the website locally: serve the `site/` folder over http (for example `Rscript -e "httpuv::runStaticServer('site')"`; browsers do not load the data files from a `file://` address), or start the Shiny app from R with `shiny::runApp("app")`. Each stage rebuilds only what changed; results are fully determined by `config/run.yml`, its seeds and the raw files.

Tests: `Rscript scripts/run_tests.R`. Check reports: `scripts/check_*.R` (each documents its stages); their results, from the development runs on data through 2026Q2, are in `reports/results/`, readable on GitHub.

The repository keeps the checks that support the data and modelling choices and those that test whether the model's results are reasonable. Further tests and checks were run during development and are not included here: of the code itself (the estimator against exact posterior moments and the `BVAR` package, the scenario engine's exact properties, the website against the model's outputs), of the computational set-up, and of the cost of running the model at full size. They are available on request.

## Data and licences

The sources are US federal agencies (BEA, BLS, Census Bureau, Federal Reserve Board) via FRED and FRED-QD, Moody's, the University of Michigan Surveys of Consumers, S&P Dow Jones Indices, the World Bank and Yahoo Finance (`docs/data.md`). The website shows every series' history with its source credited, but leaves the observed history of the Moody's yields, the S&P 500 and its realized volatility out of its tables and downloads, under their publishers' terms (`docs/data.md`, section 6).

This product uses the FRED® API but is not endorsed or certified by the Federal Reserve Bank of St. Louis.

The code and documents of this repository are released under the MIT licence (`LICENSE`). The licence does not cover the data, which remain under their publishers' terms, or the third-party files in `site/vendor/`, which keep their own licences (listed in `site/vendor/NOTICE.txt`).

## Known limitations

See `docs/methodology.md`, section 16: roots slightly above 1 in every draw and long-run drift of a model in levels, the influence of the COVID treatment, scenarios as conditional forecasts, the identification of impulse responses, and the website's speed (the overview's charts take several seconds to appear).

## Reference

Crump, R. K., S. Eusepi, D. Giannone, E. Qian and A. Sbordone (2021). "A Large Bayesian VAR of the United States Economy." Federal Reserve Bank of New York Staff Report 976.
