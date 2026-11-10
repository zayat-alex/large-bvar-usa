# Data

The model uses a quarterly dataset of 31 variables of the US economy, following Table 1 of Crump, Eusepi, Giannone, Qian and Sbordone, "A Large Bayesian VAR of the United States Economy" (Federal Reserve Bank of New York Staff Report 976, 2021). This document defines each series, where it comes from, how it is constructed, and how it is checked. Where a series differs from the paper's, the reason and the supporting evidence are given (section 3, summarized in section 5).

The sources and retrieval dates of every raw file are listed in `data/raw/SOURCE.md`. The code that builds the dataset is in `R/data_*.R` and runs as the first stage of the `targets` pipeline (`_targets.R`).

---

## 1. Principles

- **Frequency:** quarterly. Every table is indexed by quarter, stored as the first day of the quarter (an R `Date`, e.g. `1973-01-01` for 1973:Q1) and printed as `1973Q1`. FRED-QD stamps a quarter with the first day of its last month (`3/1/1973` = 1973:Q1), so dates are converted on load.
- **Sample start:** 1973:Q1. The sample end is a configuration parameter; by default it is the last quarter in which all series are available.
- **Data sources:**
  - **QD CSV:** FRED-QD files (McCracken and Ng), already quarterly (monthly series are quarterly averages), in levels. QD's own transformation codes are not applied.
  - **FRED API:** series not in QD, taken from ALFRED (FRED's archive of past vintages) at native frequency, as the full vintage history from 2018-05-01 (the first stored QD vintage). Monthly series are averaged to quarterly in project code, under the whole-quarter rule.
  - **Other:** a World Bank spreadsheet downloaded by hand (row 30) and daily S&P 500 closes from Yahoo Finance (row 31).
- **Aggregation rule:** a quarter is used only if the whole quarter's data are available. For monthly series averaged in project code, a quarter with fewer than 3 months is dropped (in practice, the last quarter).
- **Transformations:** applied to levels after assembly. `log` means `100 * ln(x)`, a level measured in log points, not a growth rate. `raw` means untransformed. There is no differencing: the model is a VAR in levels.
- **Prior mean:** the Minnesota prior is centred on a random walk (RW) for every variable except realized volatility, which is centred on white noise (WN). Listed here because it is set per variable.
- **Vintage consistency:** the QD CSV is a specific vintage, and API series are taken from ALFRED as of a date set in the run configuration (`config/run.yml`, `api_realtime_date`). The date depends on the use:
  - **Published forecasts:** the retrieval date, so the forecast uses every release available when it is made. The QD vintage is then the latest published one, which may predate some releases and revisions that reach the API series; this is the usual compromise of a live forecast.
  - **Backtests and pseudo-real-time exercises:** the QD vintage date, so base years and revisions are consistent across sources and the data match what was known at the forecast origin.

## 2. Series table

Source key: **QD** = QD CSV column; **API** = FRED API; **CALC** = constructed; **MANUAL** = file downloaded by hand. IDs are the exact column names or series IDs used.

| # | Short name | Variable | Source | ID | Transform | Prior | Notes |
|---|---|---|---|---|---|---|---|
| 1 | rgdp | Real GDP | QD | GDPC1 | log | RW | |
| 2 | rpce | Real personal consumption expenditures | QD | PCECC96 | log | RW | |
| 3 | rresinv | Real private residential fixed investment | QD | PRFIx | log | RW | |
| 4 | rnresinv | Real private nonresidential fixed investment | QD | PNFIx | log | RW | |
| 5 | rexp | Real exports of goods & services | QD | EXPGSC1 | log | RW | |
| 6 | rimp | Real imports of goods & services | QD | IMPGSC1 | log | RW | |
| 7 | rgov | Real govt. consumption expenditures & gross investment | QD | GCEC1 | log | RW | |
| 8 | rfedgov | Real federal consumption expenditures & gross investment (chain-type quantity index) | API | B823RA3Q086SBEA | log | RW | Quarterly, no aggregation. See 3.5 |
| 9 | gdpdefl | GDP chain price index | QD | GDPCTPI | log | RW | |
| 10 | pce | PCE chain price index | QD | PCECTPI | log | RW | |
| 11 | corepce | PCE less food & energy chain price index | QD | PCEPILFE | log | RW | |
| 12 | cpi | CPI-U, all items | QD | CPIAUCSL | log | RW | |
| 13 | corecpi | CPI-U, all items less food & energy | QD | CPILFESL | log | RW | |
| 14 | comphr | Business sector: compensation per hour (nominal) | API | HCOMPBS | log | RW | Quarterly, no aggregation |
| 15 | payems | All employees, total nonfarm | QD | PAYEMS | log | RW | |
| 16 | urate | Civilian unemployment rate, 16+ | QD | UNRATE | raw | RW | |
| 17 | indpro | Industrial production index | QD | INDPRO | log | RW | |
| 18 | caputil | Capacity utilization, manufacturing (SIC) | QD | CUMFNS | raw | RW | |
| 19 | houst | Housing starts | QD | HOUST | log | RW | |
| 20 | rdpi | Real disposable personal income | QD | DPIC96 | log | RW | |
| 21 | umich | U. Michigan consumer sentiment | QD | UMCSENTx | raw | RW | |
| 22 | gs2 | 2-year Treasury yield, constant maturity | API | GS2 | raw | RW | Monthly: quarterly average. Before 1976:Q3, GSW 2-year par yield; see 3.6 |
| 23 | gs5 | 5-year Treasury yield, constant maturity | QD | GS5 | raw | RW | |
| 24 | gs10 | 10-year Treasury yield, constant maturity | QD | GS10 | raw | RW | |
| 25 | aaa | Moody's seasoned Aaa corporate yield | QD | AAA | raw | RW | |
| 26 | baa | Moody's seasoned Baa corporate yield | QD | BAA | raw | RW | |
| 27 | usd | Nominal trade-weighted US$ index, advanced foreign economies (extended back by QD with the major-currencies index) | QD | TWEXAFEGSMTHx | log | RW | See 3.1 |
| 28 | sp500 | S&P 500 composite | QD | `S&P 500` | log | RW | Column header `S&P 500` in the QD CSV |
| 29 | wti | Spot oil price, WTI (nominal $/barrel) | API | WTISPLC | log | RW | Monthly: quarterly average. See 3.2 |
| 30 | gsci_nonenergy | Non-energy commodity prices (proxy for S&P GSCI non-energy) | MANUAL | World Bank Pink Sheet non-energy price index | log | RW | **Proxy.** See 3.3 |
| 31 | rvol | S&P 500 realized volatility (annualized, %) | CALC | Yahoo `^GSPC` daily closes | raw | WN | See 3.4 |

Counts: 22 log, 9 raw. 25 rows come straight from QD CSV columns, rows 8, 14, 22 and 29 come from the API, row 30 is a manual download and row 31 is computed from daily Yahoo data.

## 3. Construction notes

### 3.1 Dollar index (row 27)
QD's `TWEXAFEGSMTHx` is the nominal advanced foreign economies dollar index (FRED TWEXAFEGSMTH, quarterly average), extended back to 1973 by QD with the discontinued major-currencies index TWEXMMTH, the index the paper uses. It is used as is, from the same QD vintage as the other columns; the project does no splicing of its own. Transform: `100 * ln`.

What the check shows (2026m09 build): up to 2005:Q4 the series is TWEXMMTH rescaled (constant log ratio); from 2006:Q1, where the AFE index starts, it is the AFE index itself. The two baskets differ slightly: over 2006–2019 quarterly growth differs by up to 0.77 pp (2008:Q4) and the log level ratio drifts by about 1.5 log points; the correlation of quarterly growth over 1973–2019 is 0.9997.

Check: comparison with TWEXMMTH from `FRED-QD_2020m03.csv`, the last stored vintage that carries TWEXMMTH (it ends in 2019:Q4), over 1973:Q1–2019:Q4: correlation of `100 * Δln`, standard deviation of the log level ratio, and both series plotted in levels.

### 3.2 Nominal WTI (row 29)
The model uses the nominal WTI spot price in dollars per barrel: FRED `WTISPLC` (monthly, not seasonally adjusted, from 1946), from ALFRED as of the run's as-of date, averaged to quarterly by the project (whole-quarter rule).

QD's OILPRICEx is the same price deflated by core PCE, so QD carries no nominal oil series. On the 2026m09 build, OILPRICEx × PCEPILFE / 100 reproduces WTISPLC to within 0.07% in every quarter since 1973, which confirms WTISPLC is QD's source. FRED's MCOILWTICO (EIA, from 1986) is a different measurement of the same price and differs from WTISPLC by up to 0.75% before 2013; it is not used.

### 3.3 Non-energy commodity price index (row 30)
The paper uses the S&P GSCI Non-Energy index, which is licensed and not publicly available. The proxy is the **World Bank Commodity Price Data ("Pink Sheet") non-energy price index**: monthly, nominal US dollars, covering food, beverages, fertilizers, raw materials and metals (the World Bank file also has an energy index, which is not used).

- **Source:** the Excel file `CMO-Historical-Data-Monthly.xlsx` from the World Bank commodity markets page. It is not on FRED and has no API, so it is downloaded by hand. The file URL changes with each release, so the pipeline reads a local path from the run configuration, and the file name and download date are recorded in `data/raw/SOURCE.md`.
- **Sheet and column:** sheet `Monthly Indices`, column headed `Non-energy **` (nominal USD, 2010 = 100, from 1960M01; checked on the October 2026 release). The header row position differs between releases, so the loader locates the column by its header text. Values are stored with one decimal.
- **Aggregation:** quarterly average of the three monthly values, whole-quarter rule.
- **Transform:** `100 * ln`, prior RW.
- **Checks:**
  1. Coverage reaches 1973:Q1 with no gaps.
  2. The log level shows the known episodes: the 1973–74 boom, the 1980s decline, the 2008 peak and trough, the 2020 dip, the 2021–22 rise.
  3. Quarterly changes (`100 * Δln`) correlate positively with QD's PPI series PPIACO and PPICMM.
- **It is a proxy** and is labelled as such in every output. It is not the same object as the GSCI: the World Bank index is built from spot prices with weights that favour food, beverages and raw materials, while the GSCI is futures-based, production-weighted and includes livestock and precious metals.
- The file is revised only marginally and no real-time vintages exist, so it is treated as unrevised.

### 3.4 S&P 500 realized volatility (row 31)
Not in QD (QD's VIXCLSx is implied volatility, a different quantity). The paper gives no formula, so the following convention is used.

- **Daily data:** Yahoo Finance `^GSPC` daily closing levels (price index, not total return), downloaded with `quantmod::getSymbols("^GSPC", src = "yahoo")` and saved as a dated raw file. The stored file starts on 1927-12-30 and has no gaps in the sample period.
- **Daily return:** `r_d = ln(P_d / P_prev)`, using the previous available trading day's close.
- **Quarterly measure** (annualized, in percent):

```
rvol_q = 100 * sqrt(252 * mean(r_d^2))     # mean over trading days d in quarter q
```

  This is the root mean square of daily returns (zero mean assumed, no demeaning), scaled by 252 so it is comparable to VIX units. Using the mean rather than the sum removes the dependence on the number of trading days in the quarter.
- **Whole-quarter rule:** a quarter is used only if it has ended in the daily data.
- **Transform and prior:** `raw`, WN.
- **Checks:**
  1. Daily returns with |r| > 10% are listed for review. Since 1973 there are four, all genuine: 1987-10-19, 2008-10-13, 2008-10-28 and 2020-03-16.
  2. The quarterly average of daily closes is compared with the QD `S&P 500` column (row 28).
  3. From 1990 on, `rvol` is compared with QD's VIXCLSx (quarterly average of the VIX): a high correlation is expected, with realized volatility usually somewhat below implied. QD dropped VXOCLSx after the 2021m12 vintage, so VIXCLSx is used.
  4. The series shows the spikes of 1987:Q4, 2008:Q4 and 2020:Q1.
- **Scale:** the model is nearly invariant to the scale of a variable, so annualization only affects readability and the units of scenario paths.

### 3.5 Real federal government spending (row 8)
The paper uses real federal consumption expenditures and gross investment in chained dollars (FRED FGCEC1). FRED publishes FGCEC1 only for recent decades (from 2007:Q1 in vintages since 2023-09-28; from 1999 or 2002 in earlier ones), so it cannot cover the sample.

- **Series used:** BEA's chain-type quantity index for the same aggregate, FRED `B823RA3Q086SBEA` (seasonally adjusted, 2017 = 100, from 1947:Q1), from ALFRED as of the run's as-of date.
- **Why it is equivalent:** BEA computes the chained-dollar series as the quantity index times the base-year nominal value divided by 100. In `100 * ln` units the two differ by a constant, which the VAR intercept absorbs; growth rates and dynamics are the same.
- **Check** (validation check 8): over the quarters where both exist, `100 * Δln` of the index and of FGCEC1 agree to rounding and the log level ratio is constant. For vintages as of 2018-05-01, 2023-01-01 and 2026-09-30: largest growth difference 0.0014 pp, standard deviation of the log level ratio 0.0003.

### 3.6 2-year Treasury yield before 1976:Q3 (row 22)
FRED's GS2 (constant maturity, from 1976-06) has its first complete quarter in 1976:Q3. Quarters from the sample start to 1976:Q2 (14 quarters) are filled with the **Gürkaynak–Sack–Wright 2-year par yield** (Federal Reserve, file `feds200628.csv`, column `SVENPY02`, daily from 1961, coupon-equivalent, percent), averaged to quarterly over complete quarters and used **without adjustment**. From 1976:Q3 on, GS2 itself is used.

- **Why this series:** it is market data for the same quantity (a 2-year Treasury par yield, the concept behind constant-maturity yields), so no parameters need fitting.
- **Evidence** (validation check 9, 2026m09 build): over 1976:Q3–2026:Q2, GSW minus GS2 has an RMSE of 0.066 pp (largest 0.33 pp) and quarterly changes correlate at 0.998. In a backcast test (fit on 1981:Q3 onward, predict 1976:Q3–1981:Q2) the unadjusted GSW series misses by 0.136 pp RMSE, against a 1.2 pp standard deviation of GS2's quarterly change in that period.
- **Alternative considered:** a regression of GS2 on GSW, GS1 and GS5 (QD columns) fits better (overlap RMSE 0.040 pp; backcast test 0.056 pp). Over the backfilled quarters it differs from the unadjusted GSW series by at most 0.09 pp. It is not used because it adds fitted parameters for 14 quarters; it is kept in check 9 for comparison.
- Configuration: `data$backfill$gs2` in `config/run.yml`; the backfilled quarters are listed in the validation report.

## 4. Validation checks

Run after assembly; results are written to the validation report (`outputs/data_<vintage>/data_validation.html`).

The checks run with every data build. Numbers quoted in this document come from the builds of the 2026m09 vintage (data through 2026:Q2), on which the model was developed; the update to the 2026m10 vintage (data through 2026:Q3) reran the checks but did not update the numbers quoted here.

1. **Coverage:** each series starts by 1973:Q1; first and last non-missing quarter listed; any gap inside the sample flagged.
2. **Last-quarter completeness:** series missing the final quarter listed; whole-quarter rule applied.
3. **Log safety:** no non-positive values in log-transformed series.
4. **Aggregation check:** the monthly ALFRED GS10 averaged to quarterly with the project's function matches QD's GS10. This validates the averaging used for GS2, WTISPLC and the World Bank index.
5. **Dollar series check** (3.1), with context: the post-2006 differences compared with the series' own growth (standard deviations, regression slope and R-squared, year-on-year differences, level drift against the series' range).
6. **Break detection:** each log level plotted; quarter-to-quarter changes above a few standard deviations flagged (catches splice or base-year errors).
7. **Metadata:** for each API series, title, units, frequency, seasonal adjustment and `last_updated` saved.
8. **Federal spending index check** (3.5): growth of `B823RA3Q086SBEA` against FGCEC1 over their overlap, for several vintages including the run's as-of date: largest difference, correlation, standard deviation of the log level ratio, and both growth paths and their difference plotted.
9. **GS2 backfill check** (3.6): the GSW 2-year par yield, and regression alternatives, against GS2 over their overlap and in a backcast test; backfilled quarters listed; GS2 plotted with the backfill.

## 5. Differences from the paper's dataset

- **Vintage and base year:** the paper's data are 2012-based, with 2005–2018 real-time vintages. Current FRED/QD data use later base years and include later revisions.
- **Federal spending (row 8):** BEA's chain-type quantity index instead of chained dollars, because the chained-dollar series is published only from 2007 (3.5). Same growth rates; the log level differs by a constant.
- **2-year yield (row 22):** 1973:Q1–1976:Q2 filled with the Gürkaynak–Sack–Wright 2-year par yield, because the constant-maturity series starts in 1976:Q3 (3.6).
- **Dollar (row 27):** QD's advanced foreign economies index, extended back by QD with the major-currencies index (3.1). The major-currencies index alone has been discontinued.
- **S&P 500 (row 28):** QD's public-source series.
- **Non-energy commodities (row 30):** World Bank non-energy index as a proxy for the S&P GSCI (3.3).
- **Realized volatility (row 31):** own computation from daily closes (3.4); the paper does not give its formula.

## 6. Redistribution

The raw data are not part of the repository; `data/raw/SOURCE.md` lists each file, its source and its retrieval date, and the pipeline rebuilds everything from them.

The website shows the observed history of every series on its charts, with the sources credited on its About page. Its tables and CSV downloads leave out the observed history of the series whose publishers restrict redistribution, and of the spreads built from them:

| Series | Publisher | Terms |
|---|---|---|
| Aaa and Baa corporate yields (rows 25, 26) | Moody's, via FRED-QD | reproduction without Moody's written consent not permitted |
| S&P 500 (row 28) | S&P Dow Jones Indices, via FRED-QD | reproduction in any form requires prior written permission |
| Realized volatility (row 31) | computed from Yahoo Finance daily closes | Yahoo's terms do not permit reproduction of its data |

The model's forecasts and impulse responses of these series are the project's own results and are included. The other series are public domain (U.S. federal agencies, via FRED; citation requested) or may be used with attribution: University of Michigan consumer sentiment and the WTI spot price (FRED: "copyrighted, citation required") and the World Bank's commodity index (CC BY 4.0). Applications using the FRED API state: "This product uses the FRED® API but is not endorsed or certified by the Federal Reserve Bank of St. Louis."
