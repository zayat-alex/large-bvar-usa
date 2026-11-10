# Raw data: sources and retrieval

Every file under `data/raw/` is immutable: never edit or overwrite one. New downloads are added as new files, and each is recorded here with its source, retrieval method and date. Data files are not committed to git (third-party data); this manifest is.

## `fred_qd_vintages/`: FRED-QD monthly vintages

All files are named `FRED-QD_<yyyy>m<mm>.csv`, where `<yyyy>m<mm>` is the vintage month. The vintages come from two separate retrievals.

**1. Vintages 2018-05 to 2025-12 (92 files): historical archive**
- Source: one zip archive from the St. Louis Fed FRED-MD/FRED-QD page:
  https://www.stlouisfed.org/-/media/project/frbstl/stlouisfed/research/fred-md/historical-vintages-of-fred-qd-2018-05-to-2025-12.zip
- Method: downloaded by hand from the page, then unzipped. The zip is kept unchanged as `data/raw/FRED-QD-vintages.zip`; the 92 stored files are byte-identical to its contents (checked 2026-11-04; three renamed, see below).
- Retrieval date: 2026-08-13 (inferred from the zip's file modification time; not recorded at download).
- File names are as in the archive, except three: the archive's `FRED-QD_2019m1.csv`, `FRED-QD_2019m2.csv` and `FRED-QD_2019m3.csv` (no leading zero) were renamed to `FRED-QD_2019m01.csv` to `FRED-QD_2019m03.csv` on 2026-10-16. Contents are unchanged; the zip keeps the original names.

**2. Vintages 2026-01 to 2026-09 (9 files): current-year monthly files**
- Source: the individual monthly vintage CSVs published on the FRED-MD/FRED-QD page (not part of the archive).
- Method: downloaded by hand, one file per month, and renamed to the project convention. Contents are unchanged.
- Retrieval date: 2026-10-16.
- Original names: `2026-<mm>-QD.csv` for 01 to 07 and 09. The August file was published as `2026-rev-08-QD.csv`, a revised August vintage; no unrevised August 2026 file was retrieved. It is stored as `FRED-QD_2026m08.csv`.

Layout (row structure checked on the 2026-08 file; all 2026 files and the 2025-12 file have the same 246 column names): row 1 column names (`sasdate` then 245 series), row 2 `factors` flags, row 3 `transform` codes, then one row per quarter. `sasdate` stamps a quarter with the first day of its last month (`3/1/1959` = 1959:Q1). The project does not use QD's transform codes (docs/data.md, section 1).

**3. Vintage 2026-10 (1 file)**
- Source: `https://www.stlouisfed.org/-/media/project/frbstl/stlouisfed/research/fred-md/quarterly/2026-10-qd.csv`, the file linked as both `2026-10.csv` and `current.csv` on the FRED-MD/FRED-QD page.
- Method: downloaded with curl on 2026-11-10 and stored as `FRED-QD_2026m10.csv` (SHA-256 `be618c35f2d321ece4d3d77aa46f4af16a16718b912377813ec3f596cbcd49cb`). Contents unchanged; the same 246 column names as the 2026-09 vintage.
- Built between 2026-10-29 and 2026-11-05: its 2026:Q3 real GDP is the advance estimate first published on 2026-10-29, and the series of the BLS productivity release published on 2026-11-05 (e.g. OPHNFB, ULCNFB) have no 2026:Q3 value. All 24 QD series the model uses have a 2026:Q3 value.

Last quarter in each 2026 vintage: 2026-01 to 2026-03 end in 2025:Q4; 2026-04 to 2026-06 in 2026:Q1; 2026-07 to 2026-09 in 2026:Q2; 2026-10 in 2026:Q3.

## `alfred/`: ALFRED vintage histories

- Source: FRED API (ALFRED real-time data), through `fredr::fredr_series_observations()` with `realtime_start = 2018-05-01` and `realtime_end` = the retrieval date (the date at the end of each file name). 2018-05-01 is the date of the first stored QD vintage.
- Method: `scripts/pull_raw.R` with `config/raw_pulls.yml`. Native frequency, no aggregation by FRED (the API aggregates only for a single real-time date). Monthly series are averaged to quarterly in project code.
- Retrieval dates: 2026-10-17 (first seven series), 2026-10-19 (`B823RA3Q086SBEA`, `WTISPLC`; the WTISPLC pull ran just after midnight local time on 2026-10-20, but FRED's date, US Central, was 2026-10-19, which is what the file name records). Second retrieval, all eight series: 2026-11-09 (FRED's date; local date 2026-11-10), files `*_to_2026-11-09.csv` and `fred_metadata_2026-11-09.csv`, used from the 2026:Q3 data update (as-of date 2026-11-09).
- File layout: one row per observation and real-time period. Columns `series_id`, `date` (first day of the month or quarter the value refers to), `value` (empty = no value published), `realtime_start`, `realtime_end` (the value was the published one over this period, inclusive). A `realtime_start` of 2018-05-01 means "already published on that date"; a `realtime_end` equal to the retrieval date means "still current on the retrieval date". The series as of date d is the set of rows with `realtime_start <= d <= realtime_end`.

| File | Native frequency | Rows | Observations | Real-time periods | Use |
|---|---|---|---|---|---|
| `B823RA3Q086SBEA_vintages_2018-05-01_to_2026-10-19.csv` | quarterly | 1395 | 1947-01 to 2026-04 | 99 | row 8 (docs/data.md 3.5) |
| `FGCEC1_vintages_2018-05-01_to_2026-10-17.csv` | quarterly | 749 | 1947-01 to 2026-04 (values from 1999 at most; from 2007 since 2023-09-28) | 100 | row 8 index check |
| `HCOMPBS_vintages_2018-05-01_to_2026-10-17.csv` | quarterly | 3558 | 1947-01 to 2026-04 | 68 | row 14 |
| `GS2_vintages_2018-05-01_to_2026-10-17.csv` | monthly | 604 | 1976-06 to 2026-09 | 102 | row 22 |
| `GDPC1_vintages_2018-05-01_to_2026-10-17.csv` | quarterly | 1397 | 1947-01 to 2026-04 | 100 | vintage-date check |
| `PCEPILFE_vintages_2018-05-01_to_2026-10-17.csv` | monthly | 3710 | 1959-01 to 2026-08 | 101 | vintage-date check |
| `GS10_vintages_2018-05-01_to_2026-10-17.csv` | monthly | 882 | 1953-04 to 2026-09 | 102 | aggregation check |
| `WTISPLC_vintages_2018-05-01_to_2026-10-19.csv` | monthly | 984 | 1946-01 to 2026-09 | 116 | row 29 (docs/data.md 3.2) |
| `MCOILWTICO_vintages_2018-05-01_to_2026-10-17.csv` | monthly | 502 | 1986-01 to 2026-09 | 115 | not used (former oil check; docs/data.md 3.2) |
| `B823RA3Q086SBEA_vintages_2018-05-01_to_2026-11-09.csv` | quarterly | 1396 | 1947-01 to 2026-07 | 100 | row 8 |
| `FGCEC1_vintages_2018-05-01_to_2026-11-09.csv` | quarterly | 750 | 1947-01 to 2026-07 | 101 | row 8 index check |
| `HCOMPBS_vintages_2018-05-01_to_2026-11-09.csv` | quarterly | 3581 | 1947-01 to 2026-07 | 69 | row 14 |
| `GS2_vintages_2018-05-01_to_2026-11-09.csv` | monthly | 605 | 1976-06 to 2026-10 | 103 | row 22 |
| `GDPC1_vintages_2018-05-01_to_2026-11-09.csv` | quarterly | 1398 | 1947-01 to 2026-07 | 101 | vintage-date check |
| `PCEPILFE_vintages_2018-05-01_to_2026-11-09.csv` | monthly | 3713 | 1959-01 to 2026-09 | 102 | vintage-date check |
| `GS10_vintages_2018-05-01_to_2026-11-09.csv` | monthly | 883 | 1953-04 to 2026-10 | 103 | aggregation check |
| `WTISPLC_vintages_2018-05-01_to_2026-11-09.csv` | monthly | 985 | 1946-01 to 2026-10 | 117 | row 29 |

- `fred_metadata_2026-10-17.csv`: current metadata of the first seven series (title, units, frequency, seasonal adjustment, observation range, `last_updated`), as of the retrieval date.
- `fred_metadata_B823RA3Q086SBEA_2026-10-19.csv`: the same for `B823RA3Q086SBEA` (Index 2017=100, seasonally adjusted, last updated 2026-09-30).
- `fred_metadata_WTISPLC_2026-10-19.csv`: the same for `WTISPLC` (Dollars per Barrel, not seasonally adjusted, last updated 2026-10-07).

## `gsw/`: Gürkaynak–Sack–Wright Treasury yield curve

- File: `feds200628_2026-10-20.csv` (SHA-256 `7286435090933f0df1dcaeb15629f3a9129c083baee4c42ae5df03d0bb89335c`).
- Source: Federal Reserve Board, "The U.S. Treasury Yield Curve: 1961 to the Present" (Gürkaynak, Sack and Wright, FEDS 2006-28), `https://www.federalreserve.gov/data/yield-curve-tables/feds200628.csv`. The file states it is a staff research product, not an official statistical release.
- Method: downloaded with curl on 2026-10-20 (local date). Contents unchanged.
- Content: daily, 1961-06-14 to 2026-10-09; zero-coupon yields (`SVENYxx`), par yields (`SVENPYxx`), forward rates and the Svensson curve parameters. A 9-line text header precedes the column names.
- Use: row 22 before 1976:Q3 (docs/data.md 3.6); `SVENPY02` is the 2-year par yield.

## `yahoo/`: Yahoo Finance daily closes

- File: `GSPC_daily_2026-10-19.csv` (SHA-256 `9e9761318e0f56324e09d0fe480cc5d61dec1a56e66841e981cbe1e4846083e2`).
- Source: Yahoo Finance, symbol `^GSPC` (S&P 500 price index, not total return), through `quantmod::getSymbols(src = "yahoo")`.
- Method: `scripts/pull_raw.R yahoo` with `config/raw_pulls.yml` (requested from 1927-01-01).
- Retrieval date: 2026-10-19.
- Content: 24,816 trading days, 1927-12-30 to 2026-10-16. Columns `date`, `open`, `high`, `low`, `close`, `volume`, `adjusted`, as returned by Yahoo. No missing or non-positive closes.
- Second file: `GSPC_daily_2026-11-09.csv` (SHA-256 `e336a9a28ffdb726f99809e15d48aaa2dfd763ef7922448854d1d5b708169a50`), retrieved 2026-11-09 (FRED's date in the same run; local date 2026-11-10) by the same method: 24,831 trading days, 1927-12-30 to 2026-11-06, no missing closes. Used from the 2026:Q3 data update.
- Checked on retrieval: two calendar gaps longer than 5 days, both real market closures (1933-03-03 to 1933-03-15, bank holiday; 2001-09-10 to 2001-09-17, after 9/11). There are no Saturday sessions before 1953, although the NYSE traded on Saturdays then; this does not affect the sample from 1973.

## `worldbank/`: World Bank commodity prices (Pink Sheet)

The World Bank replaces the monthly workbook with each release and offers no API or archive of past releases; each yearly document folder on thedocs.worldbank.org keeps only its latest workbook. The series is treated as unrevised (docs/data.md 3.3). The file used in the run is set in `config/run.yml` (`worldbank_file`).

**1. January 2026 release (stored, not used)**
- File: `CMO-Historical-Data-Monthly_2026-10-19.xlsx` (SHA-256 `83a320affafc83512d4af4c7198af6ce086326c8b6923d9f64dd905172f86290`).
- Source: the workbook in the 2025 document folder, `https://thedocs.worldbank.org/en/doc/18675f1d1639c7a34d463f59263ba0a2-0050012025/related/CMO-Historical-Data-Monthly.xlsx` (server last-modified 2026-01-14; same file size). Downloaded by hand.
- Retrieval date: 2026-10-19 (from the file modification time). Original name `CMO-Historical-Data-Monthly.xlsx`; renamed with the download date. Contents unchanged.
- Content as found: sheets `AFOSHEET`, `Monthly Prices`, `Monthly Indices`, `Description`, `Index Weights`. Sheet `Monthly Indices` is headed "monthly indices based on nominal US dollars, 2010=100, 1960 to present" and "Updated on January 06, 2025" (a wrong year: the data run to 2025M12); monthly rows run from 1960M01 to 2025M12, at full precision. The non-energy index is the column headed `Non-energy **`.

**2. October 2026 release (used)**
- File: `CMO-Historical-Data-Monthly_release-2026-10-02_retrieved-2026-10-19.xlsx` (SHA-256 `ea1c350827878ea3bbe30e3cda16a13fd3bd5b409b8458940dc94a36b5a33154`).
- Source: the current workbook linked from https://www.worldbank.org/en/research/commodity-markets, `https://thedocs.worldbank.org/en/doc/74e8be41ceb20fa0da750cda2f6b9e4e-0050012026/related/CMO-Historical-Data-Monthly.xlsx` (server last-modified 2026-10-02).
- Method: downloaded with curl on 2026-10-19; renamed with release and retrieval dates. Contents unchanged. A copy downloaded by hand the same day was identical (same SHA-256) and was deleted.
- Content as found: sheets `Mismatch Details`, `Monthly Prices`, `Monthly Indices`, `Description`, `Index Weights`. `Monthly Indices` is headed "Updated on October 02, 2026", rows 1960M01 to 2026M09. Compared with the January 2026 release: the title row is gone (headers one row higher), index values are stored with one decimal instead of full precision, and a `Mismatch Details` sheet (an internal comparison of Barley and Sorghum cells from 2020M09, one row apart) has replaced `AFOSHEET`.
- Non-energy index compared with the January 2026 release over 1960M01–2025M12 (792 months): 787 months equal the January values rounded to one decimal; 2025M08 to 2025M12 differ by 0.3% to 0.5% (revisions).

**3. November 2026 release (used from the 2026:Q3 data update)**
- File: `CMO-Historical-Data-Monthly_release-2026-11-03_retrieved-2026-11-10.xlsx` (SHA-256 `3c8d9d8dbb48342fdbf24777c39a8ce8964cc2fd0958f72f6779e6e4ad99cde2`).
- Source: the same address as the October release, now linked as "November 2026" from https://www.worldbank.org/en/research/commodity-markets (server last-modified 2026-11-03).
- Method: downloaded with curl on 2026-11-10; renamed with release and retrieval dates. Contents unchanged.
- Content as found: the same sheets as the October release; `Monthly Indices` is headed "Updated on November 03, 2026", rows 1960M01 to 2026M10. The non-energy index equals the October release in all 801 months they share; the new release adds 2026M10.

## `spf/`: Survey of Professional Forecasters, median forecasts

- Files: `medianLevel_2026-11-01.xlsx` (SHA-256 `d296e737331c406f7680414504a0d005472743709387f8fda6d969a1d2200977`) and `medianGrowth_2026-11-01.xlsx` (SHA-256 `fdcf9ca0299de022c58f8fc35fd20761d8e52b45e038781ae3c4826d1a3cbb5c`).
- Source: Federal Reserve Bank of Philadelphia, Survey of Professional Forecasters, median forecasts for all variables, `https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/median-forecasts` (files `historical-data/medianLevel.xlsx` and `historical-data/medianGrowth.xlsx`). File layout: the SPF documentation (`spf-documentation.pdf`, same site).
- Method: downloaded with curl on 2026-11-01; renamed with the retrieval date. Contents unchanged.
- Content as found: `medianLevel` has 58 sheets (one per variable or derived series, e.g. `RGDP`, `UNEMP`, `EMP`, `CPI`, `CORECPI`, `PCE`, `COREPCE`, `TBOND`, `HOUSING`, `INDPROD`); `medianGrowth` has 13 (annualized growth of the level series). One row per survey, 1968Q4 to 2026Q3 (latest survey released 2026-08-14). Columns `<VAR>1` to `<VAR>6` are the quarter before the survey, the survey quarter and the next four quarters; `<VAR>A` onward are annual values. The site notes that historical Moody's Aaa and Baa yields are proprietary and excluded from the levels file.
- Use: comparison of the model's forecasts with the latest survey median.
