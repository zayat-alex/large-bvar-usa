# data_validate.R: the data validation checks (docs/data.md, sections 3 and 4).
#
# Why it exists: every number in the dataset should be traceable and checked
# before any model is estimated on it. Each check_* function returns a list
# with `stats` (a one-row or small data.frame of headline numbers) and `data`
# (the series behind the check), and has a matching plot_* function that draws
# the evidence with base graphics. Checks report; they never change the data.
# Data decisions (replacement series, proxies, release and vintage choices,
# averaging rules) each have a check here, so their evidence is regenerated
# with the data.
#
# Inputs:  the levels / model_data tables, the series list, raw files via the
#          run config.
# Outputs: lists returned to scripts/build_data.R, which saves the plots and
#          a results file read by reports/data_validation.qmd.


# ---- Helpers ----

#' 100 times the first difference of the natural log (quarterly growth, percent).
dlog100 <- function(x) c(NA, 100 * diff(log(x)))

#' Read QD columns used only by the checks.
read_qd_check_columns <- function(cfg, columns, file = cfg$data$qd_vintage_file) {
  read_qd_vintage(raw_path(cfg, file), columns)
}

#' ALFRED series as of the run's as-of date (native frequency).
alfred_series <- function(cfg, id, as_of = cfg$data$api_realtime_date) {
  alfred_as_of(read_alfred_vintages(raw_path(cfg, cfg$data$alfred_files[[id]])), as_of)
}


# ---- 1-3: coverage, last quarter, log safety ----

#' Check 1, coverage: first and last observation and gaps inside the sample.
#'
#' @param levels Levels table (sample quarters only).
#' @param series_list Output of build_series_list() (full histories).
#' @param series Series config.
#' @return list(stats, data): stats has one row per variable with the first
#'   and last quarter available in the source, whether the series starts by
#'   the sample start, and the number of missing quarters inside the sample.
check_coverage <- function(levels, series_list, series) {
  start <- min(levels$quarter)
  stats <- do.call(rbind, lapply(series$short_name, function(nm) {
    x <- series_list[[nm]]
    in_sample <- levels[[nm]]
    first_in  <- which(!is.na(in_sample))[1]
    data.frame(
      series          = nm,
      first_available = format_quarter(min(x$quarter)),
      last_available  = format_quarter(max(x$quarter)),
      starts_by_start = min(x$quarter) <= start,
      missing_in_sample  = sum(is.na(in_sample)),
      gaps_after_start   = sum(is.na(in_sample[first_in:length(in_sample)])),
      stringsAsFactors = FALSE
    )
  }))
  list(stats = stats, data = NULL)
}

#' Check 2, last-quarter completeness: which series end before or after the
#' sample end (the sample ends at the last quarter all series share).
check_last_quarter <- function(levels, series_list) {
  end <- max(levels$quarter)
  last <- vapply(series_list, function(x) as.numeric(max(x$quarter)), numeric(1))
  stats <- data.frame(series = names(series_list),
                      last_available = format_quarter(as.Date(last, origin = "1970-01-01")),
                      relative_to_sample_end = ifelse(last > end, "later", ifelse(last < end, "EARLIER", "same")),
                      stringsAsFactors = FALSE)
  list(stats = stats, data = data.frame(sample_end = format_quarter(end)))
}

#' Check 3, log safety: smallest value of each log series (must be > 0).
check_log_safety <- function(levels, series) {
  logs <- series$short_name[series$transform == "log"]
  stats <- data.frame(series = logs,
                      min_value = vapply(logs, function(nm) min(levels[[nm]], na.rm = TRUE), numeric(1)),
                      n_nonpositive = vapply(logs, function(nm) sum(levels[[nm]] <= 0, na.rm = TRUE), numeric(1)),
                      row.names = NULL)
  list(stats = stats, data = NULL)
}


# ---- 4: aggregation (project-side monthly averaging) ----

#' Check 4: GS10 averaged by the project vs QD's GS10.
#'
#' Validates monthly_to_quarterly(), which is also used for GS2 (row 22),
#' WTISPLC (row 29) and the World Bank index (row 30).
#' Both are quarterly averages of the same monthly series, so they should
#' agree up to QD's rounding.
check_aggregation_gs10 <- function(cfg) {
  ours <- monthly_to_quarterly(alfred_series(cfg, "GS10"))
  qd   <- read_qd_check_columns(cfg, "GS10")
  d <- merge(ours, data.frame(quarter = qd$quarter, qd = qd$GS10), by = "quarter")
  d <- d[!is.na(d$qd), ]
  d$diff <- d$value - d$qd
  list(stats = data.frame(quarters = nrow(d), first = format_quarter(min(d$quarter)),
                          last = format_quarter(max(d$quarter)),
                          max_abs_diff_pp = max(abs(d$diff)), mean_diff_pp = mean(d$diff)),
       data = d)
}

plot_aggregation_gs10 <- function(chk) {
  d <- chk$data
  par(mfrow = c(2, 1), mar = c(3, 5, 3, 1))
  plot(d$quarter, d$qd, type = "l", lwd = 4, col = "grey70", xlab = "", ylab = "% p.a.",
       main = "GS10 quarterly average: QD (grey) vs project averaging of ALFRED monthly (black)")
  lines(d$quarter, d$value, lwd = 1.2)
  plot(d$quarter, d$diff, type = "h", xlab = "", ylab = "pp", main = "Difference, project minus QD")
  abline(h = 0, lty = 3)
}


# ---- 5: dollar index ----

#' Check 5 (docs/data.md 3.1): QD's TWEXAFEGSMTHx (run vintage) vs TWEXMMTH
#' (last QD vintage carrying it), 1973:Q1-2019:Q4.
#'
#' Up to 2005 QD's series is TWEXMMTH rescaled; from 2006 it is the AFE index,
#' a slightly different basket. `context` sizes the post-2006 differences
#' against the series' own movements, so a reader can judge whether they matter.
#'
#' @param cfg Parsed run config.
#' @param split First quarter of the AFE-based segment.
#' @return list(stats, context, data).
check_dollar <- function(cfg, split = as.Date("2006-04-01")) {
  afe <- read_qd_check_columns(cfg, "TWEXAFEGSMTHx")
  mm  <- read_qd_check_columns(cfg, "TWEXMMTH", file = cfg$data$qd_twexmmth_file)
  d <- merge(afe, mm, by = "quarter")
  d <- d[!is.na(d$TWEXAFEGSMTHx) & !is.na(d$TWEXMMTH), ]
  d$g_afe <- dlog100(d$TWEXAFEGSMTHx); d$g_mm <- dlog100(d$TWEXMMTH)
  d$logratio <- 100 * log(d$TWEXAFEGSMTHx / d$TWEXMMTH)
  ok <- !is.na(d$g_afe)

  # Context, AFE segment only (first growth rate fully inside it onward).
  p   <- d[d$quarter >= split, ]
  dd  <- p$g_afe - p$g_mm
  fit <- stats::lm(g_afe ~ g_mm, p)
  yoy <- function(z) 100 * (log(z[-(1:4)]) - log(z[1:(length(z) - 4)]))
  py  <- d[d$quarter >= seq(split, by = "-1 year", length.out = 2)[2], ]   # four quarters of lead-in for YoY
  ydiff <- yoy(py$TWEXAFEGSMTHx) - yoy(py$TWEXMMTH)
  context <- data.frame(
    measure = c("quarters", "sd of q/q growth, AFE (pp)", "sd of q/q growth difference (pp)",
                "ratio: sd difference / sd growth", "mean |q/q growth|, AFE (pp)", "mean |q/q growth difference| (pp)",
                "slope, AFE growth on MM growth", "R-squared, AFE growth on MM growth",
                "sd of YoY growth, AFE (pp)", "max |YoY growth difference| (pp)",
                "level drift over the segment (log points)", "range of 100 x ln(AFE) over the segment"),
    value = c(nrow(p), sd(p$g_afe), sd(dd), sd(dd) / sd(p$g_afe), mean(abs(p$g_afe)), mean(abs(dd)),
              unname(stats::coef(fit)[2]), summary(fit)$r.squared,
              sd(yoy(py$TWEXAFEGSMTHx)), max(abs(ydiff)),
              p$logratio[nrow(p)] - d$logratio[d$quarter == quarter_start(split - 1)], diff(range(100 * log(p$TWEXAFEGSMTHx))))
  )
  list(stats = data.frame(quarters = nrow(d), first = format_quarter(min(d$quarter)),
                          last = format_quarter(max(d$quarter)),
                          corr_growth = cor(d$g_afe[ok], d$g_mm[ok]),
                          max_abs_growth_diff_pp = max(abs(d$g_afe[ok] - d$g_mm[ok])),
                          sd_log_ratio = sd(d$logratio)),
       context = context,
       data = d)
}

plot_dollar <- function(chk) {
  d <- chk$data
  par(mfrow = c(3, 1), mar = c(3, 5, 3, 1))
  plot(d$quarter, d$TWEXMMTH, type = "l", lwd = 4, col = "grey70", log = "y", xlab = "", ylab = "index (log scale)",
       ylim = range(c(d$TWEXMMTH, d$TWEXAFEGSMTHx)),
       main = "Dollar: TWEXMMTH, QD 2020m03 (grey) vs TWEXAFEGSMTHx, run vintage (black)")
  lines(d$quarter, d$TWEXAFEGSMTHx, lwd = 1.2)
  plot(d$quarter, d$g_afe - d$g_mm, type = "h", xlab = "", ylab = "pp", main = "Growth difference (100 x dln), AFE minus MM")
  abline(h = 0, lty = 3)
  plot(d$quarter, d$logratio, type = "l", xlab = "", ylab = "100 x ln(AFE / MM)", main = "Log level ratio")
}


# ---- 6: break detection ----

#' Check 6: flag unusually large quarter-to-quarter changes.
#'
#' For each series, the change is 100*dln for log series and the first
#' difference for raw series; a quarter is flagged when the change is more than
#' `k` standard deviations from the series' mean change. Flags are for review:
#' real events (1980 credit controls, 2008-09, the 2020 pandemic quarters) are
#' expected to appear alongside any splice or base-year error.
check_breaks <- function(model_data, series, k = 4) {
  flags <- do.call(rbind, lapply(seq_len(nrow(series)), function(i) {
    nm <- series$short_name[i]
    ch <- c(NA, diff(model_data[[nm]]))      # model_data is already 100*ln for log series
    z  <- (ch - mean(ch, na.rm = TRUE)) / sd(ch, na.rm = TRUE)
    hit <- which(abs(z) > k)
    if (length(hit) == 0) return(NULL)
    data.frame(series = nm, quarter = format_quarter(model_data$quarter[hit]),
               change = round(ch[hit], 2), z = round(z[hit], 1))
  }))
  if (is.null(flags)) flags <- data.frame(series = character(), quarter = character(), change = numeric(), z = numeric())
  list(stats = data.frame(threshold_sd = k, n_flags = nrow(flags)), data = flags)
}


# ---- 7: metadata ----

#' Check 7: FRED metadata of the API series, as saved at retrieval.
check_metadata <- function(cfg) {
  m <- do.call(rbind, lapply(cfg$data$alfred_metadata_files, function(f) utils::read.csv(raw_path(cfg, f), stringsAsFactors = FALSE)))
  list(stats = m[, c("id", "title", "units", "frequency", "seasonal_adjustment", "last_updated")], data = NULL)
}


# ---- 8: federal spending index vs chained dollars ----

#' Check 8 (docs/data.md 3.5): B823RA3Q086SBEA vs FGCEC1, several vintages.
#'
#' If the index measures the same BEA quantity, growth rates agree to
#' rounding and the log level ratio is constant within each vintage.
check_federal_index <- function(cfg, vintages = c("2018-05-01", "2023-01-01", cfg$data$api_realtime_date)) {
  idx <- read_alfred_vintages(raw_path(cfg, cfg$data$alfred_files[["B823RA3Q086SBEA"]]))
  fg  <- read_alfred_vintages(raw_path(cfg, cfg$data$alfred_files[["FGCEC1"]]))
  per <- lapply(vintages, function(v) {
    d <- merge(alfred_as_of(idx, v), alfred_as_of(fg, v), by = "date", suffixes = c(".idx", ".fg"))
    d$g_idx <- dlog100(d$value.idx); d$g_fg <- dlog100(d$value.fg)
    d$diff <- d$g_idx - d$g_fg
    d$logratio <- 100 * log(d$value.fg / d$value.idx)
    d
  })
  names(per) <- vintages
  stats <- do.call(rbind, lapply(vintages, function(v) {
    d <- per[[v]]; ok <- !is.na(d$diff)
    data.frame(as_of = v, first = format_quarter(min(d$date)), last = format_quarter(max(d$date)),
               quarters = nrow(d), max_abs_growth_diff_pp = max(abs(d$diff[ok])),
               corr_growth = cor(d$g_idx[ok], d$g_fg[ok]), sd_log_ratio = sd(d$logratio))
  }))
  list(stats = stats, data = per)
}

plot_federal_index <- function(chk) {
  per <- chk$data; latest <- per[[length(per)]]; ok <- !is.na(latest$diff)
  par(mfrow = c(2, 1), mar = c(3, 5, 3, 1))
  plot(latest$date[ok], latest$g_fg[ok], type = "l", lwd = 4, col = "grey70", xlab = "", ylab = "100 x dln, % q/q",
       main = paste("Federal spending growth as of", names(per)[length(per)], ": FGCEC1 (grey) vs quantity index (black)"))
  lines(latest$date[ok], latest$g_idx[ok], lwd = 1.2); abline(h = 0, lty = 3)
  cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")[seq_along(per)]
  xr <- as.Date(range(unlist(lapply(per, function(d) as.character(d$date)))))
  m  <- max(0.01, max(abs(unlist(lapply(per, function(d) d$diff))), na.rm = TRUE))
  plot(xr, c(-m, m), type = "n", xlab = "", ylab = "pp", main = "Growth difference, index minus FGCEC1, by vintage")
  for (i in seq_along(per)) lines(per[[i]]$date, per[[i]]$diff, col = cols[i], type = "o", pch = 16, cex = 0.5)
  abline(h = 0, lty = 3)
  legend("topleft", paste("as of", names(per)), col = cols, pch = 16, bty = "n", horiz = TRUE)
}


# ---- 9: GS2 backfill ----

#' Check 9 (docs/data.md 3.6): GS2 against the GSW 2-year par yield used to
#' backfill it, and against regression alternatives.
#'
#' Over the overlap (GS2 from 1976:Q3), each candidate is compared with GS2;
#' then a backcast test fits each candidate on 1981:Q3 onward and predicts the
#' first five years of GS2, mimicking the backfill. The unadjusted GSW series
#' is the one used; the regression on GSW + GS1 + GS5 fits better and is
#' reported as the documented alternative.
#'
#' @param cfg Parsed run config.
#' @param test_split First quarter of the fitting window in the backcast test.
#' @return list(stats = one row per candidate, data = quarterly series).
check_gs2_backfill <- function(cfg, test_split = as.Date("1981-07-01")) {
  b   <- cfg$data$backfill$gs2
  gs2 <- monthly_to_quarterly(alfred_series(cfg, "GS2"))
  gsw <- daily_to_quarterly(read_gsw_daily(raw_path(cfg, b$file), b$column), as.Date(b$complete_through))
  qd  <- read_qd_check_columns(cfg, c("GS1", "GS5"))
  d <- merge(merge(setNames(gsw, c("quarter", "gsw")), qd, by = "quarter"),
             setNames(gs2, c("quarter", "gs2")), by = "quarter", all.x = TRUE)
  d <- d[d$quarter >= as.Date("1973-01-01"), ]
  o <- d[!is.na(d$gs2), ]
  tr <- o[o$quarter >= test_split, ]; te <- o[o$quarter < test_split, ]

  formulas <- list("OLS: GS2 on GSW par" = gs2 ~ gsw,
                   "OLS: GS2 on GS1 + GS5" = gs2 ~ GS1 + GS5,
                   "OLS: GS2 on GSW par + GS1 + GS5" = gs2 ~ gsw + GS1 + GS5)
  row <- function(name, e_all, e_test) {
    data.frame(candidate = name, rmse_overlap_pp = sqrt(mean(e_all^2)), max_abs_overlap_pp = max(abs(e_all)),
               rmse_backcast_test_pp = sqrt(mean(e_test^2)), max_abs_backcast_test_pp = max(abs(e_test)))
  }
  stats <- row("GSW par, unadjusted (used)", o$gs2 - o$gsw, te$gs2 - te$gsw)
  for (nm in names(formulas)) {
    stats <- rbind(stats, row(nm, stats::resid(stats::lm(formulas[[nm]], o)),
                              te$gs2 - stats::predict(stats::lm(formulas[[nm]], tr), te)))
  }
  full <- stats::lm(gs2 ~ gsw + GS1 + GS5, o)
  d$ols_best <- stats::predict(full, d)
  context <- data.frame(
    measure = c("overlap", "correlation of q/q changes, GS2 vs GSW", "sd of GS2 q/q change, overlap to 1985Q4 (pp)",
                "backcast test window", "max |GSW - best OLS| over the backfilled quarters (pp)"),
    value = c(paste(format_quarter(min(o$quarter)), "to", format_quarter(max(o$quarter))),
              format(round(cor(diff(o$gs2), diff(o$gsw)), 4)),
              format(round(sd(diff(o$gs2[o$quarter <= as.Date("1985-10-01")])), 2)),
              paste(format_quarter(min(te$quarter)), "to", format_quarter(max(te$quarter))),
              format(round(max(abs(d$gsw - d$ols_best)[is.na(d$gs2)]), 3))))
  list(stats = stats, context = context, data = d, coef_best = stats::coef(full))
}

plot_gs2_backfill <- function(chk) {
  d <- chk$data; s <- d[d$quarter <= as.Date("1990-10-01"), ]; o <- d[!is.na(d$gs2), ]
  par(mfrow = c(2, 1), mar = c(3, 5, 3, 1))
  plot(s$quarter, s$gsw, type = "l", col = "red", lwd = 2, ylim = range(s[, c("gs2", "gsw", "ols_best")], na.rm = TRUE),
       xlab = "", ylab = "% p.a.", main = "GS2 (black) and the GSW 2-year par yield (red) used before 1976Q3, 1973-1990")
  lines(s$quarter, s$ols_best, col = "blue", lty = 2)
  lines(s$quarter, s$gs2, lwd = 2.5)
  abline(v = min(o$quarter), lty = 3)
  legend("topleft", c("GS2 (FRED, from 1976Q3)", "GSW par 2y, unadjusted (backfill)", "OLS on GSW + GS1 + GS5 (alternative)"),
         col = c("black", "red", "blue"), lty = c(1, 1, 2), lwd = c(2.5, 2, 1), bty = "n", cex = 0.85)
  plot(o$quarter, o$gs2 - o$ols_best, type = "p", pch = 16, cex = 0.45, col = "blue", ylim = c(-0.5, 0.5),
       xlab = "", ylab = "pp", main = "GS2 minus GSW par (red bars) and minus the OLS alternative (blue dots), overlap")
  segments(o$quarter, 0, o$quarter, o$gs2 - o$gsw, col = "red")
  abline(h = 0, lty = 3)
}


# ---- 3.3: World Bank non-energy index ----

#' docs/data.md 3.3 checks: correlation of quarterly changes with QD's PPI
#' series, and the comparison of the stored releases (release choice).
check_worldbank <- function(levels, cfg) {
  qd <- read_qd_check_columns(cfg, c("PPIACO", "PPICMM"))
  d  <- merge(levels[, c("quarter", "gsci_nonenergy")], qd, by = "quarter")
  d$g_wb <- dlog100(d$gsci_nonenergy); d$g_ppiaco <- dlog100(d$PPIACO); d$g_ppicmm <- dlog100(d$PPICMM)
  ok <- !is.na(d$g_wb) & !is.na(d$g_ppiaco) & !is.na(d$g_ppicmm)

  new <- read_worldbank_index(raw_path(cfg, cfg$data$worldbank_file))
  old <- read_worldbank_index(raw_path(cfg, cfg$data$worldbank_previous_file))
  rel <- merge(old, new, by = "date", suffixes = c(".jan2026", ".oct2026"))
  rel$pct_diff <- 100 * (rel$value.oct2026 / rel$value.jan2026 - 1)

  list(stats = data.frame(corr_growth_PPIACO = cor(d$g_wb[ok], d$g_ppiaco[ok]),
                          corr_growth_PPICMM = cor(d$g_wb[ok], d$g_ppicmm[ok]),
                          release_overlap_months = nrow(rel),
                          release_months_beyond_rounding = sum(abs(rel$value.oct2026 - rel$value.jan2026) > 0.05 + 1e-9),
                          release_max_abs_pct_diff = max(abs(rel$pct_diff))),
       data = list(ppi = d, releases = rel))
}

plot_worldbank <- function(chk) {
  d <- chk$data$ppi; rel <- chk$data$releases
  par(mfrow = c(3, 1), mar = c(3, 5, 3, 1))
  plot(d$quarter, 100 * log(d$gsci_nonenergy), type = "l", xlab = "", ylab = "100 x ln",
       main = "World Bank non-energy index (proxy for S&P GSCI non-energy), log level")
  plot(d$quarter, d$g_wb, type = "l", xlab = "", ylab = "100 x dln", ylim = range(c(d$g_wb, d$g_ppicmm), na.rm = TRUE),
       main = "Quarterly growth: World Bank non-energy (black) vs QD PPICMM, metals (red)")
  lines(d$quarter, d$g_ppicmm, col = "red"); abline(h = 0, lty = 3)
  plot(rel$date, rel$pct_diff, type = "h", xlab = "", ylab = "%",
       main = "Release comparison: October 2026 vs January 2026 release, monthly, percent difference")
  abline(h = c(-0.2, 0, 0.2), lty = c(2, 3, 2), col = c("grey50", "black", "grey50"))
}


# ---- 3.4: realized volatility ----

#' docs/data.md 3.4 checks: large daily returns, average close vs QD S&P 500,
#' realized vs implied volatility (VIXCLSx, from 1990).
check_rvol <- function(levels, cfg, threshold = 0.10) {
  daily <- read_yahoo_daily(raw_path(cfg, cfg$data$yahoo_file))
  daily <- daily[daily$date >= min(levels$quarter) - 7, ]
  r <- data.frame(date = daily$date[-1], r = diff(log(daily$close)))
  big <- r[abs(r$r) > threshold, ]

  avg_close <- aggregate(close ~ quarter, transform(daily, quarter = quarter_start(date)), mean)
  qd <- read_qd_check_columns(cfg, c("S&P 500", "VIXCLSx"))
  sp <- merge(avg_close, data.frame(quarter = qd$quarter, qd = qd[["S&P 500"]]), by = "quarter")
  sp <- sp[sp$quarter >= min(levels$quarter) & sp$quarter <= max(levels$quarter) & !is.na(sp$qd), ]
  sp$pct_diff <- 100 * (sp$close / sp$qd - 1)

  v <- merge(levels[, c("quarter", "rvol")], data.frame(quarter = qd$quarter, vix = qd$VIXCLSx), by = "quarter")
  v <- v[v$quarter >= as.Date("1990-01-01") & !is.na(v$vix), ]

  top <- levels[order(-levels$rvol), c("quarter", "rvol")][1:5, ]
  top$quarter <- format_quarter(top$quarter)

  list(stats = data.frame(daily_moves_over_10pct = nrow(big),
                          sp500_max_abs_pct_diff = max(abs(sp$pct_diff)), sp500_mean_pct_diff = mean(sp$pct_diff),
                          corr_rvol_vix_1990 = cor(v$rvol, v$vix), mean_rvol_minus_vix = mean(v$rvol - v$vix)),
       data = list(big_moves = transform(big, r = round(r, 4)), sp500 = sp, vix = v, top = top))
}

plot_rvol <- function(chk) {
  sp <- chk$data$sp500; v <- chk$data$vix
  par(mfrow = c(3, 1), mar = c(3, 5, 3, 1))
  plot(v$quarter, v$vix, type = "l", col = "red", xlab = "", ylab = "% annualized", ylim = range(c(v$vix, v$rvol)),
       main = "Realized volatility from daily ^GSPC (black) vs implied, QD VIXCLSx (red), from 1990")
  lines(v$quarter, v$rvol)
  plot(sp$quarter, sp$pct_diff, type = "h", xlab = "", ylab = "%",
       main = "Quarterly average of Yahoo ^GSPC closes vs QD 'S&P 500', percent difference")
  abline(h = 0, lty = 3)
  plot(chk$data$vix$quarter, chk$data$vix$rvol - chk$data$vix$vix, type = "h", xlab = "", ylab = "pp",
       main = "Realized minus implied volatility")
  abline(h = 0, lty = 3)
}


# ---- Vintage-date match (choice of api_realtime_date) ----

#' Evidence for the ALFRED as-of date: distance between QD's GDPC1 and
#' PCEPILFE and the ALFRED series as of each candidate date.
#'
#' The chosen date should be the first from which both match up to QD's
#' rounding (3 decimals for PCEPILFE).
check_vintage_date <- function(cfg, from = "2026-08-01", to = NULL) {
  gdp <- read_alfred_vintages(raw_path(cfg, cfg$data$alfred_files[["GDPC1"]]))
  pce <- read_alfred_vintages(raw_path(cfg, cfg$data$alfred_files[["PCEPILFE"]]))
  qd  <- read_qd_check_columns(cfg, c("GDPC1", "PCEPILFE"))
  if (is.null(to)) to <- min(max(gdp$realtime_end), max(pce$realtime_end))
  days <- seq(as.Date(from), as.Date(to), by = "day")
  dist <- function(q, col) {
    m <- merge(q, data.frame(quarter = qd$quarter, qd = qd[[col]]), by = "quarter")
    m <- m[!is.na(m$qd), ]
    max(abs(m$value / m$qd - 1)) * 100
  }
  d <- data.frame(as_of = days,
                  gdp_max_abs_pct = vapply(days, function(x) { s <- alfred_as_of(gdp, x); dist(data.frame(quarter = s$date, value = s$value), "GDPC1") }, numeric(1)),
                  pce_max_abs_pct = vapply(days, function(x) dist(monthly_to_quarterly(alfred_as_of(pce, x)), "PCEPILFE"), numeric(1)))
  chosen <- as.Date(cfg$data$api_realtime_date)
  list(stats = data.frame(as_of = format(chosen),
                          gdp_max_abs_pct = d$gdp_max_abs_pct[d$as_of == chosen],
                          pce_max_abs_pct = d$pce_max_abs_pct[d$as_of == chosen],
                          gdp_day_before = d$gdp_max_abs_pct[d$as_of == chosen - 1]),
       data = d)
}

plot_vintage_date <- function(chk, chosen) {
  d <- chk$data
  par(mfrow = c(1, 1), mar = c(4, 5, 3, 1))
  # Exact matches (0) cannot be drawn on a log scale; they are shown at 1e-6.
  y <- pmax(cbind(d$gdp_max_abs_pct, d$pce_max_abs_pct), 1e-6)
  matplot(d$as_of, y, type = "s", lty = c(1, 2), lwd = 2.5, col = c("black", "red"), log = "y", xaxt = "n",
          xlab = "ALFRED as-of date", ylab = "max |ALFRED / QD - 1|, % (log scale)",
          main = "Which ALFRED date reproduces the QD vintage?")
  axis.Date(1, d$as_of)
  abline(v = as.Date(chosen), lty = 3, col = "grey40")
  abline(h = 1e-6, lty = 3, col = "grey70")
  legend("bottomleft", c("GDPC1 (quarterly)", "PCEPILFE (monthly, averaged)", "1e-6 = exact match", "chosen as-of date"),
         col = c("black", "red", "grey70", "grey40"), lty = c(1, 2, 3, 3), lwd = c(2.5, 2.5, 1, 1), bty = "n")
}


# ---- Run all ----

#' Run every check and return the results as one list.
#'
#' @param levels,model_data Data tables for the sample.
#' @param series_list Full-history series from build_series_list().
#' @param series Series config.
#' @param cfg Parsed run config.
run_data_checks <- function(levels, model_data, series_list, series, cfg) {
  list(
    coverage      = check_coverage(levels, series_list, series),
    last_quarter  = check_last_quarter(levels, series_list),
    log_safety    = check_log_safety(levels, series),
    aggregation   = check_aggregation_gs10(cfg),
    dollar        = check_dollar(cfg),
    breaks        = check_breaks(model_data, series),
    metadata      = check_metadata(cfg),
    federal_index = check_federal_index(cfg),
    gs2_backfill  = c(check_gs2_backfill(cfg), list(audit = attr(series_list, "backfill"))),
    worldbank     = check_worldbank(levels, cfg),
    rvol          = check_rvol(levels, cfg),
    vintage_date  = check_vintage_date(cfg)
  )
}

#' Save the plot of each check as a PNG.
#'
#' @param checks Output of run_data_checks().
#' @param dir Output folder for the PNG files.
#' @param chosen_date The run's ALFRED as-of date (marked on the plot).
#' @return Named character vector of file paths, invisibly.
save_check_plots <- function(checks, dir, chosen_date) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  plots <- list(
    aggregation   = function() plot_aggregation_gs10(checks$aggregation),
    dollar        = function() plot_dollar(checks$dollar),
    federal_index = function() plot_federal_index(checks$federal_index),
    gs2_backfill  = function() plot_gs2_backfill(checks$gs2_backfill),
    worldbank     = function() plot_worldbank(checks$worldbank),
    rvol          = function() plot_rvol(checks$rvol),
    vintage_date  = function() plot_vintage_date(checks$vintage_date, chosen_date)
  )
  paths <- character()
  for (nm in names(plots)) {
    p <- file.path(dir, paste0(nm, ".png"))
    grDevices::png(p, width = 1400, height = if (nm == "vintage_date") 700 else 1300, res = 130)
    plots[[nm]]()
    grDevices::dev.off()
    paths[nm] <- p
  }
  invisible(paths)
}


# ---- Exploratory PDF ----

#' One-page-per-six-series PDF of every variable over the sample.
#'
#' For a quick visual scan only: model units (100*ln for log series, natural
#' units for raw series), no transformations beyond those of the docs/data.md.
#'
#' @param model_data Model-units table.
#' @param series Series config (titles and transforms).
#' @param file Output PDF path.
#' @param subtitle Text printed at the bottom of each page (vintage, sample).
#' @return file (the path written).
plot_exploratory_pdf <- function(model_data, series, file, subtitle) {
  grDevices::pdf(file, width = 11, height = 8.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  per_page <- 6
  for (start in seq(1, nrow(series), by = per_page)) {
    par(mfrow = c(3, 2), mar = c(3, 4.5, 3, 1), oma = c(2, 0, 0, 0))
    for (i in start:min(start + per_page - 1, nrow(series))) {
      nm <- series$short_name[i]
      ylab <- if (series$transform[i] == "log") "100 x ln(level)" else "level (natural units)"
      plot(model_data$quarter, model_data[[nm]], type = "l", xlab = "", ylab = ylab,
           main = sprintf("%d. %s: %s", series$row[i], nm, series$variable[i]), cex.main = 0.85)
      grid(col = "grey85")
    }
    mtext(subtitle, side = 1, outer = TRUE, cex = 0.7, line = 0.5)
  }
  file
}
