# build_example_scenarios.R: write the example scenario files in scenarios/
# from published bank outlooks and from stated oil-price assumptions.
#
# Run from the project root (so .Rprofile activates renv):
#   Rscript scripts/build_example_scenarios.R
#
# Why it exists: published outlooks are stated as growth rates and year-end
# levels (for example "core PCE inflation of 3.4% at year end"), while a
# scenario file holds levels in model units. This script records how each
# published number is mapped to a cell of a scenario file, so the files can
# be checked and rebuilt when the data change. The mapping is rough by
# design: the scenarios borrow a few headline numbers from an outlook; they
# are not the banks' forecasts, and every cell the outlook does not pin down
# is left to the model.
#
# Conventions used in the mapping:
#   - an annual growth or inflation rate for year Y is read as the change
#     from Y-1 Q4 to Y Q4 (Q4 over Q4), placed in quarter Y Q4;
#   - a December or year-end value is placed in Q4, although the data are
#     quarterly averages;
#   - an annual average (for example an oil price for 2027) is placed in
#     each quarter of that year.
#
# Inputs:  config/run.yml (data output folder, horizon), config/series.csv,
#          <data output folder>/model_data.csv (the last observed quarters
#          are the base of the growth rates), <data output folder>/levels.csv
#          (oil prices in dollars for the oil scenarios).
# Outputs: scenarios/<id>.csv for the scenarios defined below; each file is
#          read back with read_scenario() to check it.

for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)

cfg     <- yaml::read_yaml("config/run.yml")
series  <- read_series_config("config/series.csv")
horizon <- cfg$forecast$horizon

md <- utils::read.csv(file.path(cfg$data$output_dir, "model_data.csv"), check.names = FALSE)
history <- as.matrix(md[, series$short_name])
rownames(history) <- md$quarter
last_obs <- md$quarter[nrow(md)]
quarters <- forecast_quarters(last_obs, horizon)
cat("data end in", last_obs, "\n")

#' Growth rates as a vector over the forecast quarters: value g in quarter q.
at <- function(...) {
  x <- c(...)
  g <- rep(NA_real_, horizon)
  g[match(names(x), quarters)] <- x
  if (anyNA(match(names(x), quarters))) stop("quarter outside the horizon: ", paste(names(x), collapse = ", "))
  g
}

#' Level path that is `value` (natural units) in the listed quarters.
level_at <- function(var, ...) {
  x <- c(...)
  out <- rep(NA_real_, horizon)
  out[match(names(x), quarters)] <- to_model_units(var, unname(x), series)
  out
}

#' Write the fixed cells (a list of paths over the forecast quarters) as a
#' scenario file: only the quarters and variables with a fixed cell.
write_scenario <- function(id, paths) {
  m <- do.call(cbind, paths)
  keep <- rowSums(!is.na(m)) > 0
  df <- data.frame(quarter = quarters[keep], round(m[keep, , drop = FALSE], 4), check.names = FALSE)
  path <- file.path("scenarios", paste0(id, ".csv"))
  # LF line endings on every system, as the repository stores them
  # (.gitattributes), so a rebuild does not change an unchanged file. A
  # binary connection, because on Windows a text connection turns "\n"
  # into "\r\n".
  con <- file(path, "wb")
  utils::write.csv(df, con, row.names = FALSE, na = "")
  close(con)
  cond <- read_scenario(path, series$short_name, last_obs, horizon)
  cat(id, ":", sum(!is.na(cond)), "fixed cells in", sum(keep), "quarters ->", path, "\n")
  print(df, row.names = FALSE)
}

# ---- J.P. Morgan, 2026 midyear outlook ----
# Source: J.P. Morgan, "2026 Economic Outlook [Midyear Update]", 9 July 2026
# (forecasts as of 2 July 2026),
# https://www.jpmorgan.com/insights/markets-and-economy/economy/economic-trends
#   - "U.S. GDP growth of 1.5-2.0% expected for 2026": midpoint 1.75%, read
#     as Q4 over Q4 (the basis is not stated);
#   - "core PCE at year end is now 3.4%": 2026Q4 over 2025Q4;
#   - "2-year and 10-year yields around 4.2% and 4.7% at year end": 2026Q4.
# The outlook's Fed funds view (hold at 3.50-3.75%) has no variable in the
# model, and it gives no unemployment forecast; both are left to the model.
write_scenario("jpm_midyear_2026", list(
  rgdp    = path_from_growth("rgdp", at("2026Q4" = 1.75), "yoy", history, series),
  corepce = path_from_growth("corepce", at("2026Q4" = 3.4), "yoy", history, series),
  gs2     = level_at("gs2", "2026Q4" = 4.2),
  gs10    = level_at("gs10", "2026Q4" = 4.7)
))

# ---- Goldman Sachs, inflation and oil, June 2026 ----
# Source: Yahoo Finance, "Goldman Expects Inflation to Ease in 2027 as AI and
# Energy Pressures Recede", 27 June 2026, reporting Goldman Sachs research,
# https://finance.yahoo.com/economy/articles/goldman-expects-inflation-ease-2027-135200615.html
#   - core PCE 3.2% year on year by December 2026 and 2.2% by end-2027:
#     2026Q4 and 2027Q4, each over the same quarter a year earlier;
#   - core CPI 2.6% by December 2026 and 2.2% by December 2027: same;
#   - oil at $80 a barrel on average in 2026Q4 and $75 on average in 2027.
#     The article does not say Brent or WTI; the values are used for WTI
#     (row 29), whose 2026Q2 average was $95.8.
write_scenario("gs_inflation_oil_2026", list(
  corepce = path_from_growth("corepce", at("2026Q4" = 3.2, "2027Q4" = 2.2), "yoy", history, series),
  corecpi = path_from_growth("corecpi", at("2026Q4" = 2.6, "2027Q4" = 2.2), "yoy", history, series),
  wti     = level_at("wti", "2026Q4" = 80, "2027Q1" = 75, "2027Q2" = 75, "2027Q3" = 75, "2027Q4" = 75)
))

# ---- Oil-price assumptions ----
# Two scenarios that fix only the WTI oil price (row 29) over 2026Q3-2027Q2
# and leave everything else, including oil after 2027Q2, to the model. Both
# use only the data the model has (through 2026Q2):
#   - "current" oil is the 2026Q2 average, the last observed quarter (the
#     price spiked from March 2026);
#   - "normal" oil is the 2025 average of the quarterly averages, the year
#     before the spike;
#   - "high" oil is $110 a barrel, a round number above the 2026 spike's
#     monthly peak ($102 in May 2026).
# The two paths are symmetric: current oil in 2026Q3, halfway (in logs, i.e.
# the geometric mean) to the end level in 2026Q4, the end level in 2027Q1
# and 2027Q2.
# Oil rises: the end level is high oil.
# Oil normalizes: the end level is normal oil.
lv <- utils::read.csv(file.path(cfg$data$output_dir, "levels.csv"), check.names = FALSE)
wti_current <- lv$wti[lv$quarter == last_obs]
wti_normal  <- mean(lv$wti[lv$quarter %in% c("2025Q1", "2025Q2", "2025Q3", "2025Q4")])
wti_high    <- 110
#' Oil path: current in 2026Q3, geometric halfway in 2026Q4, `end` after.
oil_path <- function(end) {
  level_at("wti", "2026Q3" = wti_current, "2026Q4" = sqrt(wti_current * end), "2027Q1" = end,
           "2027Q2" = end)
}
cat(sprintf("WTI ($/barrel): current (%s) %.2f; normal (2025 average) %.2f, halfway %.2f; high %.2f, halfway %.2f\n",
            last_obs, wti_current, wti_normal, sqrt(wti_current * wti_normal), wti_high,
            sqrt(wti_current * wti_high)))
write_scenario("oil_rises", list(wti = oil_path(wti_high)))
write_scenario("oil_normalizes", list(wti = oil_path(wti_normal)))
