# settings.R: what the website shows and how it names things.
#
# Presentation only: which scenarios appear and in what colour, the short
# labels and units of the variables, the topic tabs and the default charts
# of the overview. The numbers all come from the bundle (bundle.R); nothing
# here changes a forecast.


#' Scenarios listed in the header when the site opens, in that order.
#'
#' Every scenario in the bundle can be listed; the reader changes the list
#' on the Settings tab. The first one is the default baseline. Which
#' forecast counts as the baseline is a judgement, so the reader can pick
#' any of them. The default baseline is the forecast tilted to the Survey of
#' Professional Forecasters, not the unconditional forecast, which can still
#' be listed. Ids not in the bundle are ignored.
APP_DEFAULT_SCENARIOS <- c("spf_2026q3_2028", "oil_rises", "oil_normalizes")

#' Fixed colours of the published scenarios.
#'
#' The default baseline keeps the site's first blue; the unconditional
#' forecast is dark grey (the observed history is near-black, #1B1F24); the
#' oil scenarios are dark red (higher prices) and dark green (lower prices);
#' the scenarios after J.P. Morgan and Goldman Sachs are drawn in colours of
#' those firms (a teal from J.P. Morgan's website, since its main colour,
#' black, would read as history; Goldman Sachs blue).
APP_SCENARIO_COLOURS <- c(spf_2026q3_2028 = "#0055A4", baseline = "#595F66", oil_rises = "#9E1B1B",
                          oil_normalizes = "#1E6B34", jpm_midyear_2026 = "#1A7B99",
                          gs_inflation_oil_2026 = "#7399C6")

#' Colours of any other scenario, assigned in this order.
#'
#' Scenarios without a fixed colour take the palette's colours in turn,
#' skipping those already fixed above, so a scenario keeps its colour
#' whatever else is switched on. The eight colours (blue, orange, plum,
#' green, ochre, sky, rust, violet) were checked in this order for
#' lightness, chroma, separation of neighbours under normal vision and the
#' common colour-vision deficiencies, and contrast on white. Once they run
#' out, a scenario gets grey rather than a recycled colour. The reader can
#' change any colour on the Settings tab.
APP_PALETTE <- c("#0055A4", "#E0702E", "#A2337A", "#14936A", "#A88200", "#2B86C2", "#9A4A1E", "#6A55C2")
APP_PALETTE_EXTRA <- "#6F7680"

#' The scenarios of a bundle, with default colours.
#'
#' @param bundle The bundle.
#' @return data.frame id, label (short label), description, type, base (the
#'   id of the tilted forecast a scenario is built on, NA if none), base_label
#'   (its short label), colour,
#'   default (TRUE if listed in the header when the site opens); default
#'   scenarios first, in APP_DEFAULT_SCENARIOS order, then the others in
#'   bundle order.
scenario_catalog <- function(bundle) {
  o <- bundle$objects
  if (!nrow(o)) stop("the bundle has no forecasts")
  label <- if (is.null(o$short_label)) o$label else o$short_label
  first <- intersect(APP_DEFAULT_SCENARIOS, o$id)
  ord <- c(match(first, o$id), setdiff(seq_len(nrow(o)), match(first, o$id)))
  ids <- o$id[ord]
  # Fixed colours first; the others take the free palette colours in order.
  colour <- unname(APP_SCENARIO_COLOURS[ids])
  free <- setdiff(APP_PALETTE, APP_SCENARIO_COLOURS)
  k <- seq_len(sum(is.na(colour)))
  colour[is.na(colour)] <- ifelse(k <= length(free), free[pmin(k, length(free))], APP_PALETTE_EXTRA)
  base <- if (is.null(o$base)) rep(NA_character_, nrow(o)) else o$base
  base_label <- label[match(base, o$id)]
  data.frame(id = ids, label = label[ord], description = o$description[ord], type = o$type[ord],
             base = base[ord], base_label = base_label[ord], colour = colour, default = ids %in% first,
             stringsAsFactors = FALSE)
}

#' Short labels and level units of every variable in the bundle.
#'
#' `unit` is the unit of the series when it is shown as a level: natural
#' units for log series (x = exp(level / 100)), the published units for
#' rates and indices. Growth rates get their own units (unit_label()).
APP_VARIABLES <- data.frame(
  variable = c("rgdp", "rpce", "rresinv", "rnresinv", "rexp", "rimp", "rgov", "rfedgov",
               "gdpdefl", "pce", "corepce", "cpi", "corecpi", "comphr", "payems", "urate",
               "indpro", "caputil", "houst", "rdpi", "umich", "gs2", "gs5", "gs10", "aaa", "baa",
               "usd", "sp500", "wti", "gsci_nonenergy", "rvol", "term_10y_2y", "baa_aaa", "baa_10y"),
  label = c("Real GDP", "Real consumption", "Real residential investment", "Real business investment",
            "Real exports", "Real imports", "Real government spending", "Real federal government spending",
            "GDP deflator", "PCE price index", "Core PCE price index", "CPI", "Core CPI",
            "Compensation per hour", "Nonfarm payrolls", "Unemployment rate", "Industrial production",
            "Capacity utilization, manufacturing", "Housing starts", "Real disposable income",
            "Consumer sentiment (Michigan)", "2y Treasury yield", "5y Treasury yield", "10y Treasury yield",
            "Aaa corporate yield", "Baa corporate yield", "Dollar index (advanced economies)", "S&P 500",
            "WTI crude oil", "Non-energy commodities (proxy)", "S&P 500 realized volatility",
            "10y − 2y Treasury spread", "Baa − Aaa spread", "Baa − 10y spread"),
  unit = c("$bn, chained", "$bn, chained", "$bn, chained", "$bn, chained", "$bn, chained", "$bn, chained",
           "$bn, chained", "index", "index", "index", "index", "index", "index", "index", "thousands", "%",
           "index", "%", "thousands, annual rate", "$bn, chained", "index", "%", "%", "%", "%", "%",
           "index", "index", "$ per barrel", "index", "%, annualized", "pp", "pp", "pp"),
  stringsAsFactors = FALSE
)

#' Source and construction of every variable, for the About tab
#' (docs/data.md, sections 2 and 3). Monthly series are quarterly averages
#' (FRED-QD's, or the project's for series taken from the FRED API).
APP_SOURCES <- data.frame(
  variable = APP_VARIABLES$variable,
  source = c(
    "BEA, via FRED-QD (GDPC1)", "BEA, via FRED-QD (PCECC96)", "BEA, via FRED-QD (PRFIx)", "BEA, via FRED-QD (PNFIx)",
    "BEA, via FRED-QD (EXPGSC1)", "BEA, via FRED-QD (IMPGSC1)", "BEA, via FRED-QD (GCEC1)",
    "BEA, via FRED (B823RA3Q086SBEA)", "BEA, via FRED-QD (GDPCTPI)", "BEA, via FRED-QD (PCECTPI)",
    "BEA, via FRED-QD (PCEPILFE)", "BLS, via FRED-QD (CPIAUCSL)", "BLS, via FRED-QD (CPILFESL)",
    "BLS, via FRED (HCOMPBS)", "BLS, via FRED-QD (PAYEMS)", "BLS, via FRED-QD (UNRATE)",
    "Federal Reserve Board, via FRED-QD (INDPRO)", "Federal Reserve Board, via FRED-QD (CUMFNS)",
    "Census Bureau, via FRED-QD (HOUST)", "BEA, via FRED-QD (DPIC96)",
    "University of Michigan Surveys of Consumers, via FRED-QD (UMCSENTx)",
    "Federal Reserve Board, via FRED (GS2)", "Federal Reserve Board, via FRED-QD (GS5)",
    "Federal Reserve Board, via FRED-QD (GS10)", "Moody's, via FRED-QD (AAA)", "Moody's, via FRED-QD (BAA)",
    "Federal Reserve Board, via FRED-QD (TWEXAFEGSMTHx)", "S&P Dow Jones Indices, via FRED-QD (S&P 500)",
    "Federal Reserve Bank of St. Louis, FRED (WTISPLC)", "World Bank Commodity Price Data (Pink Sheet)",
    "Computed from Yahoo Finance daily S&P 500 closes (^GSPC)",
    "Computed: 10y minus 2y yield", "Computed: Baa minus Aaa yield", "Computed: Baa minus 10y yield"),
  note = c(
    "", "", "", "", "", "", "",
    "Chain-type quantity index (2017 = 100), used in place of chained dollars, which FRED does not publish for the whole sample",
    "", "", "", "", "", "", "", "", "", "", "", "", "", "Before 1976Q3: Gürkaynak–Sack–Wright 2-year par yield (Federal Reserve)",
    "", "", "", "",
    "Advanced foreign economies index, extended back by FRED-QD with the major-currencies index",
    "", "Nominal, dollars per barrel", "Non-energy index, nominal; a proxy for the S&P GSCI non-energy index used in the paper",
    "Root mean square of daily log returns, annualized, percent", "", "", ""),
  stringsAsFactors = FALSE
)

#' Variables whose observed history is left out of tables and downloads
#' (the publishers' terms restrict redistribution: Moody's, S&P Dow Jones
#' Indices, Yahoo Finance). Their charts show the history, with the source
#' credited; their forecasts are the model's and are included.
APP_NO_HISTORY_DOWNLOAD <- c("aaa", "baa", "sp500", "rvol", "baa_aaa", "baa_10y")

#' Topic tabs: id, title, variables (in display order).
APP_TOPICS <- list(
  list(id = "output", title = "Output and demand",
       variables = c("rgdp", "rpce", "rexp", "rimp", "rgov", "rfedgov", "rdpi")),
  list(id = "investment", title = "Investment and housing", variables = c("rresinv", "rnresinv", "houst")),
  list(id = "labor", title = "Labor and activity", variables = c("payems", "urate", "comphr", "indpro", "caputil")),
  list(id = "prices", title = "Prices", variables = c("gdpdefl", "pce", "corepce", "cpi", "corecpi")),
  list(id = "rates", title = "Rates and credit",
       variables = c("gs2", "gs5", "gs10", "aaa", "baa", "term_10y_2y", "baa_aaa", "baa_10y")),
  list(id = "markets", title = "Markets, commodities, sentiment",
       variables = c("sp500", "rvol", "usd", "wti", "gsci_nonenergy", "umich"))
)

#' Default charts of the overview, and the most it can hold.
APP_OVERVIEW_DEFAULT <- c("rgdp", "pce", "corepce", "urate", "gs2", "gs10")
APP_OVERVIEW_MAX <- 9

#' History windows (quarters of history shown before the forecast).
APP_HISTORY <- c("5y" = 20, "10y" = 40, "20y" = 80, "All" = Inf)

#' Unit label of a chart.
#'
#' @param variable Variable name.
#' @param transform Transform shown ("level", "natural", "yoy", "qoq_ann").
#' @param difference TRUE when the chart shows differences from the baseline.
#' @return Character string, e.g. "% y/y" or "pp".
unit_label <- function(variable, transform, difference = FALSE) {
  unit <- APP_VARIABLES$unit[match(variable, APP_VARIABLES$variable)]
  if (transform %in% c("yoy", "qoq_ann")) {
    growth <- if (transform == "yoy") "% y/y" else "% q/q ann."
    return(if (difference) paste0("pp (", growth, ")") else growth)
  }
  if (difference && unit %in% c("%", "%, annualized", "pp")) return("pp")
  unit
}

#' Short label of a variable.
var_label <- function(variable) APP_VARIABLES$label[match(variable, APP_VARIABLES$variable)]
