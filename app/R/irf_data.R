# irf_data.R: impulse responses, taken from the bundle.
#
# An impulse response shows how the model's forecast of each variable moves
# after a one-time shock to one variable, compared with no shock (Crump et
# al. 2021, Sections 2.3 and 3). The bundle holds, per identification,
# quantiles over posterior draws of the response to a unit shock: 1 log
# point (about 1%) for a log series, 1 percentage point for a rate. Two
# kinds of identification:
#   generalized: the shocked variable is ordered first in a Cholesky
#     factorization, so the other variables move on impact as they
#     typically do together with it (a forecast revision, not a structural
#     shock);
#   recursive: the variables are ordered in blocks (config/run.yml, irf),
#     and those ordered before the shocked variable cannot move on impact.
# The shocks listed on the site are those in the bundle (irf_shocks()).
# Responses are linear in the shock, so a shock of size s scales every
# quantile by s (and swaps the lower and upper bands when s < 0). The bands
# cover parameter uncertainty only.


#' The shock selected when the site opens (if the bundle has it).
APP_IRF_DEFAULT_SHOCK <- "recursive_gs2"

#' Responses shown by default: chosen sets for some shocks; any other shock
#' (for example one added later to config/run.yml) shows its own variable
#' followed by the general set, up to six charts.
APP_IRF_DEFAULT_RESPONSES <- list(
  recursive_gs2 = c("rgdp", "corepce", "urate", "gs2", "gs10", "sp500"),
  girf_gs2 = c("rgdp", "corepce", "urate", "gs2", "gs10", "sp500"),
  girf_rgdp = c("rgdp", "payems", "urate", "corepce", "gs2", "sp500"),
  girf_wti = c("wti", "pce", "corepce", "cpi", "rgdp", "gs2"),
  girf_sp500 = c("sp500", "rgdp", "rnresinv", "urate", "gs2", "baa_10y"),
  girf_baa = c("baa", "baa_10y", "rgdp", "urate", "sp500", "gs2")
)
APP_IRF_GENERAL_RESPONSES <- c("rgdp", "corepce", "urate", "gs2", "gs10", "sp500")

APP_IRF_HORIZONS <- c(8, 12, 16, 20)

#' The shocks of a bundle, named for readers.
#'
#' Every identification in the bundle's irf table is listed: the default
#' first, then the recursive ones, then the generalized ones, each in
#' variable order. The name is the shocked variable's label; the label adds
#' the kind of identification.
#'
#' @param bundle The bundle.
#' @return data.frame id, scheme, shock (variable), name, kind, label.
irf_shocks <- function(bundle) {
  x <- unique(bundle$irf[, c("identification", "scheme", "shock")])
  names(x)[1] <- "id"
  scheme_name <- c(recursive = "recursive", girf = "generalized")
  name <- var_label(x$shock)
  name[is.na(name)] <- x$shock[is.na(name)]
  kind <- ifelse(x$scheme %in% names(scheme_name), scheme_name[x$scheme], x$scheme)
  x$name <- name
  x$kind <- kind
  x$label <- paste0(name, " (", kind, ")")
  ord <- order(x$id != APP_IRF_DEFAULT_SHOCK, x$scheme != "recursive", match(x$shock, bundle$variables$variable))
  out <- x[ord, ]
  rownames(out) <- NULL
  out
}

#' Default responses of a shock.
#'
#' @param bundle The bundle.
#' @param shock_id Identification id.
#' @return Variable names (in the bundle), at most six.
irf_default_responses <- function(bundle, shock_id) {
  set <- APP_IRF_DEFAULT_RESPONSES[[shock_id]]
  if (is.null(set)) set <- c(sub("^[a-z]+_", "", shock_id), APP_IRF_GENERAL_RESPONSES)
  utils::head(intersect(unique(set), bundle$variables$variable), 6)
}

#' Unit of a response.
#'
#' @param vars The bundle's variables table.
#' @param variable Variable name.
#' @param transform "level" or "yoy" (4-quarter change, log series only).
#' @return Character, e.g. "%", "pp", "pp, 4-quarter change".
irf_unit <- function(vars, variable, transform) {
  if (vars$model_transform[match(variable, vars$variable)] == "log") {
    return(if (transform == "yoy") "pp, 4-quarter change" else "%")
  }
  unit <- APP_VARIABLES$unit[match(variable, APP_VARIABLES$variable)]
  if (isTRUE(unit == "index")) "index points" else "pp"
}

#' Size unit of a shock, for the shock-size hint.
#'
#' @param vars The bundle's variables table.
#' @param shock_id Identification id (e.g. "girf_wti").
#' @return "log points (about %)" or "pp".
irf_shock_unit <- function(vars, shock_id) {
  shock <- sub("^[a-z]+_", "", shock_id)
  if (vars$model_transform[match(shock, vars$variable)] == "log") "log points (about %)" else "pp"
}

#' One response, ready to plot.
#'
#' A variable without a 4-quarter change (a rate) is shown in levels.
#'
#' @param bundle The bundle.
#' @param shock_id Identification id.
#' @param variable Response variable.
#' @param transform "level" or "yoy".
#' @param horizon Last horizon shown (quarters after the shock).
#' @param size Shock size (multiplies the unit response).
#' @return data.frame horizon, q05, q16, q50, q84, q95, with attribute
#'   `transform` (the one used).
irf_series <- function(bundle, shock_id, variable, transform, horizon, size = 1) {
  x <- bundle$irf[bundle$irf$identification == shock_id & bundle$irf$variable == variable, ]
  if (!nrow(x)) stop("no response of ", variable, " to ", shock_id)
  if (!transform %in% x$transform) transform <- "level"
  x <- x[x$transform == transform & x$horizon <= horizon, c("horizon", "q05", "q16", "q50", "q84", "q95")]
  x <- x[order(x$horizon), ]
  q <- c("q05", "q16", "q50", "q84", "q95")
  x[q] <- x[q] * size
  # A negative shock turns the bands upside down: lower and upper swap.
  if (size < 0) x[q] <- x[rev(q)]
  rownames(x) <- NULL
  attr(x, "transform") <- transform
  x
}

#' The numbers of the IRF charts, as one table.
#'
#' @param bundle,shock_id,transform,horizon,size As in irf_series().
#' @param variables Response variables.
#' @return data.frame: Variable, Unit, Quarters after shock, Median and the
#'   68% and 90% bands.
irf_data <- function(bundle, shock_id, variables, transform, horizon, size = 1) {
  rows <- lapply(variables, function(v) {
    x <- irf_series(bundle, shock_id, v, transform, horizon, size)
    data.frame(Variable = var_label(v), Unit = irf_unit(bundle$variables, v, attr(x, "transform")),
               `Quarters after shock` = x$horizon, Median = x$q50, `68% low` = x$q16, `68% high` = x$q84,
               `90% low` = x$q05, `90% high` = x$q95, check.names = FALSE)
  })
  if (!length(rows)) return(data.frame(Variable = character(0)))
  do.call(rbind, rows)
}
