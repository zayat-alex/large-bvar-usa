# bundle.R: load the export bundle, the only data the app reads.
#
# The bundle (bundle/bundle.rds) is written by the project's export stage:
# quantile tables of every forecast, the observed history and the impulse
# responses, already in display units. The app runs no model code; if a
# table or column it needs is missing, it stops with a message saying which.


#' Tables and columns the app relies on.
BUNDLE_SCHEMA <- list(
  variables = c("variable", "label", "kind", "model_transform", "display"),
  objects   = c("id", "label", "description", "type"),
  history   = c("variable", "transform", "quarter", "date", "value"),
  forecasts = c("object", "variable", "transform", "quarter", "h", "date", "q05", "q16", "q50", "q84", "q95"),
  irf       = c("identification", "scheme", "shock", "variable", "transform", "horizon", "q05", "q16", "q50",
                "q84", "q95"),
  meta      = c("run_id", "data", "forecast")
)

#' Read and check the bundle.
#'
#' @param path Path of bundle.rds.
#' @return The bundle (a list of tables plus meta).
load_bundle <- function(path) {
  if (!file.exists(path)) {
    stop("bundle not found at ", path, ": run scripts/publish_bundle.R from the project root", call. = FALSE)
  }
  b <- readRDS(path)
  check_bundle(b)
  b
}

#' Check that a bundle has every table and column the app uses.
#'
#' @param b A bundle.
#' @return TRUE invisibly; stops naming the first missing table or columns.
check_bundle <- function(b) {
  for (tab in names(BUNDLE_SCHEMA)) {
    if (is.null(b[[tab]])) stop("the bundle has no table '", tab, "'", call. = FALSE)
    miss <- setdiff(BUNDLE_SCHEMA[[tab]], names(b[[tab]]))
    if (length(miss)) stop("bundle table '", tab, "' lacks: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}
