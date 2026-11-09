# covid.R: the COVID treatment stage (pipeline stage 2), between the data and
# the estimation.
#
# Why it exists: the 2020 pandemic quarters are far outside anything else in
# the 1973-2019 data (real GDP fell by 7.9 percent in 2020:Q2 and rose by 7.8
# percent in 2020:Q3, quarter on quarter, QD vintage 2026m09). A VAR with constant parameters and a
# constant innovation covariance Sigma treats them as ordinary, very large
# shocks, which inflates the estimated Sigma and moves the hyperparameters.
# How to handle them is a modelling choice the paper does not face (its sample
# ends in 2019). This stage isolates that choice: it decides which quarters
# the model is estimated on and which history the forecasts start from, and
# nothing downstream needs to know which method was used.
#
# Methods (config/run.yml, covid$method):
#   "none"            : estimate on the full sample; forecast from the latest data.
#   "precovid_params" : estimate on the sample up to covid$cutoff (default
#                       2019Q4), so the parameters reflect pre-2020 dynamics;
#                       forecast from the latest data.
#   "exclude_window"  : estimate on the full sample but leave out the
#                       equations of the quarters in covid$window (default
#                       2020Q1-2021Q4), so the parameters reflect the data
#                       before and after the pandemic; forecast from the
#                       latest data. As in Schorfheide and Song (2024), who
#                       exclude the pandemic observations from estimation.
#
# Input:  model_data, a data.frame with a Date column `quarter` (first day of
#         each quarter) and the 31 model variables in model units (100*ln or
#         raw; docs/data.md, section 1), one row per quarter, no missing values;
#         and the run configuration (the list read from config/run.yml).
# Output: a list with exactly four elements:
#   est_data    : data.frame, same columns as model_data; the quarters used to
#                 estimate the model.
#   origin_data : data.frame, same columns as model_data; the actual history
#                 from which forecasts start (the engine uses its last `lags`
#                 rows; the export uses it for the history shown with forecasts).
#   est_weights : NULL, or a numeric vector with one entry per row of est_data:
#                 1 for the quarters whose equation is used in the
#                 estimation, 0 for those left out (method exclude_window).
#                 The left-out quarters stay in est_data, since they still
#                 serve as lags of later quarters.
#   notes       : data.frame with one row per quarter that a method dropped or
#                 changed (columns quarter, quarter_label, method, action, reason).
#                 Zero rows when nothing was changed.


# ---- Main entry point ----

#' Apply the configured COVID treatment to the model data.
#'
#' @param model_data data.frame: `quarter` (Date) plus the 31 model variables.
#' @param cfg Run configuration (list from config/run.yml). Uses
#'   cfg$covid$method, cfg$covid$cutoff and cfg$model$lags.
#' @return List with est_data, origin_data, est_weights, notes (see file header).
#'   Stops with an error if the configuration or the result is invalid.
covid_treatment <- function(model_data, cfg) {
  # Step 1: read the settings this stage uses.
  method <- cfg$covid$method
  lags   <- cfg$model$lags

  # Step 2: check the settings against the data before doing anything.
  check_covid_config(cfg$covid, model_data, lags)

  # Step 3: apply the chosen method. Each method is its own function so that
  # its rule can be read in one place.
  if (method == "none") {
    out <- covid_none(model_data)
  } else if (method == "precovid_params") {
    cutoff <- parse_quarter(cfg$covid$cutoff)
    out <- covid_precovid_params(model_data, cutoff)
  } else if (method == "exclude_window") {
    window <- parse_quarter(unlist(cfg$covid$window))
    out <- covid_exclude_window(model_data, window)
  }

  # Step 4: check that the result has the agreed shape before passing it on.
  check_covid_output(out, model_data, lags)

  out
}


# ---- Methods ----

#' Method "none": no treatment.
#'
#' The model is estimated on every quarter, including 2020-2021, and forecasts
#' start from the same data.
#'
#' @param model_data data.frame: `quarter` plus the model variables.
#' @return The four-element list described in the file header.
covid_none <- function(model_data) {
  # Step 1: the estimation data are the full sample.
  est_data <- model_data

  # Step 2: the forecast history is the full sample as well.
  origin_data <- model_data

  # Step 3: nothing was dropped or changed, so the notes table is empty.
  notes <- make_covid_notes(quarters = as.Date(character(0)), method = "none",
                            action = character(0), reason = character(0))

  list(est_data = est_data, origin_data = origin_data,
       est_weights = NULL, notes = notes)
}

#' Method "precovid_params": estimate on data up to a cutoff quarter.
#'
#' The parameters (coefficients, Sigma, hyperparameters) are estimated only on
#' quarters up to and including `cutoff`, so the 2020-2021 observations do not
#' affect them. Forecasts still start from the latest actual data, so the
#' pre-2020 dynamics are applied to the current state of the economy. The
#' price of this choice is that information after the cutoff is not used to
#' estimate the parameters at all.
#'
#' @param model_data data.frame: `quarter` plus the model variables.
#' @param cutoff Date: first day of the last quarter used for estimation.
#' @return The four-element list described in the file header.
covid_precovid_params <- function(model_data, cutoff) {
  # Step 1: mark the quarters kept for estimation (on or before the cutoff).
  keep <- model_data$quarter <= cutoff

  # Step 2: the estimation data are those quarters only.
  est_data <- model_data[keep, , drop = FALSE]
  rownames(est_data) <- NULL

  # Step 3: the forecast history is the full actual sample.
  origin_data <- model_data

  # Step 4: record every quarter left out of the estimation data.
  dropped <- model_data$quarter[!keep]
  notes <- make_covid_notes(
    quarters = dropped,
    method   = "precovid_params",
    action   = rep("dropped from est_data", length(dropped)),
    reason   = rep(paste0("after the estimation cutoff ", format_quarter(cutoff)),
                   length(dropped))
  )

  list(est_data = est_data, origin_data = origin_data,
       est_weights = NULL, notes = notes)
}

#' Method "exclude_window": leave the pandemic quarters' equations out.
#'
#' The equation of quarter t relates y_t to the p quarters before it.
#' The equations of the quarters in the window (by default 2020Q1-2021Q4) are
#' left out of the estimation; all other equations are used, including those
#' after the window. A left-out quarter's data still enter later equations as
#' lags: with p = 5 and the default window, the first equation used after the
#' window (2022Q1) has the lags 2020Q4-2021Q4, so the most extreme quarters
#' (2020Q2-Q3) appear in no equation used. Excluding observations from the
#' estimation is the treatment of Schorfheide and Song (2024); Lenza and
#' Primiceri (2022) note that it suits the estimation of the parameters, while
#' it understates uncertainty in forecasts made during the pandemic.
#'
#' Forecasts start from the latest actual data.
#'
#' @param model_data data.frame: `quarter` plus the model variables.
#' @param window Date vector of length 2: first and last quarter whose
#'   equations are left out.
#' @return The four-element list described in the file header.
covid_exclude_window <- function(model_data, window) {
  # Step 1: mark the quarters whose equation is left out.
  out_of_sample <- model_data$quarter >= window[1] & model_data$quarter <= window[2]

  # Step 2: the estimation data are the full sample (the left-out quarters
  # are still needed as lags); the weights say which equations are used.
  est_data <- model_data
  rownames(est_data) <- NULL
  est_weights <- as.numeric(!out_of_sample)

  # Step 3: the forecast history is the full actual sample.
  origin_data <- model_data

  # Step 4: record every quarter whose equation is left out.
  excluded <- model_data$quarter[out_of_sample]
  notes <- make_covid_notes(
    quarters = excluded,
    method   = "exclude_window",
    action   = rep("equation excluded from the estimation (data kept as lags)", length(excluded)),
    reason   = rep(paste0("inside the window ", format_quarter(window[1]), "-", format_quarter(window[2])),
                   length(excluded))
  )

  list(est_data = est_data, origin_data = origin_data,
       est_weights = est_weights, notes = notes)
}

#' Which equations a COVID treatment uses in the estimation.
#'
#' @param covid Output of covid_treatment().
#' @return NULL when every equation is used (no weights), else a logical
#'   vector with one entry per row of est_data (TRUE = equation used).
covid_keep <- function(covid) {
  if (is.null(covid$est_weights)) return(NULL)
  covid$est_weights > 0
}


# ---- Helpers ----

#' Build the notes table that records what a method did to each quarter.
#'
#' @param quarters Date vector of the quarters concerned (may be empty).
#' @param method Name of the method (one string).
#' @param action,reason Character vectors, one entry per quarter.
#' @return data.frame with columns quarter (Date), quarter_label ("2020Q1"),
#'   method, action, reason; one row per quarter.
make_covid_notes <- function(quarters, method, action, reason) {
  data.frame(
    quarter       = as.Date(quarters),
    quarter_label = format_quarter(quarters),
    method        = rep(method, length(quarters)),
    action        = action,
    reason        = reason,
    stringsAsFactors = FALSE
  )
}


# ---- Checks ----

#' Check the COVID settings against the data.
#'
#' Fails loudly, with a message naming the problem, rather than letting a
#' wrong setting produce a silently different estimation sample.
#'
#' @param covid_cfg The `covid` section of the run configuration.
#' @param model_data data.frame: `quarter` plus the model variables.
#' @param lags Number of VAR lags (cfg$model$lags).
#' @return TRUE invisibly; stops on the first problem found.
check_covid_config <- function(covid_cfg, model_data, lags) {
  implemented <- c("none", "precovid_params", "exclude_window")
  planned     <- c("interpolate_window", "scaling")
  method      <- covid_cfg$method

  # Check 1: the method is one string and is implemented.
  if (is.null(method) || length(method) != 1) {
    stop("covid$method must be one string, one of: ", paste(implemented, collapse = ", "))
  }
  if (method %in% planned) {
    stop("covid$method '", method, "' is not implemented yet")
  }
  if (!method %in% implemented) {
    stop("unknown covid$method '", method, "'; expected one of: ",
         paste(implemented, collapse = ", "))
  }

  # Check 2: the number of lags is a positive whole number.
  if (is.null(lags) || length(lags) != 1 || lags < 1 || lags != round(lags)) {
    stop("model$lags must be a positive whole number")
  }

  # Check 3 (precovid_params only): the cutoff is a valid quarter inside the
  # sample, before its last quarter, and leaves more than `lags` quarters,
  # since the first `lags` quarters only serve as initial lags.
  if (method == "precovid_params") {
    if (is.null(covid_cfg$cutoff)) stop("covid$cutoff is required for method precovid_params")
    cutoff <- parse_quarter(covid_cfg$cutoff)
    first  <- min(model_data$quarter)
    last   <- max(model_data$quarter)
    if (cutoff < first) {
      stop("covid$cutoff ", format_quarter(cutoff), " is before the sample start ",
           format_quarter(first))
    }
    if (cutoff >= last) {
      stop("covid$cutoff ", format_quarter(cutoff), " is not before the last quarter ",
           format_quarter(last), "; nothing would be dropped (use method none)")
    }
    n_kept <- sum(model_data$quarter <= cutoff)
    if (n_kept <= lags) {
      stop("covid$cutoff ", format_quarter(cutoff), " leaves ", n_kept,
           " quarters, not more than the ", lags, " lags")
    }
  }

  # Check 4 (exclude_window only): the window is two quarters in order,
  # after the initial lags and inside the sample, and leaves equations to
  # estimate on both sides.
  if (method == "exclude_window") {
    w <- unlist(covid_cfg$window)
    if (is.null(w) || length(w) != 2) stop("covid$window must give two quarters (first and last)")
    w <- parse_quarter(w)
    if (w[1] > w[2]) stop("covid$window must start before it ends")
    first <- min(model_data$quarter)
    last  <- max(model_data$quarter)
    if (sum(model_data$quarter < w[1]) <= lags) {
      stop("covid$window starts within the first ", lags, " quarters (the initial lags)")
    }
    if (w[2] >= last) stop("covid$window must end before the last quarter ", format_quarter(last))
  }

  invisible(TRUE)
}

#' Check that a COVID treatment result has the agreed shape.
#'
#' These are the properties the estimation stage and the engine rely on, so
#' they are checked once here instead of being assumed downstream.
#'
#' @param out Result of a method function.
#' @param model_data The model_data the method received.
#' @param lags Number of VAR lags.
#' @return TRUE invisibly; stops on the first problem found.
check_covid_output <- function(out, model_data, lags) {
  # Check 1: exactly the four agreed elements, in order.
  expected <- c("est_data", "origin_data", "est_weights", "notes")
  if (!identical(names(out), expected)) {
    stop("COVID treatment output must have elements ", paste(expected, collapse = ", "))
  }

  # Check 2: est_data and origin_data have the columns of model_data, in the
  # same order (the order fixes which column is which variable in the VAR),
  # and no missing values.
  for (name in c("est_data", "origin_data")) {
    x <- out[[name]]
    if (!is.data.frame(x)) stop(name, " must be a data.frame")
    if (!identical(names(x), names(model_data))) {
      stop(name, " must have the columns of model_data in the same order")
    }
    if (anyNA(x)) stop(name, " has missing values")
  }

  # Check 3: est_data is an unbroken run of quarters starting at the sample
  # start (a VAR needs consecutive observations), with more rows than lags.
  est_q    <- out$est_data$quarter
  expected_q <- seq(min(model_data$quarter), by = "3 months", length.out = length(est_q))
  if (!identical(as.numeric(est_q), as.numeric(expected_q))) {
    stop("est_data must be consecutive quarters starting at ",
         format_quarter(min(model_data$quarter)))
  }
  if (nrow(out$est_data) <= lags) stop("est_data must have more than ", lags, " rows")

  # Check 4: origin_data ends at the last actual quarter and has at least
  # `lags` rows, so a forecast can start from it.
  if (max(out$origin_data$quarter) != max(model_data$quarter)) {
    stop("origin_data must end at the last quarter of model_data, ",
         format_quarter(max(model_data$quarter)))
  }
  if (nrow(out$origin_data) < lags) stop("origin_data must have at least ", lags, " rows")

  # Check 5: weights, if any, are 0 or 1, one per row of est_data, and leave
  # more equations than lags (a method that down-weights quarters gradually
  # would need a different check).
  if (!is.null(out$est_weights)) {
    w <- out$est_weights
    if (!is.numeric(w) || length(w) != nrow(out$est_data) || anyNA(w) || !all(w %in% c(0, 1))) {
      stop("est_weights must be 0 or 1, one per row of est_data")
    }
    if (sum(w[-seq_len(lags)]) <= lags) stop("est_weights leave too few equations to estimate")
  }

  # Check 6: the notes table has the agreed columns.
  notes_cols <- c("quarter", "quarter_label", "method", "action", "reason")
  if (!is.data.frame(out$notes) || !identical(names(out$notes), notes_cols)) {
    stop("notes must be a data.frame with columns ", paste(notes_cols, collapse = ", "))
  }

  invisible(TRUE)
}
