# ui.R: the page layout.
#
# A full-height side panel (title and tabs, Forecasts grouping the topic
# tabs) and, beside it, a header with the settings, which stay the same
# across the forecast tabs: which scenarios to show, with or without
# bands, which one is the baseline, forecasts or differences from the
# baseline, the growth view and the history window.
# On the impulse-response tab the header holds the response settings
# instead. Below the header, the current tab's charts.


#' A short line in a scenario's colour, for the scenario list.
swatch <- function(colour) {
  htmltools::HTML(sprintf(paste0('<svg width="18" height="10" aria-hidden="true" class="swatch">',
                                 '<line x1="0" y1="5" x2="18" y2="5" stroke="%s" stroke-width="2"/></svg>'), colour))
}

#' Text read by screen readers only.
sr_only <- function(text) htmltools::span(class = "visually-hidden", text)

#' The scenario list: Show, Bands and Baseline columns, one row per scenario.
#'
#' Three standard inputs side by side (two checkbox groups and a radio
#' group), aligned row by row by the stylesheet. Built by the server, since
#' the reader chooses on the Settings tab which scenarios it lists and in
#' which colours.
#'
#' Each scenario's name carries its description as a tooltip.
#'
#' @param scenarios The listed scenarios (rows of scenario_catalog(), with
#'   the current names and colours).
#' @param show,bands,baseline Ids selected in each column.
#' @param show_bands FALSE to leave out the Bands column (difference view).
scenario_controls <- function(scenarios, show, bands, baseline, show_bands = TRUE) {
  ids <- scenarios$id
  if (!length(ids)) {
    return(htmltools::tags$fieldset(class = "setting", htmltools::tags$legend("Scenarios"),
                                    htmltools::p(class = "empty-note", "No scenario listed: choose some on the Settings tab.")))
  }
  htmltools::tags$fieldset(
    class = "setting",
    htmltools::tags$legend("Scenarios"),
    htmltools::div(
      class = "scenario-grid",
      shiny::checkboxGroupInput("show", "Show", selected = show, choiceValues = ids,
                                choiceNames = lapply(seq_along(ids), function(k) {
                                  desc <- if (is.null(scenarios$description)) NULL else scenarios$description[k]
                                  htmltools::span(class = "scenario-name", title = desc,
                                                  swatch(scenarios$colour[k]), scenarios$label[k])
                                })),
      if (show_bands) shiny::checkboxGroupInput("bands", "Bands", selected = bands, choiceValues = ids,
                                choiceNames = lapply(scenarios$label, function(l) sr_only(paste("Bands for", l)))),
      shiny::radioButtons("baseline", "Baseline", selected = baseline, choiceValues = ids,
                          choiceNames = lapply(scenarios$label, function(l) sr_only(paste(l, "as baseline"))))
    )
  )
}

#' The settings header: the forecast settings, or on the impulse-response
#' tab the response settings. Both sets keep their values when hidden.
settings_header <- function(bundle) {
  shocks <- irf_shocks(bundle)
  htmltools::tags$header(
    class = "settings",
    # The rows sit inside the conditional panels: a panel is shown with
    # display: block, which would undo a flex layout set on the panel itself.
    shiny::conditionalPanel(
      "input.page !== 'irf'",
      htmltools::div(class = "settings-row",
        shiny::uiOutput("scenario_controls", inline = FALSE),
        htmltools::div(class = "setting",
          shiny::radioButtons("mode", "Show", choiceValues = c("forecast", "difference"),
                              choiceNames = c("Forecasts", "Difference from baseline (medians)"))),
        htmltools::div(class = "setting",
          shiny::radioButtons("growth", "Growth", choiceValues = c("yoy", "qoq_ann", "natural"),
                              choiceNames = c("Year over year", "Quarterly, annualized", "Level"))),
        htmltools::div(class = "setting setting-history",
          shiny::radioButtons("history", "History", choices = names(APP_HISTORY), inline = TRUE)))
    ),
    shiny::conditionalPanel(
      "input.page === 'irf'",
      htmltools::div(class = "settings-row",
        # A drop-down grouped by identification: one shock per variable is too
        # long a list for radio buttons.
        htmltools::div(class = "setting setting-shock",
          shiny::selectInput("irf_shock", "Shock", selected = shocks$id[1], selectize = TRUE, width = "300px",
                             choices = lapply(split(shocks, factor(tools::toTitleCase(shocks$kind),
                                                                   levels = unique(tools::toTitleCase(shocks$kind)))),
                                              function(g) stats::setNames(g$id, g$name)))),
        htmltools::div(class = "setting",
          shiny::radioButtons("irf_units", "Response of log series", choiceValues = c("level", "yoy"),
                              choiceNames = c("Level (%)", "4-quarter change (pp)"))),
        htmltools::div(class = "setting setting-history",
          shiny::radioButtons("irf_horizon", "Quarters shown", choices = APP_IRF_HORIZONS,
                              selected = max(APP_IRF_HORIZONS), inline = TRUE)),
        htmltools::div(class = "setting setting-size",
          shiny::numericInput("irf_size", "Shock size", value = 1, step = 0.25, width = "90px"),
          shiny::uiOutput("irf_size_unit", class = "size-unit")))
    )
  )
}

#' One tab in the side panel: a radio input of the `page` group.
#'
#' The tabs are written by hand rather than with radioButtons() so they can
#' be grouped; Shiny reads every radio named "page" inside the group's
#' container, wherever it sits.
#'
#' @param value,label The tab's id and name.
#' @param selected TRUE for the tab shown first.
#' @param class "sub" for a tab inside a group, "top" for a top-level tab.
page_tab <- function(value, label, selected = FALSE, class = "top") {
  htmltools::div(class = paste("radio", class),
    htmltools::tags$label(
      htmltools::tags$input(type = "radio", name = "page", value = value, checked = if (selected) NA),
      htmltools::span(label)))
}

#' The side panel: title and tabs. Forecasts groups the overview and the
#' topic tabs; Impulse responses sits at the same level; Settings and
#' About below a rule.
side_panel <- function() {
  topics <- lapply(APP_TOPICS, function(t) page_tab(t$id, t$title, class = "sub"))
  htmltools::tags$aside(
    class = "side",
    htmltools::h1("US Economy - BVAR Forecasts"),
    htmltools::tags$nav(`aria-label` = "Tabs",
      htmltools::div(id = "page", class = "form-group shiny-input-radiogroup shiny-input-container",
                     role = "radiogroup", `aria-label` = "Tabs",
        htmltools::div(class = "shiny-options-group",
          htmltools::div(class = "nav-group", "Forecasts"),
          page_tab("overview", "Overview", selected = TRUE, class = "sub"),
          topics,
          page_tab("irf", "Impulse responses"),
          htmltools::div(class = "nav-rule"),
          page_tab("settings", "Settings"),
          page_tab("about", "About the model and data"))))
  )
}

#' The impulse-response tab: response selector, charts, a note and data.
irf_panel <- function(bundle) {
  vars <- intersect(APP_VARIABLES$variable, bundle$variables$variable)
  bslib::nav_panel_hidden(
    "irf",
    htmltools::div(
      class = "chart-picker",
      shiny::selectizeInput("irf_vars", "Responses", choices = stats::setNames(vars, var_label(vars)),
                            selected = irf_default_responses(bundle, irf_shocks(bundle)$id[1]), multiple = TRUE,
                            options = list(maxItems = APP_OVERVIEW_MAX, plugins = list("remove_button")),
                            width = "100%"),
      shiny::actionLink("irf_reset", "Reset to the shock's default responses")
    ),
    shiny::uiOutput("grid_irf"),
    htmltools::p(class = "irf-note",
      "Each chart shows how the model's forecast of a variable moves after a one-time shock, compared with ",
      "no shock: the median and 68% and 90% credible bands over parameter uncertainty. A generalized shock ",
      "orders the shocked variable first, so the other variables move on impact as they typically do with it; ",
      "it describes co-movement, not a structural cause. A recursive shock orders the variables in blocks: ",
      "those ordered before the shocked variable cannot respond on impact, those after it can. Responses ",
      "therefore depend on the identification. A unit shock is 1 log point (about 1%) for a log series and ",
      "1 percentage point for a rate; responses are linear, so the shock size scales them exactly."),
    data_section("irf")
  )
}

#' The overview tab: chart selector and charts.
overview_panel <- function(bundle) {
  vars <- intersect(APP_VARIABLES$variable, bundle$variables$variable)
  bslib::nav_panel_hidden(
    "overview",
    htmltools::div(
      class = "chart-picker",
      shiny::selectizeInput("overview_vars", "Charts", choices = stats::setNames(vars, var_label(vars)),
                            selected = APP_OVERVIEW_DEFAULT, multiple = TRUE,
                            options = list(maxItems = APP_OVERVIEW_MAX, plugins = list("remove_button")),
                            width = "100%"),
      shiny::actionLink("overview_reset", "Reset to default")
    ),
    shiny::uiOutput("grid_overview"),
    data_section("overview")
  )
}

#' One row of the scenario table on the Settings tab: list it in the
#' header, its colour, its name and what it assumes.
#'
#' The colour box is the browser's own colour picker; its value reaches the
#' server as input colour_<id> when the reader changes it. The name box
#' renames the scenario everywhere on the site (input name_<id>).
#'
#' @param row One row of scenario_catalog().
#' @param listed TRUE if the scenario is listed in the header.
#' @param colour,name Its current colour and name.
scenario_settings_row <- function(row, listed, colour, name) {
  colour_id <- paste0("colour_", row$id)
  type <- switch(row$type, baseline = "unconditional forecast", tilted = "tilted to targets (soft conditions)",
                 hard = "conditional forecast (fixed paths)", row$type)
  htmltools::tags$tr(
    htmltools::tags$td(class = "col-narrow",
      shiny::checkboxInput(paste0("list_", row$id), sr_only(paste("List", row$label)), listed)),
    htmltools::tags$td(class = "col-narrow",
      htmltools::tags$label(class = "visually-hidden", `for` = colour_id, paste("Colour of", row$label)),
      htmltools::tags$input(type = "color", id = colour_id, value = colour, class = "colour-input",
                            oninput = "Shiny.setInputValue(this.id, this.value)")),
    htmltools::tags$td(class = "col-name",
      shiny::textInput(paste0("name_", row$id), sr_only(paste("Name of", row$label)), name),
      htmltools::div(class = "scenario-type", type)),
    htmltools::tags$td(class = "scenario-description", row$description)
  )
}

#' The Settings tab. For now one section: the scenarios the header lists,
#' their names and colours.
settings_panel <- function() {
  bslib::nav_panel_hidden(
    "settings",
    htmltools::div(class = "settings-page",
      htmltools::h2("Scenarios"),
      htmltools::p("Choose the scenarios listed in the header, their names and their colours. ",
                   "The choices last for this visit."),
      shiny::uiOutput("scenario_table"),
      shiny::actionLink("reset_scenarios", "Reset scenarios to default")
    )
  )
}

#' About tab: what the model is and where the data come from.
about_panel <- function(bundle) {
  m <- bundle$meta
  v <- bundle$variables
  model_units <- ifelse(v$model_transform == "log", "100 × ln(level)", "as published")
  shown_as <- c(growth = "growth rate or level (Growth setting)", natural = "level, in natural units",
                level = "level")[v$display]
  src <- APP_SOURCES[match(v$variable, APP_SOURCES$variable), ]
  data_table <- htmltools::tags$table(
    class = "about-table",
    htmltools::tags$thead(htmltools::tags$tr(lapply(c("Variable", "Description", "Source", "In the model",
                                                      "Shown as", "Notes"), htmltools::tags$th))),
    htmltools::tags$tbody(lapply(seq_len(nrow(v)), function(k) {
      htmltools::tags$tr(htmltools::tags$td(var_label(v$variable[k])), htmltools::tags$td(v$label[k]),
                         htmltools::tags$td(src$source[k]), htmltools::tags$td(model_units[k]),
                         htmltools::tags$td(shown_as[[k]]), htmltools::tags$td(src$note[k]))
    }))
  )
  # Counts, sample, window and horizon come from the bundle, so the text
  # follows the published run.
  n_draws <- m$forecast$n_draws %||% NA
  n_paths <- m$forecast$paths_per_draw %||% NA
  n_series <- sum(v$kind == "series")
  window <- unlist(m$model$covid_window)
  window_text <- if (length(window) == 2) {
    sprintf(" The equations of %s–%s are left out of the estimation because of the pandemic's extreme movements.",
            window[1], window[2])
  } else ""
  bslib::nav_panel_hidden(
    "about",
    htmltools::div(class = "about",
      htmltools::h2("The model"),
      htmltools::p(sprintf("A Bayesian vector autoregression (BVAR) of %d quarterly US series with %s lags, following ",
                           n_series, m$model$lags %||% "?"),
                   "Crump, Eusepi, Giannone, Qian and Sbordone (2021), “A Large Bayesian VAR of the ",
                   "United States Economy” (Federal Reserve Bank of New York Staff Report 976). The series ",
                   "enter in levels. The prior combines a Minnesota prior with a sum-of-coefficients prior; how ",
                   "tight they are is estimated from the data (Giannone, Lenza and Primiceri 2015)."),
      htmltools::p(sprintf("Estimated on %s–%s.%s Forecasts start from the latest data and the site shows %s ",
                           m$data$first_quarter, m$data$last_observed, window_text,
                           m$forecast$horizon_exported %||% "?"),
                   "quarters. Beyond two to three years a model in levels mostly extrapolates trends in its ",
                   "sample, so longer horizons are not shown."),
      htmltools::p(sprintf("Each forecast is simulated from %s posterior draws of the parameters, with %s paths ",
                           format(n_draws, big.mark = ","), n_paths),
                   "per draw. The bands are 68% and 90% Bayesian credible intervals: they cover parameter and ",
                   "shock uncertainty under the model, not the uncertainty about the model itself."),
      htmltools::h2("Scenarios"),
      htmltools::p("A scenario is a conditional forecast: what the model expects given the scenario's ",
                   "assumptions. Some scenarios fix the path of a few variables exactly (marked on the charts ",
                   "with hollow circles); others reweight the forecast so that its median meets a target, such ",
                   "as a survey's inflation forecast (Crump et al. 2021, Section 4.2; targets marked with ",
                   "hollow diamonds). Scenarios are listed, renamed and recoloured on the Settings tab. Scenarios are not predictions. The difference view shows one scenario's ",
                   "median minus the baseline's, a difference between two conditional forecasts, not a causal ",
                   "effect. Hover over a scenario's name for its assumptions."),
      htmltools::h2("Data"),
      htmltools::p(sprintf("Quarterly, %s–%s. FRED-QD vintage %s; series from the FRED API as of %s. ",
                           m$data$first_quarter, m$data$last_observed,
                           sub("^FRED-QD_(.*)\\.csv$", "\\1", basename(m$data$qd_vintage_file %||% "")),
                           m$data$api_realtime_date %||% ""),
                   "Monthly series are quarterly averages."),
      htmltools::div(class = "about-table-wrap", data_table),
      htmltools::h2("Sources and credits"),
      htmltools::p("Sources: U.S. Bureau of Economic Analysis, U.S. Bureau of Labor Statistics, U.S. Census ",
                   "Bureau, Board of Governors of the Federal Reserve System, via FRED and FRED-QD (Federal ",
                   "Reserve Bank of St. Louis); Moody's (corporate yields); University of Michigan Surveys of ",
                   "Consumers (sentiment); S&P Dow Jones Indices (S&P 500); World Bank Commodity Price Data ",
                   "(non-energy index, CC BY 4.0; used as a proxy for the S&P GSCI non-energy index); ",
                   "realized volatility computed from daily S&P 500 closes (Yahoo Finance)."),
      htmltools::p("This product uses the FRED® API but is not endorsed or certified by the Federal ",
                   "Reserve Bank of St. Louis.")
    )
  )
}

#' The whole page.
#'
#' @param bundle The bundle.
app_ui <- function(bundle) {
  m <- bundle$meta
  vintage <- sub("^FRED-QD_(.*)\\.csv$", "\\1", basename(m$data$qd_vintage_file))
  topic_panels <- lapply(APP_TOPICS, function(t) {
    bslib::nav_panel_hidden(t$id, shiny::uiOutput(paste0("grid_", t$id)), data_section(t$id))
  })
  panels <- c(list(overview_panel(bundle)), topic_panels,
              list(irf_panel(bundle),
                   settings_panel(), about_panel(bundle)))
  bslib::page(
    title = "US Economy - BVAR Forecasts",
    theme = bslib::bs_theme(version = 5, primary = "#0055A4", "font-family-sans-serif" = "'Public Sans', sans-serif"),
    htmltools::tags$head(
      htmltools::tags$link(rel = "stylesheet",
                           href = "https://fonts.googleapis.com/css2?family=Public+Sans:wght@400;500;600;700&display=swap"),
      # The file's modification time in the link, so a browser never keeps an old copy.
      htmltools::tags$link(rel = "stylesheet",
                           href = paste0("app.css?v=", as.integer(file.mtime(file.path("www", "app.css")))))
    ),
    htmltools::div(
      class = "app",
      side_panel(),
      htmltools::div(
        class = "main",
        settings_header(bundle),
        htmltools::tags$main(class = "content", do.call(bslib::navset_hidden, c(list(id = "pages"), panels))),
        htmltools::tags$footer(
          class = "footer",
          htmltools::span(paste("Data through", m$data$last_observed)),
          htmltools::span(paste("FRED-QD vintage", vintage)),
          htmltools::span(paste("Run", m$run_id)),
          htmltools::span("Credible bands, not confidence intervals; scenarios are conditional forecasts."),
          htmltools::span("This product uses the FRED® API but is not endorsed or certified by the ",
                          "Federal Reserve Bank of St. Louis.")
        )
      )
    )
  )
}
