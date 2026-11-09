# server.R: what happens when the reader changes a setting.
#
# Every chart is drawn from the bundle by fan_chart(); the server only
# picks the variable, transform and scenarios from the current settings.
# The settings live in the header, so they apply to every tab. Which
# scenarios the header lists, and their colours, are chosen on the
# Settings tab and last for the visit.


#' One chart, from the current settings.
#'
#' @param bundle The bundle.
#' @param variable Variable name.
#' @param s Settings (the server's settings()).
#' @param scenarios The listed scenarios (id, label, colour).
#' @return A plotly widget.
make_chart <- function(bundle, variable, s, scenarios) {
  tr <- display_transform(bundle$variables, variable, s$growth)
  # The difference view needs a baseline among the listed scenarios.
  difference <- identical(s$mode, "difference") && !is.null(s$baseline) && s$baseline %in% scenarios$id
  shown <- scenarios[scenarios$id %in% s$show, , drop = FALSE]
  # The difference view needs the baseline's median even if its line is off.
  ids <- if (difference) union(shown$id, s$baseline) else shown$id
  series <- chart_series(bundle, variable, tr, ids, APP_HISTORY[[s$history]])
  markers <- if (difference) NULL else chart_markers(bundle, variable, tr, shown$id)
  fan_chart(series, shown, s$bands, s$baseline, difference, markers)
}

#' The Data section under a tab's charts: a table and a CSV download.
#'
#' Collapsed by default; its table is computed only when opened.
#'
#' @param page The tab's id.
data_section <- function(page) {
  note <- "Numbers as the charts show them, in the current settings."
  if (page != "irf" && length(APP_NO_HISTORY_DOWNLOAD)) {
    # The restricted series, named from the settings (not fixed text).
    labels <- var_label(APP_NO_HISTORY_DOWNLOAD)
    listing <- if (length(labels) > 1) {
      paste0(paste(labels[-length(labels)], collapse = ", "), " and ", labels[length(labels)])
    } else labels
    note <- paste0(note, " The observed history of ", listing, " is left out, under their publishers' terms.")
  }
  htmltools::tags$details(
    class = "data-section",
    htmltools::tags$summary("Data"),
    htmltools::div(class = "data-actions",
      shiny::downloadButton(paste0("dl_", page), "Download CSV", class = "btn-download"),
      htmltools::span(class = "data-note", note)),
    htmltools::div(class = "data-table", shiny::tableOutput(paste0("tbl_", page)))
  )
}

#' A small table icon, for the chart's data button.
table_icon <- function() {
  htmltools::HTML(paste0(
    '<svg width="16" height="14" viewBox="0 0 16 14" aria-hidden="true" focusable="false">',
    '<rect x="0.75" y="0.75" width="14.5" height="12.5" rx="1" fill="none" stroke="currentColor" stroke-width="1.5"/>',
    '<line x1="0.75" y1="4.75" x2="15.25" y2="4.75" stroke="currentColor" stroke-width="1.5"/>',
    '<line x1="0.75" y1="9" x2="15.25" y2="9" stroke="currentColor" stroke-width="1"/>',
    '<line x1="6" y1="4.75" x2="6" y2="13.25" stroke="currentColor" stroke-width="1"/></svg>'))
}

#' A card holding one chart: title, unit, a data-table button and the plot.
#'
#' @param page,variable The tab and the variable (the plot's output id is
#'   c_<page>_<variable>, its table button data_<page>_<variable>).
#' @param unit The unit label shown beside the title.
chart_card <- function(page, variable, unit) {
  label <- var_label(variable)
  htmltools::tags$section(
    class = "chart-card",
    htmltools::div(class = "chart-head",
      htmltools::h2(label, htmltools::span(class = "unit", unit)),
      shiny::actionButton(paste0("data_", page, "_", variable), label = NULL, icon = table_icon(),
                          class = "chart-data-btn", title = paste("Data table:", label),
                          `aria-label` = paste("Data table for", label))),
    plotly::plotlyOutput(paste0("c_", page, "_", variable), height = "250px")
  )
}

#' The pop-up with one variable's numbers: a table and a CSV download.
#'
#' @param variable Variable name.
#' @param restricted TRUE if the table leaves out the observed history.
data_modal <- function(variable, restricted = FALSE) {
  shiny::modalDialog(
    title = var_label(variable), size = "l", easyClose = TRUE, footer = shiny::modalButton("Close"),
    htmltools::div(class = "data-actions",
      shiny::downloadButton("modal_dl", "Download CSV", class = "btn-download"),
      htmltools::span(class = "data-note", "Numbers as the chart shows them, in the current settings.",
                      if (restricted) " The observed history is left out, under its publisher's terms.")),
    htmltools::div(class = "data-table", shiny::tableOutput("modal_tbl"))
  )
}

#' The server function.
#'
#' @param bundle The bundle.
app_server <- function(bundle) {
  catalog <- scenario_catalog(bundle)
  default_listed <- catalog$id[catalog$default]
  default_colours <- stats::setNames(catalog$colour, catalog$id)
  default_shock <- irf_shocks(bundle)$id[1]

  function(input, output, session) {
    shiny::observeEvent(input$page, bslib::nav_select("pages", input$page, session = session))
    shiny::observeEvent(input$overview_reset, {
      shiny::updateSelectizeInput(session, "overview_vars", selected = APP_OVERVIEW_DEFAULT)
    })

    # ---- Settings tab: which scenarios are listed, their names and colours ----
    default_names <- stats::setNames(catalog$label, catalog$id)
    picks <- shiny::reactiveValues(listed = default_listed, colours = default_colours, names = default_names,
                                   version = 0)
    for (sid in catalog$id) local({
      id <- sid
      shiny::observeEvent(input[[paste0("list_", id)]], {
        on <- isTRUE(input[[paste0("list_", id)]])
        listed <- if (on) union(picks$listed, id) else setdiff(picks$listed, id)
        picks$listed <- catalog$id[catalog$id %in% listed]          # keep the catalog's order
      })
      shiny::observeEvent(input[[paste0("colour_", id)]], {
        colours <- picks$colours
        colours[[id]] <- input[[paste0("colour_", id)]]
        picks$colours <- colours
      })
      # A blank name falls back to the scenario's own short label.
      shiny::observeEvent(input[[paste0("name_", id)]], {
        name <- trimws(input[[paste0("name_", id)]])
        nm <- picks$names
        nm[[id]] <- if (nzchar(name)) name else default_names[[id]]
        picks$names <- nm
      })
    })
    shiny::observeEvent(input$reset_scenarios, {
      picks$listed <- default_listed
      picks$colours <- default_colours
      picks$names <- default_names
      picks$version <- picks$version + 1
    })
    # Names change as the reader types: wait for a pause before redrawing.
    names_now <- shiny::debounce(shiny::reactive(picks$names), 500)
    listed <- shiny::reactive({
      x <- catalog[catalog$id %in% picks$listed, c("id", "label", "colour", "description"), drop = FALSE]
      x$colour <- unname(picks$colours[x$id])
      x$label <- unname(names_now()[x$id])
      x
    })

    output$scenario_table <- shiny::renderUI({
      picks$version                                # rebuilt on reset only
      on <- shiny::isolate(picks$listed)
      colours <- shiny::isolate(picks$colours)
      nm <- shiny::isolate(picks$names)
      htmltools::tags$table(
        class = "scenario-table",
        htmltools::tags$thead(htmltools::tags$tr(htmltools::tags$th(class = "col-narrow", "Listed"),
                                                 htmltools::tags$th(class = "col-narrow", "Colour"),
                                                 htmltools::tags$th(class = "col-name", "Name"),
                                                 htmltools::tags$th("Assumptions"))),
        htmltools::tags$tbody(lapply(seq_len(nrow(catalog)), function(k) {
          id <- catalog$id[k]
          scenario_settings_row(catalog[k, ], id %in% on, colours[[id]], nm[[id]])
        }))
      )
    })

    # The header's scenario list, rebuilt when the listed scenarios or their
    # colours change; the reader's ticks are kept, and a newly listed
    # scenario starts shown.
    known <- shiny::reactiveVal(NULL)
    output$scenario_controls <- shiny::renderUI({
      x <- listed()
      before <- shiny::isolate(known())
      show <- if (is.null(before)) x$id else union(intersect(shiny::isolate(input$show), x$id), setdiff(x$id, before))
      bands <- if (is.null(before)) x$id[1] else intersect(shiny::isolate(input$bands), x$id)
      baseline <- shiny::isolate(input$baseline)
      if (is.null(baseline) || !baseline %in% x$id) baseline <- x$id[1]
      known(x$id)
      # Bands are not drawn in the difference view, so the column is hidden there.
      scenario_controls(x, show, bands, baseline, show_bands = !identical(input$mode, "difference"))
    })

    settings <- shiny::reactive({
      list(show = input$show, bands = input$bands, baseline = input$baseline %||% listed()$id[1],
           mode = input$mode %||% "forecast", growth = input$growth %||% "yoy",
           history = input$history %||% names(APP_HISTORY)[1])
    })

    # ---- Charts ----
    pages <- c(list(overview = shiny::reactive(input$overview_vars)),
               stats::setNames(lapply(APP_TOPICS, function(t) shiny::reactive(t$variables)),
                               vapply(APP_TOPICS, `[[`, "", "id")))
    all_vars <- intersect(APP_VARIABLES$variable, bundle$variables$variable)

    # ---- Impulse-response settings (the header on the IRF tab) ----
    irf_settings <- shiny::reactive({
      size <- input$irf_size
      list(shock = input$irf_shock %||% default_shock, units = input$irf_units %||% "level",
           horizon = as.numeric(input$irf_horizon %||% max(APP_IRF_HORIZONS)),
           size = if (is.numeric(size) && length(size) == 1 && is.finite(size) && size != 0) size else 1)
    })
    irf_defaults <- function() {
      shiny::updateSelectizeInput(session, "irf_vars",
                                  selected = irf_default_responses(bundle, irf_settings()$shock))
    }
    shiny::observeEvent(input$irf_shock, irf_defaults(), ignoreInit = TRUE)
    shiny::observeEvent(input$irf_reset, irf_defaults())
    output$irf_size_unit <- shiny::renderUI({
      htmltools::span(irf_shock_unit(bundle$variables, irf_settings()$shock))
    })

    # The table pop-up: one variable's numbers and CSV, in the current settings
    # (a forecast chart's or a response chart's).
    modal <- shiny::reactiveVal(NULL)                       # list(variable, kind)
    modal_variable <- shiny::reactive(modal()$variable)
    modal_data <- shiny::reactive({
      shiny::req(modal())
      m <- modal()
      if (m$kind == "irf") {
        s <- irf_settings()
        irf_data(bundle, s$shock, m$variable, s$units, s$horizon, s$size)
      } else {
        tab_data(bundle, m$variable, settings(), listed())
      }
    })
    output$modal_tbl <- shiny::renderTable(modal_data(), digits = 3, na = "", striped = TRUE)
    output$modal_dl <- shiny::downloadHandler(
      filename = function() paste0("bvar_", modal()$kind, "_", modal_variable(), "_", bundle$meta$data$last_observed,
                                   ".csv"),
      content = function(file) utils::write.csv(modal_data(), file, row.names = FALSE, na = ""))

    # ---- Impulse-response charts ----
    output$grid_irf <- shiny::renderUI({
      s <- irf_settings()
      htmltools::div(class = "chart-grid", lapply(input$irf_vars, function(v) {
        tr <- if (bundle$variables$model_transform[match(v, bundle$variables$variable)] == "log") s$units else "level"
        chart_card("irf", v, irf_unit(bundle$variables, v, tr))
      }))
    })
    for (vv in all_vars) local({
      v <- vv
      output[[paste0("c_irf_", v)]] <- plotly::renderPlotly({
        s <- irf_settings()
        irf_chart(irf_series(bundle, s$shock, v, s$units, s$horizon, s$size))
      })
      shiny::observeEvent(input[[paste0("data_irf_", v)]], {
        modal(list(variable = v, kind = "irf"))
        shiny::showModal(data_modal(v))
      })
    })
    output$tbl_irf <- shiny::renderTable({
      s <- irf_settings()
      irf_data(bundle, s$shock, input$irf_vars, s$units, s$horizon, s$size)
    }, digits = 3, na = "", striped = TRUE)
    output$dl_irf <- shiny::downloadHandler(
      filename = function() paste0("bvar_irf_", irf_settings()$shock, "_", bundle$meta$data$last_observed, ".csv"),
      content = function(file) {
        s <- irf_settings()
        utils::write.csv(irf_data(bundle, s$shock, input$irf_vars, s$units, s$horizon, s$size), file,
                         row.names = FALSE, na = "")
      })

    # ---- Forecast charts ----
    for (pg in names(pages)) local({
      page <- pg
      vars_of_page <- if (page == "overview") all_vars else APP_TOPICS[[match(page, names(pages)) - 1]]$variables
      # The grid of cards: rebuilt when the charts or their units change.
      output[[paste0("grid_", page)]] <- shiny::renderUI({
        vars <- pages[[page]]()
        difference <- identical(input$mode, "difference")
        growth <- input$growth %||% "yoy"
        htmltools::div(class = "chart-grid", lapply(vars, function(v) {
          chart_card(page, v, unit_label(v, display_transform(bundle$variables, v, growth), difference))
        }))
      })
      # One plot output per variable the page can show; only the visible ones
      # are computed. Its table button opens the variable's numbers.
      for (vv in vars_of_page) local({
        v <- vv
        output[[paste0("c_", page, "_", v)]] <- plotly::renderPlotly(make_chart(bundle, v, settings(), listed()))
        shiny::observeEvent(input[[paste0("data_", page, "_", v)]], {
          modal(list(variable = v, kind = "forecast"))
          shiny::showModal(data_modal(v, restricted = v %in% APP_NO_HISTORY_DOWNLOAD))
        })
      })
      # The page's numbers: a table (computed when the Data section is open) and a CSV.
      output[[paste0("tbl_", page)]] <- shiny::renderTable(
        tab_data(bundle, pages[[page]](), settings(), listed()), digits = 3, na = "", striped = TRUE)
      output[[paste0("dl_", page)]] <- shiny::downloadHandler(
        filename = function() paste0("bvar_", page, "_data_", bundle$meta$data$last_observed, ".csv"),
        content = function(file) {
          utils::write.csv(tab_data(bundle, pages[[page]](), settings(), listed()), file, row.names = FALSE, na = "")
        })
    })
  }
}
