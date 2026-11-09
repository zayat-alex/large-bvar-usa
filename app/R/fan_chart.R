# fan_chart.R: the one chart function the site uses.
#
# Draws a chart_series() (chart_data.R) with plotly. Forecast view: the
# observed history, then each scenario's median; with bands switched on,
# the baseline's 68% and 90% bands as shaded areas and every other
# scenario's as lines (solid 68%, dotted 90%), so several sets of bands
# stay readable. Difference view: each scenario's median minus the
# baseline's, around a zero line. A dashed vertical line marks the last
# observed quarter. Each chart is its own plotly widget, so it can be
# zoomed and panned on its own (double-click resets it).


#' A colour as rgba() with transparency.
#'
#' @param hex Colour as "#RRGGBB".
#' @param alpha Opacity between 0 and 1.
#' @return Character "rgba(r,g,b,alpha)".
hex_alpha <- function(hex, alpha) {
  rgb <- grDevices::col2rgb(hex)[, 1]
  sprintf("rgba(%d,%d,%d,%.2f)", rgb[1], rgb[2], rgb[3], alpha)
}

#' Hover text of a line: "<label> 2027Q4: 2.49".
hover_text <- function(label, quarter, value) sprintf("%s %s: %.2f", label, quarter, value)

#' Fan chart of one variable.
#'
#' @param series Output of chart_series().
#' @param scenarios The rows of APP_SCENARIOS to draw (id, label, colour),
#'   in drawing order.
#' @param bands Ids of the scenarios whose bands are drawn.
#' @param baseline Id of the baseline.
#' @param difference TRUE for the difference view.
#' @param markers NULL or chart_markers(): cells a scenario fixes (hollow
#'   circles) and soft targets (hollow diamonds), drawn in the forecast view.
#' @return A plotly htmlwidget.
fan_chart <- function(series, scenarios, bands, baseline, difference = FALSE, markers = NULL) {
  axis_style <- list(showline = TRUE, linecolor = "#6F7680", linewidth = 1, ticks = "outside",
                     tickcolor = "#6F7680", ticklen = 4, gridcolor = "#ECEEF1", zeroline = FALSE,
                     tickfont = list(size = 11, color = "#4F5660"))
  last_date <- series$last$date
  p <- plotly::plot_ly()
  shapes <- list(list(type = "line", xref = "x", yref = "paper", x0 = last_date, x1 = last_date, y0 = 0, y1 = 1,
                      line = list(color = "#8A9099", width = 1, dash = "dot")))

  if (!difference) {
    # Bands first, so the medians are drawn on top of them.
    for (k in seq_len(nrow(scenarios))) {
      id <- scenarios$id[k]
      if (!id %in% bands) next
      f <- series$forecasts[[id]]
      col <- scenarios$colour[k]
      if (id == baseline) {
        p <- plotly::add_ribbons(p, data = f, x = ~date, ymin = ~q05, ymax = ~q95, fillcolor = hex_alpha(col, 0.14),
                                 line = list(width = 0), hoverinfo = "skip", showlegend = FALSE)
        p <- plotly::add_ribbons(p, data = f, x = ~date, ymin = ~q16, ymax = ~q84, fillcolor = hex_alpha(col, 0.28),
                                 line = list(width = 0), hoverinfo = "skip", showlegend = FALSE)
      } else {
        for (q in c("q05", "q95", "q16", "q84")) {
          dash <- if (q %in% c("q05", "q95")) "dot" else "solid"
          p <- plotly::add_lines(p, x = f$date, y = f[[q]], line = list(color = col, width = 1, dash = dash),
                                 hoverinfo = "skip", showlegend = FALSE)
        }
      }
    }
    h <- series$history
    p <- plotly::add_lines(p, x = h$date, y = h$value, line = list(color = "#1B1F24", width = 1.8),
                           text = hover_text("Observed", h$quarter, h$value), hoverinfo = "text", showlegend = FALSE)
    for (k in seq_len(nrow(scenarios))) {
      f <- series$forecasts[[scenarios$id[k]]]
      p <- plotly::add_lines(p, x = f$date, y = f$q50, line = list(color = scenarios$colour[k], width = 2),
                             text = hover_text(scenarios$label[k], f$quarter, f$q50), hoverinfo = "text",
                             showlegend = FALSE)
    }
    # What each scenario assumes: fixed cells and soft targets, on top.
    for (k in seq_len(nrow(scenarios))) {
      m <- if (is.null(markers)) NULL else markers[markers$object == scenarios$id[k], , drop = FALSE]
      if (is.null(m) || !nrow(m)) next
      col <- scenarios$colour[k]
      text <- ifelse(m$kind == "fixed",
                     sprintf("%s %s: %.2f, fixed by the scenario", scenarios$label[k], m$quarter, m$value),
                     sprintf("%s %s: target median %.2f", scenarios$label[k], m$quarter, m$value))
      p <- plotly::add_markers(p, x = m$date, y = m$value, text = text, hoverinfo = "text", showlegend = FALSE,
                               marker = list(symbol = ifelse(m$kind == "fixed", "circle-open", "diamond-open"),
                                             size = 9, color = col, line = list(color = col, width = 1.8)))
    }
  } else {
    d <- median_differences(series, baseline)
    shapes <- c(shapes, list(list(type = "line", xref = "paper", yref = "y", x0 = 0, x1 = 1, y0 = 0, y1 = 0,
                                  line = list(color = "#6F7680", width = 1))))
    for (k in seq_len(nrow(scenarios))) {
      x <- d[[scenarios$id[k]]]
      if (is.null(x)) next
      p <- plotly::add_lines(p, x = x$date, y = x$value, line = list(color = scenarios$colour[k], width = 2),
                             text = hover_text(scenarios$label[k], x$quarter, x$value), hoverinfo = "text",
                             showlegend = FALSE)
    }
  }

  p <- plotly::layout(p, shapes = shapes, hovermode = "closest", showlegend = FALSE,
                      margin = list(l = 44, r = 10, t = 8, b = 28),
                      font = list(family = "Public Sans, sans-serif", size = 11, color = "#1B1F24"),
                      xaxis = c(axis_style, list(type = "date", title = list(text = ""))),
                      yaxis = c(axis_style, list(title = list(text = ""))),
                      paper_bgcolor = "#FFFFFF", plot_bgcolor = "#FFFFFF")
  plotly::config(p, displaylogo = FALSE,
                 modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d", "hoverClosestCartesian",
                                            "hoverCompareCartesian", "toggleSpikelines"))
}

#' Impulse-response chart: median with shaded 68% and 90% bands.
#'
#' @param x Output of irf_series().
#' @param colour Colour of the response.
#' @return A plotly htmlwidget; x axis in quarters after the shock, a line
#'   at zero (no response).
irf_chart <- function(x, colour = "#0055A4") {
  axis_style <- list(showline = TRUE, linecolor = "#6F7680", linewidth = 1, ticks = "outside",
                     tickcolor = "#6F7680", ticklen = 4, gridcolor = "#ECEEF1", zeroline = FALSE,
                     tickfont = list(size = 11, color = "#4F5660"))
  p <- plotly::plot_ly()
  p <- plotly::add_ribbons(p, data = x, x = ~horizon, ymin = ~q05, ymax = ~q95, fillcolor = hex_alpha(colour, 0.14),
                           line = list(width = 0), hoverinfo = "skip", showlegend = FALSE)
  p <- plotly::add_ribbons(p, data = x, x = ~horizon, ymin = ~q16, ymax = ~q84, fillcolor = hex_alpha(colour, 0.28),
                           line = list(width = 0), hoverinfo = "skip", showlegend = FALSE)
  p <- plotly::add_lines(p, x = x$horizon, y = x$q50, line = list(color = colour, width = 2), hoverinfo = "text",
                         text = sprintf("Quarter %d: %.3f (68%%: %.3f to %.3f)", x$horizon, x$q50, x$q16, x$q84),
                         showlegend = FALSE)
  zero <- list(type = "line", xref = "paper", yref = "y", x0 = 0, x1 = 1, y0 = 0, y1 = 0,
               line = list(color = "#6F7680", width = 1))
  p <- plotly::layout(p, shapes = list(zero), hovermode = "closest", showlegend = FALSE,
                      margin = list(l = 44, r = 10, t = 8, b = 36),
                      font = list(family = "Public Sans, sans-serif", size = 11, color = "#1B1F24"),
                      xaxis = c(axis_style, list(title = list(text = "quarters after the shock", font = list(size = 11)),
                                                 dtick = 4)),
                      yaxis = c(axis_style, list(title = list(text = ""))),
                      paper_bgcolor = "#FFFFFF", plot_bgcolor = "#FFFFFF")
  plotly::config(p, displaylogo = FALSE,
                 modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d", "hoverClosestCartesian",
                                            "hoverCompareCartesian", "toggleSpikelines"))
}
