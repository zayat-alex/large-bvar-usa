// site.js: the static website, drawn in the browser.
//
// The same site as the Shiny app (app/), with the same numbers: the app's
// R code reads the bundle and builds plotly charts on a server; here the
// same steps run in the browser on data/core.json and data/irf.json,
// written from the bundle by build_site.R. Each function below names the
// app function it follows (app/R/chart_data.R, fan_chart.R, settings.R,
// ui.R, server.R), so the two can be compared line by line.
//
// What a chart shows (as in the app): one variable in one transform, its
// recent history, then for each scenario switched on the forecast median
// and optionally the 68% and 90% credible bands. Every forecast starts at
// the last observed quarter, at the observed value. In the difference view
// each line is a scenario's median minus the baseline's median: a
// difference between two conditional forecasts, not a causal effect.

"use strict";

// ---- Data and state ----

let CORE = null;          // data/core.json
let IRF = null;           // data/irf.json, loaded when the IRF tab first opens

// The reader's choices, for the visit (the app's inputs).
const state = {
  page: "overview",
  listed: [],             // scenarios listed in the header (Settings tab)
  known: [],              // the listed scenarios when the header was last built
  colours: {},            // scenario id -> colour
  names: {},              // scenario id -> name
  show: [], bands: [], baseline: null,
  mode: "forecast", growth: "yoy", history: "5y",
  overviewVars: [],
};

// ---- Small helpers ----

// Element builder: el("div", {class: "x"}, child, ...).
function el(tag, attrs, ...children) {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v === null || v === undefined || v === false) continue;
    if (k.startsWith("on")) e.addEventListener(k.slice(2), v);
    else if (v === true) e.setAttribute(k, "");
    else e.setAttribute(k, v);
  }
  for (const c of children.flat()) {
    if (c === null || c === undefined || c === false) continue;
    e.append(c instanceof Node ? c : document.createTextNode(String(c)));
  }
  return e;
}

// Quarter label to the first day of the quarter ("2026Q3" -> "2026-07-01"),
// as quarter_date() and the bundle's date column.
function quarterDate(q) {
  const m = 3 * Number(q.slice(5, 6)) - 2;
  return `${q.slice(0, 4)}-${String(m).padStart(2, "0")}-01`;
}

// Consecutive quarters from a first one ("1973Q1", n) -> labels.
function quarterSeq(start, n) {
  let k = Number(start.slice(0, 4)) * 4 + Number(start.slice(5, 6)) - 1;
  const out = [];
  for (let i = 0; i < n; i++, k++) out.push(`${Math.floor(k / 4)}Q${(k % 4) + 1}`);
  return out;
}

// A colour as rgba() with transparency (hex_alpha()).
function hexAlpha(hex, alpha) {
  const n = parseInt(hex.slice(1), 16);
  return `rgba(${(n >> 16) & 255},${(n >> 8) & 255},${n & 255},${alpha.toFixed(2)})`;
}

const variableRow = (v) => CORE.variables.find((x) => x.variable === v);
const varLabel = (v) => variableRow(v).label;

// ---- Settings (app/R/settings.R) ----

// Unit label of a chart (unit_label()).
function unitLabel(variable, transform, difference) {
  const unit = variableRow(variable).unit;
  if (transform === "yoy" || transform === "qoq_ann") {
    const growth = transform === "yoy" ? "% y/y" : "% q/q ann.";
    return difference ? `pp (${growth})` : growth;
  }
  if (difference && ["%", "%, annualized", "pp"].includes(unit)) return "pp";
  return unit;
}

// History window in quarters (APP_HISTORY); Infinity for all.
function historyQuarters(name) {
  const h = CORE.settings.history.find((x) => x.name === name);
  return h.quarters === null ? Infinity : h.quarters;
}

// ---- Chart data (app/R/chart_data.R) ----

// Transform in which a variable is shown (display_transform()).
function displayTransform(variable, growth) {
  const d = variableRow(variable).display;
  return d === "growth" ? growth : d;
}

// History and forecasts of one variable, ready to plot (chart_series()).
// Returns {history: [{quarter, date, value}], last, forecasts: {id: rows}},
// each forecast starting with the last observation as its anchor.
function chartSeries(variable, transform, objects, historyQ) {
  const h = CORE.history[variable][transform];
  const quarters = quarterSeq(h.start, h.values.length);
  let history = quarters.map((q, i) => ({ quarter: q, date: quarterDate(q), value: h.values[i] }));
  if (Number.isFinite(historyQ)) history = history.slice(-historyQ);
  const last = history[history.length - 1];
  const fq = CORE.meta.forecast_quarters;
  const forecasts = {};
  for (const obj of objects) {
    const f = CORE.forecasts[obj][variable][transform];
    const rows = [{ quarter: last.quarter, date: last.date, q05: last.value, q16: last.value, q50: last.value,
                    q84: last.value, q95: last.value }];
    fq.forEach((q, i) => rows.push({ quarter: q, date: quarterDate(q), q05: f.q05[i], q16: f.q16[i],
                                     q50: f.q50[i], q84: f.q84[i], q95: f.q95[i] }));
    forecasts[obj] = rows;
  }
  return { history, last, forecasts };
}

// Cells a scenario fixes and median targets it meets, on this chart
// (chart_markers()): [{object, quarter, date, value, kind}].
function chartMarkers(variable, transform, ids) {
  const pick = (rows, kind) => rows
    .filter((r) => r.variable === variable && r.transform === transform && ids.includes(r.object))
    .map((r) => ({ object: r.object, quarter: r.quarter, date: quarterDate(r.quarter), value: r.value, kind }));
  return [...pick(CORE.conditions, "fixed"), ...pick(CORE.targets, "target")];
}

// Medians minus the baseline's median (median_differences()).
function medianDifferences(series, baseline) {
  const base = series.forecasts[baseline];
  const out = {};
  for (const [obj, rows] of Object.entries(series.forecasts)) {
    if (obj === baseline) continue;
    out[obj] = rows.map((r, i) => ({ quarter: r.quarter, date: r.date, value: r.q50 - base[i].q50 }));
  }
  return out;
}

// ---- Charts (app/R/fan_chart.R) ----

// automargin is not set in fan_chart.R, but the R plotly package turns it
// on for every axis by default, so it is set here to draw the same charts.
const AXIS_STYLE = { showline: true, linecolor: "#6F7680", linewidth: 1, ticks: "outside", tickcolor: "#6F7680",
                     ticklen: 4, gridcolor: "#ECEEF1", zeroline: false,
                     tickfont: { size: 11, color: "#4F5660" }, automargin: true };
const PLOT_CONFIG = { displaylogo: false, responsive: true,
                      modeBarButtonsToRemove: ["lasso2d", "select2d", "autoScale2d", "hoverClosestCartesian",
                                               "hoverCompareCartesian", "toggleSpikelines"] };

// Hover text of a line: "<label> 2027Q4: 2.49" (hover_text()).
const hoverText = (label, quarter, value) =>
  `${label} ${quarter}: ${value === null ? "NA" : value.toFixed(2)}`;

// A band between two quantiles as a filled polygon (plotly's add_ribbons()).
function ribbon(rows, lo, hi, fill) {
  return { type: "scatter", mode: "lines", x: [...rows.map((r) => r.date), ...rows.map((r) => r.date).reverse()],
           y: [...rows.map((r) => r[hi]), ...rows.map((r) => r[lo]).reverse()], fill: "toself", fillcolor: fill,
           line: { width: 0 }, hoverinfo: "skip", showlegend: false };
}

function line(x, y, style, text) {
  const t = { type: "scatter", mode: "lines", x, y, line: style, showlegend: false };
  if (text) { t.text = text; t.hoverinfo = "text"; } else t.hoverinfo = "skip";
  return t;
}

// Fan chart of one variable (fan_chart()): traces and layout.
function fanChart(series, scenarios, bands, baseline, difference, markers) {
  const traces = [];
  const lastDate = series.last.date;
  const shapes = [{ type: "line", xref: "x", yref: "paper", x0: lastDate, x1: lastDate, y0: 0, y1: 1,
                    line: { color: "#8A9099", width: 1, dash: "dot" } }];
  if (!difference) {
    // Bands first, so the medians are drawn on top of them.
    for (const s of scenarios) {
      if (!bands.includes(s.id)) continue;
      const f = series.forecasts[s.id];
      if (s.id === baseline) {
        traces.push(ribbon(f, "q05", "q95", hexAlpha(s.colour, 0.14)));
        traces.push(ribbon(f, "q16", "q84", hexAlpha(s.colour, 0.28)));
      } else {
        for (const q of ["q05", "q95", "q16", "q84"]) {
          const dash = q === "q05" || q === "q95" ? "dot" : "solid";
          traces.push(line(f.map((r) => r.date), f.map((r) => r[q]), { color: s.colour, width: 1, dash }));
        }
      }
    }
    const h = series.history;
    traces.push(line(h.map((r) => r.date), h.map((r) => r.value), { color: "#1B1F24", width: 1.8 },
                     h.map((r) => hoverText("Observed", r.quarter, r.value))));
    for (const s of scenarios) {
      const f = series.forecasts[s.id];
      traces.push(line(f.map((r) => r.date), f.map((r) => r.q50), { color: s.colour, width: 2 },
                       f.map((r) => hoverText(s.label, r.quarter, r.q50))));
    }
    // What each scenario assumes: fixed cells and soft targets, on top.
    for (const s of scenarios) {
      const m = (markers || []).filter((r) => r.object === s.id);
      if (!m.length) continue;
      traces.push({
        type: "scatter", mode: "markers", x: m.map((r) => r.date), y: m.map((r) => r.value), hoverinfo: "text",
        showlegend: false,
        text: m.map((r) => r.kind === "fixed"
          ? `${s.label} ${r.quarter}: ${r.value.toFixed(2)}, fixed by the scenario`
          : `${s.label} ${r.quarter}: target median ${r.value.toFixed(2)}`),
        marker: { symbol: m.map((r) => (r.kind === "fixed" ? "circle-open" : "diamond-open")), size: 9,
                  color: s.colour, line: { color: s.colour, width: 1.8 } },
      });
    }
  } else {
    const d = medianDifferences(series, baseline);
    shapes.push({ type: "line", xref: "paper", yref: "y", x0: 0, x1: 1, y0: 0, y1: 0,
                  line: { color: "#6F7680", width: 1 } });
    for (const s of scenarios) {
      const x = d[s.id];
      if (!x) continue;
      traces.push(line(x.map((r) => r.date), x.map((r) => r.value), { color: s.colour, width: 2 },
                       x.map((r) => hoverText(s.label, r.quarter, r.value))));
    }
  }
  const layout = {
    shapes, hovermode: "closest", showlegend: false, margin: { l: 44, r: 10, t: 8, b: 28 },
    font: { family: "Public Sans, sans-serif", size: 11, color: "#1B1F24" },
    xaxis: { ...AXIS_STYLE, type: "date", title: { text: "" } },
    yaxis: { ...AXIS_STYLE, title: { text: "" } },
    paper_bgcolor: "#FFFFFF", plot_bgcolor: "#FFFFFF",
  };
  return { traces, layout };
}

// ---- Current settings (server.R) ----

// The listed scenarios with their current names and colours (listed()).
function listedScenarios() {
  return CORE.scenarios.filter((s) => state.listed.includes(s.id))
    .map((s) => ({ id: s.id, label: state.names[s.id], colour: state.colours[s.id], description: s.description }));
}

// The difference view needs a baseline among the listed scenarios.
function isDifference() {
  return state.mode === "difference" && state.baseline !== null && state.listed.includes(state.baseline);
}

// One chart from the current settings (make_chart()).
function makeChart(variable) {
  const tr = displayTransform(variable, state.growth);
  const difference = isDifference();
  const shown = listedScenarios().filter((s) => state.show.includes(s.id));
  const ids = shown.map((s) => s.id);
  if (difference && !ids.includes(state.baseline)) ids.push(state.baseline);
  const series = chartSeries(variable, tr, ids, historyQuarters(state.history));
  const markers = difference ? null : chartMarkers(variable, tr, shown.map((s) => s.id));
  return fanChart(series, shown, state.bands, state.baseline, difference, markers);
}

// ---- Tables and downloads (chart_data.R tab_data(), server.R) ----

const BAND_COLUMNS = ["68% low", "68% high", "90% low", "90% high"];
const NUMERIC_COLUMNS = ["Value", "Median", "Quarters after shock", ...BAND_COLUMNS];
const INTEGER_COLUMNS = ["Quarters after shock"];               // printed without decimals, as R's integers

// The numbers of a tab's charts, as one table (tab_data()). Forecast view:
// for each variable the observed history in the chart's window, then each
// shown scenario's median and bands. Difference view: each shown
// scenario's median minus the baseline's. The observed history of the
// variables in no_history_download is left out (the publishers' terms).
// Returns {columns, rows}, rows being objects keyed by column.
function tabData(variables) {
  const difference = isDifference();
  const scenarios = listedScenarios();
  const shown = scenarios.filter((s) => state.show.includes(s.id));
  const restricted = CORE.settings.no_history_download;
  const rows = [];
  for (const v of variables) {
    const tr = displayTransform(v, state.growth);
    const ids = shown.map((s) => s.id);
    if (difference && !ids.includes(state.baseline)) ids.push(state.baseline);
    const series = chartSeries(v, tr, ids, historyQuarters(state.history));
    const base = { Variable: varLabel(v), Unit: unitLabel(v, tr, difference) };
    if (difference) {
      const baseLabel = scenarios.find((s) => s.id === state.baseline).label;
      const d = medianDifferences(series, state.baseline);
      for (const s of shown) {
        if (!d[s.id]) continue;
        for (const r of d[s.id].slice(1)) {                      // drop the anchor (zero)
          rows.push({ ...base, Series: `${s.label} minus ${baseLabel}`, Quarter: r.quarter, Value: r.value });
        }
      }
    } else {
      if (!restricted.includes(v)) {
        for (const r of series.history) {
          rows.push({ ...base, Series: "Observed", Quarter: r.quarter, Value: r.value, "68% low": null,
                      "68% high": null, "90% low": null, "90% high": null });
        }
      }
      for (const s of shown) {
        for (const r of series.forecasts[s.id].slice(1)) {       // drop the anchor (the last observation)
          rows.push({ ...base, Series: `${s.label} (median)`, Quarter: r.quarter, Value: r.q50, "68% low": r.q16,
                      "68% high": r.q84, "90% low": r.q05, "90% high": r.q95 });
        }
      }
    }
  }
  const columns = ["Variable", "Unit", "Series", "Quarter", "Value"];
  if (!difference && rows.length) columns.push(...BAND_COLUMNS);
  return { columns, rows };
}

const isMissing = (x) => x === null || x === undefined || Number.isNaN(x);

// A number with three decimals, as R's formatC(x, format = "f", digits = 3)
// prints it: like toFixed(3), except that R keeps the sign of a negative
// zero (a zero response scaled by a negative shock size).
const fixed3 = (x) => (Object.is(x, -0) ? "-" : "") + x.toFixed(3);

// A number as R's write.csv() writes it: 15 significant digits, in fixed
// or scientific notation, whichever is shorter (fixed on a tie); the
// exponent has a sign and at least two digits ("1e-05", "2.5e+07").
function rNumber(x) {
  if (x === 0) return "0";
  const p = Number(x.toPrecision(15));
  const [mantissa, exp] = p.toExponential().split("e");
  const e = Number(exp);
  const sci = `${mantissa}e${e < 0 ? "-" : "+"}${String(Math.abs(e)).padStart(2, "0")}`;
  const fixed = String(p);
  if (fixed.includes("e")) return sci;                          // beyond 1e21 or below 1e-6: fixed is longer
  return fixed.length <= sci.length ? fixed : sci;
}

// A table as Shiny's renderTable(digits = 3, na = "", striped = TRUE)
// draws it: numbers with three decimals, right-aligned; missing values empty.
function tableElement({ columns, rows }) {
  const numeric = (c) => NUMERIC_COLUMNS.includes(c);
  const head = el("tr", {}, columns.map((c) => el("th", { style: `text-align: ${numeric(c) ? "right" : "left"};` }, c)));
  const body = rows.map((r) => el("tr", {}, columns.map((c) => {
    if (!numeric(c)) return el("td", {}, r[c]);
    if (isMissing(r[c])) return el("td", { align: "right", class: "NA" });
    return el("td", { align: "right" }, INTEGER_COLUMNS.includes(c) ? String(r[c]) : fixed3(r[c]));
  })));
  return el("table", { class: "table shiny-table table-striped spacing-s", style: "width:auto;" },
            el("thead", {}, head), el("tbody", {}, body));
}

// A table as CSV, as utils::write.csv(row.names = FALSE, na = ""): text
// quoted, numbers as rNumber() writes them, missing values empty.
function toCsv({ columns, rows }) {
  const quote = (s) => `"${String(s).replace(/"/g, '""')}"`;
  const lines = [columns.map(quote).join(",")];
  for (const r of rows) {
    lines.push(columns.map((c) => {
      if (NUMERIC_COLUMNS.includes(c)) return isMissing(r[c]) ? "" : rNumber(r[c]);
      return quote(r[c]);
    }).join(","));
  }
  return lines.join("\n") + "\n";
}

// Save a table as a CSV file in the reader's browser.
function downloadCsv(filename, table) {
  const url = URL.createObjectURL(new Blob([toCsv(table)], { type: "text/csv;charset=utf-8" }));
  const a = el("a", { href: url, download: filename });
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

// The download button (shiny::downloadButton()).
function downloadButton(id, onClick) {
  return el("a", { class: "btn btn-default shiny-download-link btn-download", download: true, href: "#", id,
                   onclick: (e) => { e.preventDefault(); onClick(); } },
            el("i", { class: "fas fa-download", role: "presentation", "aria-label": "download icon" }), " Download CSV");
}

// Labels joined as "a, b and c".
function listing(labels) {
  return labels.length > 1 ? `${labels.slice(0, -1).join(", ")} and ${labels[labels.length - 1]}` : labels.join("");
}

// The Data section under a tab's charts (data_section()): collapsed; its
// table is computed when opened, and redrawn while open.
function dataSection(page) {
  let note = "Numbers as the charts show them, in the current settings.";
  const restricted = CORE.settings.no_history_download;
  if (page !== "irf" && restricted.length) {
    note += ` The observed history of ${listing(restricted.map(varLabel))} is left out, under their publishers' terms.`;
  }
  const details = el("details", { class: "data-section", id: `data-section-${page}` },
    el("summary", {}, "Data"),
    el("div", { class: "data-actions" },
       downloadButton(`dl_${page}`, () => downloadCsv(pageCsvName(page), pageTable(page))),
       el("span", { class: "data-note" }, note)),
    el("div", { class: "data-table" }, el("div", { id: `tbl_${page}`, class: "shiny-html-output shiny-table-output" })));
  details.addEventListener("toggle", () => refreshTable(page));
  return details;
}

// The numbers behind a page (its Data section and CSV).
function pageTable(page) {
  return page === "irf" ? irfTable(state.irf.vars) : tabData(pageVariables(page));
}

function pageCsvName(page) {
  const last = CORE.meta.last_observed;
  return page === "irf" ? `bvar_irf_${state.irf.shock}_${last}.csv` : `bvar_${page}_data_${last}.csv`;
}

// Draw a page's table if its Data section is open.
function refreshTable(page) {
  const details = document.getElementById(`data-section-${page}`);
  if (!details || !details.open) return;
  document.getElementById(`tbl_${page}`).replaceChildren(tableElement(pageTable(page)));
}

// The pop-up with one variable's numbers (data_modal()), a Bootstrap
// modal as shiny::modalDialog(size = "l", easyClose = TRUE) draws it.
function openDataModal(variable, kind) {
  const table = kind === "irf" ? irfTable([variable]) : tabData([variable]);
  const restricted = kind !== "irf" && CORE.settings.no_history_download.includes(variable);
  const filename = kind === "irf"
    ? `bvar_irf_${variable}_${CORE.meta.last_observed}.csv`
    : `bvar_forecast_${variable}_${CORE.meta.last_observed}.csv`;
  document.getElementById("shiny-modal")?.remove();
  const modal = el("div", { class: "modal fade", id: "shiny-modal", tabindex: "-1" },
    el("div", { class: "modal-dialog modal-lg" },
       el("div", { class: "modal-content" },
          el("div", { class: "modal-header" }, el("h4", { class: "modal-title" }, varLabel(variable))),
          el("div", { class: "modal-body" },
             el("div", { class: "data-actions" },
                downloadButton("modal_dl", () => downloadCsv(filename, table)),
                el("span", { class: "data-note" }, "Numbers as the chart shows them, in the current settings.",
                   restricted ? " The observed history is left out, under its publisher's terms." : null)),
             el("div", { class: "data-table" }, tableElement(table))),
          el("div", { class: "modal-footer" },
             el("button", { type: "button", class: "btn btn-default", "data-dismiss": "modal", "data-bs-dismiss": "modal" },
                "Close")))));
  document.body.append(modal);
  modal.addEventListener("hidden.bs.modal", () => modal.remove());
  new bootstrap.Modal(modal).show();
}

// ---- Impulse responses (app/R/irf_data.R, fan_chart.R irf_chart()) ----
//
// How the model's forecast of each variable moves after a one-time shock
// to one variable, compared with no shock: quantiles over posterior draws
// of the response to a unit shock (1 log point for a log series, 1
// percentage point for a rate). Responses are linear in the shock, so a
// shock of size s scales every quantile by s, and a negative size swaps
// the lower and upper bands. The bands cover parameter uncertainty only.

const isLogSeries = (v) => variableRow(v).model_transform === "log";

// Default responses of a shock (irf_default_responses()): a chosen set,
// or the shock's own variable followed by the general set; at most six.
function irfDefaultResponses(shockId) {
  let set = CORE.settings.irf_default_responses[shockId];
  if (!set) set = [shockId.replace(/^[a-z]+_/, ""), ...CORE.settings.irf_general_responses];
  const vars = CORE.variables.map((v) => v.variable);
  return [...new Set(set)].filter((v) => vars.includes(v)).slice(0, 6);
}

// Unit of a response (irf_unit()).
function irfUnit(variable, transform) {
  if (isLogSeries(variable)) return transform === "yoy" ? "pp, 4-quarter change" : "%";
  return variableRow(variable).unit === "index" ? "index points" : "pp";
}

// Size unit of a shock (irf_shock_unit()).
function irfShockUnit(shockId) {
  return isLogSeries(shockId.replace(/^[a-z]+_/, "")) ? "log points (about %)" : "pp";
}

// One response, ready to plot (irf_series()): rows {horizon, q05..q95},
// in the transform used (a rate has no 4-quarter change: shown in levels).
function irfSeries(shockId, variable, transform, horizon, size) {
  const byTransform = IRF.responses[shockId][variable];
  if (!byTransform[transform]) transform = "level";
  const x = byTransform[transform];
  const q = ["q05", "q16", "q50", "q84", "q95"];
  const rows = [];
  for (let h = 0; h <= horizon; h++) {
    const r = { horizon: h };
    for (const k of q) r[k] = x[k][h] * size;
    // A negative shock turns the bands upside down: lower and upper swap.
    if (size < 0) { [r.q05, r.q95] = [r.q95, r.q05]; [r.q16, r.q84] = [r.q84, r.q16]; }
    rows.push(r);
  }
  return { rows, transform };
}

// The numbers of the IRF charts, as one table (irf_data()).
function irfTable(variables) {
  const s = state.irf;
  const rows = [];
  for (const v of variables) {
    const x = irfSeries(s.shock, v, s.units, s.horizon, s.size);
    for (const r of x.rows) {
      rows.push({ Variable: varLabel(v), Unit: irfUnit(v, x.transform), "Quarters after shock": r.horizon,
                  Median: r.q50, "68% low": r.q16, "68% high": r.q84, "90% low": r.q05, "90% high": r.q95 });
    }
  }
  const columns = variables.length
    ? ["Variable", "Unit", "Quarters after shock", "Median", ...BAND_COLUMNS] : ["Variable"];
  return { columns, rows };
}

// Impulse-response chart: median with shaded 68% and 90% bands (irf_chart()).
function irfChart(rows, colour = "#0055A4") {
  const traces = [
    ribbon(rows.map((r) => ({ ...r, date: r.horizon })), "q05", "q95", hexAlpha(colour, 0.14)),
    ribbon(rows.map((r) => ({ ...r, date: r.horizon })), "q16", "q84", hexAlpha(colour, 0.28)),
    line(rows.map((r) => r.horizon), rows.map((r) => r.q50), { color: colour, width: 2 },
         rows.map((r) => `Quarter ${r.horizon}: ${fixed3(r.q50)} (68%: ${fixed3(r.q16)} to ${fixed3(r.q84)})`)),
  ];
  const layout = {
    shapes: [{ type: "line", xref: "paper", yref: "y", x0: 0, x1: 1, y0: 0, y1: 0, line: { color: "#6F7680", width: 1 } }],
    hovermode: "closest", showlegend: false, margin: { l: 44, r: 10, t: 8, b: 36 },
    font: { family: "Public Sans, sans-serif", size: 11, color: "#1B1F24" },
    xaxis: { ...AXIS_STYLE, title: { text: "quarters after the shock", font: { size: 11 } }, dtick: 4 },
    yaxis: { ...AXIS_STYLE, title: { text: "" } },
    paper_bgcolor: "#FFFFFF", plot_bgcolor: "#FFFFFF",
  };
  return { traces, layout };
}

// Load data/irf.json the first time the IRF tab opens.
async function loadIrf() {
  if (IRF) return;
  const res = await fetch("data/irf.json");
  if (!res.ok) throw new Error(`data/irf.json: HTTP ${res.status}`);
  IRF = await res.json();
  const first = IRF.shocks[0].id;                               // the default shock is listed first
  state.irf = { shock: first, units: "level", horizon: Math.max(...CORE.settings.irf_horizons), size: 1,
                vars: irfDefaultResponses(first) };
  renderIrfHeader();
  buildIrfPanel();
}

// The shock size as the app reads it: a finite non-zero number, else 1.
function readSize(value) {
  const x = Number(value);
  return value !== "" && Number.isFinite(x) && x !== 0 ? x : 1;
}

// The response settings header (settings_header(), the IRF part).
function renderIrfHeader() {
  const s = state.irf;
  // The shock drop-down, grouped by identification (selectInput(selectize = TRUE)).
  const shockSelect = el("select", { id: "irf_shock" });
  const shock = el("div", { class: "setting setting-shock" },
    el("div", { class: "form-group shiny-input-container", style: "width:300px;" },
       el("label", { class: "control-label", id: "irf_shock-label", for: "irf_shock-selectized" }, "Shock"),
       el("div", {}, shockSelect)));
  const sizeInput = el("input", { id: "irf_size", type: "number", class: "shiny-input-number form-control", value: "1",
                                  step: "0.25",
                                  oninput: (e) => { state.irf.size = readSize(e.target.value); renderIrf(false); } });
  const sizeUnit = el("div", { id: "irf_size_unit", class: "shiny-html-output size-unit" }, el("span", {}, irfShockUnit(s.shock)));
  document.getElementById("header-irf").replaceChildren(
    shock,
    el("div", { class: "setting" }, radioGroup("irf_units", "Response of log series",
      [["level", "Level (%)"], ["yoy", "4-quarter change (pp)"]], s.units,
      (v) => { state.irf.units = v; renderIrf(true); })),
    el("div", { class: "setting setting-history" }, radioGroup("irf_horizon", "Quarters shown",
      CORE.settings.irf_horizons.map((h) => [String(h), String(h)]), String(s.horizon),
      (v) => { state.irf.horizon = Number(v); renderIrf(false); }, true)),
    el("div", { class: "setting setting-size" },
       el("div", { class: "form-group shiny-input-container", style: "width:90px;" },
          el("label", { class: "control-label", id: "irf_size-label", for: "irf_size" }, "Shock size"), sizeInput),
       sizeUnit),
  );
  const kinds = [...new Set(IRF.shocks.map((x) => x.kind))];
  const title = (k) => k.charAt(0).toUpperCase() + k.slice(1);
  $(shockSelect).selectize({
    plugins: ["selectize-plugin-a11y"], valueField: "id", labelField: "name", searchField: ["name"],
    optgroupField: "group", optgroups: kinds.map((k) => ({ value: title(k), label: title(k) })),
    options: IRF.shocks.map((x) => ({ id: x.id, name: x.name, group: title(x.kind) })), items: [s.shock],
    lockOptgroupOrder: true,
    onChange: (value) => {
      if (!value) return;
      state.irf.shock = value;
      sizeUnit.replaceChildren(el("span", {}, irfShockUnit(value)));
      setIrfVars(irfDefaultResponses(value));                    // a new shock shows its default responses
    },
  });
}

// The IRF tab's content (irf_panel()).
function buildIrfPanel() {
  const vars = CORE.variables.map((v) => v.variable);
  const select = el("select", { id: "irf_vars", multiple: true });
  document.getElementById("panel-irf").replaceChildren(
    el("div", { class: "chart-picker" },
       el("div", { class: "form-group shiny-input-container", style: "width: 100%;" },
          el("label", { class: "control-label", id: "irf_vars-label", for: "irf_vars-selectized" }, "Responses"),
          el("div", {}, select)),
       el("a", { href: "#", class: "action-button", id: "irf_reset",
                 onclick: (e) => { e.preventDefault(); setIrfVars(irfDefaultResponses(state.irf.shock)); } },
          "Reset to the shock's default responses")),
    el("div", { id: "grid_irf", class: "shiny-html-output" }),
    el("p", { class: "irf-note" },
       "Each chart shows how the model's forecast of a variable moves after a one-time shock, compared with ",
       "no shock: the median and 68% and 90% credible bands over parameter uncertainty. A generalized shock ",
       "orders the shocked variable first, so the other variables move on impact as they typically do with it; ",
       "it describes co-movement, not a structural cause. A recursive shock orders the variables in blocks: ",
       "those ordered before the shocked variable cannot respond on impact, those after it can. Responses ",
       "therefore depend on the identification. A unit shock is 1 log point (about 1%) for a log series and ",
       "1 percentage point for a rate; responses are linear, so the shock size scales them exactly."),
    dataSection("irf"));
  $(select).selectize({
    plugins: ["remove_button", "selectize-plugin-a11y"], maxItems: CORE.settings.overview_max,
    valueField: "value", labelField: "label", searchField: ["label"],
    options: vars.map((v) => ({ value: v, label: varLabel(v) })), items: state.irf.vars,
    onChange: (value) => { state.irf.vars = value || []; renderIrf(true); },
  });
}

function setIrfVars(vars) {
  document.getElementById("irf_vars").selectize.setValue(vars, true);
  state.irf.vars = [...vars];
  renderIrf(true);
}

// Draw the IRF tab; with rebuild, also its cards (their units depend on
// the response settings).
function renderIrf(rebuild = true) {
  if (!IRF || state.page !== "irf") return;
  const s = state.irf;
  if (rebuild) {
    document.getElementById("grid_irf").replaceChildren(el("div", { class: "chart-grid" }, s.vars.map((v) =>
      chartCard("irf", v, irfUnit(v, isLogSeries(v) ? s.units : "level")))));
  }
  for (const v of s.vars) {
    const { traces, layout } = irfChart(irfSeries(s.shock, v, s.units, s.horizon, s.size).rows);
    Plotly.react(`c_irf_${v}`, traces, layout, PLOT_CONFIG);
  }
  refreshTable("irf");
}

// ---- About (ui.R about_panel()) ----

function aboutPanel() {
  const m = CORE.meta;
  const shownAs = { growth: "growth rate or level (Growth setting)", natural: "level, in natural units", level: "level" };
  const nSeries = CORE.variables.filter((v) => v.kind === "series").length;
  const w = m.covid_window;
  const windowText = w.length === 2
    ? ` The equations of ${w[0]}–${w[1]} are left out of the estimation because of the pandemic's extreme movements.`
    : "";
  const table = el("table", { class: "about-table" },
    el("thead", {}, el("tr", {}, ["Variable", "Description", "Source", "In the model", "Shown as", "Notes"]
      .map((h) => el("th", {}, h)))),
    el("tbody", {}, CORE.variables.map((v) => el("tr", {},
      el("td", {}, v.label), el("td", {}, v.description), el("td", {}, v.source),
      el("td", {}, v.model_transform === "log" ? "100 × ln(level)" : "as published"),
      el("td", {}, shownAs[v.display]), el("td", {}, v.note)))));
  return el("div", { class: "about" },
    el("h2", {}, "The model"),
    el("p", {}, `A Bayesian vector autoregression (BVAR) of ${nSeries} quarterly US series with ${m.lags} lags, following `,
       "Crump, Eusepi, Giannone, Qian and Sbordone (2021), “A Large Bayesian VAR of the ",
       "United States Economy” (Federal Reserve Bank of New York Staff Report 976). The series ",
       "enter in levels. The prior combines a Minnesota prior with a sum-of-coefficients prior; how ",
       "tight they are is estimated from the data (Giannone, Lenza and Primiceri 2015)."),
    el("p", {}, `Estimated on ${m.first_quarter}–${m.last_observed}.${windowText} Forecasts start from the latest data and the site shows ${m.horizon_exported} `,
       "quarters. Beyond two to three years a model in levels mostly extrapolates trends in its ",
       "sample, so longer horizons are not shown."),
    el("p", {}, `Each forecast is simulated from ${m.n_draws.toLocaleString("en-US")} posterior draws of the parameters, with ${m.paths_per_draw} paths `,
       "per draw. The bands are 68% and 90% Bayesian credible intervals: they cover parameter and ",
       "shock uncertainty under the model, not the uncertainty about the model itself."),
    el("h2", {}, "Scenarios"),
    el("p", {}, "A scenario is a conditional forecast: what the model expects given the scenario's ",
       "assumptions. Some scenarios fix the path of a few variables exactly (marked on the charts ",
       "with hollow circles); others reweight the forecast so that its median meets a target, such ",
       "as a survey's inflation forecast (Crump et al. 2021, Section 4.2; targets marked with ",
       "hollow diamonds). Scenarios are listed, renamed and recoloured on the Settings tab. Scenarios are not predictions. The difference view shows one scenario's ",
       "median minus the baseline's, a difference between two conditional forecasts, not a causal ",
       "effect. Hover over a scenario's name for its assumptions."),
    el("h2", {}, "Data"),
    el("p", {}, `Quarterly, ${m.first_quarter}–${m.last_observed}. FRED-QD vintage ${m.vintage}; series from the FRED API as of ${m.api_realtime_date}. `,
       "Monthly series are quarterly averages."),
    el("div", { class: "about-table-wrap" }, table),
    el("h2", {}, "Sources and credits"),
    el("p", {}, "Sources: U.S. Bureau of Economic Analysis, U.S. Bureau of Labor Statistics, U.S. Census ",
       "Bureau, Board of Governors of the Federal Reserve System, via FRED and FRED-QD (Federal ",
       "Reserve Bank of St. Louis); Moody's (corporate yields); University of Michigan Surveys of ",
       "Consumers (sentiment); S&P Dow Jones Indices (S&P 500); World Bank Commodity Price Data ",
       "(non-energy index, CC BY 4.0; used as a proxy for the S&P GSCI non-energy index); ",
       "realized volatility computed from daily S&P 500 closes (Yahoo Finance)."),
    el("p", {}, "This product uses the FRED® API but is not endorsed or certified by the Federal ",
       "Reserve Bank of St. Louis."));
}

// ---- Settings tab (ui.R settings_panel(), server.R) ----

const SCENARIO_TYPES = { baseline: "unconditional forecast", tilted: "tilted to targets (soft conditions)",
                         hard: "conditional forecast (fixed paths)" };

// One row of the scenario table (scenario_settings_row()).
function scenarioSettingsRow(s) {
  const id = s.id;
  const label = s.label;                                         // the scenario's own short label
  const listBox = el("input", { id: `list_${id}`, type: "checkbox", class: "shiny-input-checkbox",
                                checked: state.listed.includes(id),
                                onchange: (e) => setListed(id, e.target.checked) });
  const colour = el("input", { type: "color", id: `colour_${id}`, value: state.colours[id], class: "colour-input",
                               oninput: (e) => { state.colours[id] = e.target.value; scenariosChanged(); } });
  // A blank name falls back to the scenario's own short label.
  const name = el("input", { id: `name_${id}`, type: "text", class: "shiny-input-text form-control",
                             value: state.names[id],
                             onchange: (e) => {
                               const v = e.target.value.trim();
                               state.names[id] = v ? v : label;
                               scenariosChanged();
                             } });
  return el("tr", {},
    el("td", { class: "col-narrow" },
       el("div", { class: "form-group shiny-input-container" },
          el("div", { class: "checkbox" }, el("label", {}, listBox, el("span", {}, srOnly(`List ${label}`)))))),
    el("td", { class: "col-narrow" },
       el("label", { class: "visually-hidden", for: `colour_${id}` }, `Colour of ${label}`), colour),
    el("td", { class: "col-name" },
       el("div", { class: "form-group shiny-input-container" },
          el("label", { class: "control-label", id: `name_${id}-label`, for: `name_${id}` }, srOnly(`Name of ${label}`)),
          name),
       el("div", { class: "scenario-type" }, SCENARIO_TYPES[s.type] ?? s.type)),
    el("td", { class: "scenario-description" }, s.description));
}

// The scenario table (output$scenario_table), rebuilt on reset only.
function renderScenarioTable() {
  document.getElementById("scenario_table").replaceChildren(
    el("table", { class: "scenario-table" },
       el("thead", {}, el("tr", {},
          el("th", { class: "col-narrow" }, "Listed"), el("th", { class: "col-narrow" }, "Colour"),
          el("th", { class: "col-name" }, "Name"), el("th", {}, "Assumptions"))),
       el("tbody", {}, CORE.scenarios.map(scenarioSettingsRow))));
}

// List or unlist a scenario, keeping the catalog's order.
function setListed(id, on) {
  const listed = on ? [...state.listed, id] : state.listed.filter((x) => x !== id);
  state.listed = CORE.scenarios.map((s) => s.id).filter((x) => listed.includes(x));
  scenariosChanged();
}

// The header's scenario list after the listed scenarios, names or colours
// change (output$scenario_controls): the reader's ticks are kept, a newly
// listed scenario starts shown, and the baseline falls back to the first
// listed scenario if it was unlisted.
function scenariosChanged() {
  const ids = state.listed;
  state.show = [...state.show.filter((x) => ids.includes(x)), ...ids.filter((x) => !state.known.includes(x))];
  state.bands = state.bands.filter((x) => ids.includes(x));
  if (!ids.includes(state.baseline)) state.baseline = ids[0] ?? null;
  state.known = [...ids];
  renderForecastHeader();
  renderCharts();
}

// Restore the default list, names and colours (input$reset_scenarios).
function resetScenarios() {
  const cat = CORE.scenarios;
  state.listed = cat.filter((s) => s.default).map((s) => s.id);
  for (const s of cat) { state.colours[s.id] = s.colour; state.names[s.id] = s.label; }
  renderScenarioTable();
  scenariosChanged();
}

// The Settings tab (settings_panel()).
function settingsPanel() {
  return el("div", { class: "settings-page" },
    el("h2", {}, "Scenarios"),
    el("p", {}, "Choose the scenarios listed in the header, their names and their colours. ",
       "The choices last for this visit."),
    el("div", { id: "scenario_table", class: "shiny-html-output" }),
    el("a", { href: "#", id: "reset_scenarios", class: "action-button",
              onclick: (e) => { e.preventDefault(); resetScenarios(); } }, "Reset scenarios to default"));
}

// ---- Layout pieces (ui.R) ----

// A short line in a scenario's colour (swatch()).
function swatch(colour) {
  const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  svg.setAttribute("width", "18"); svg.setAttribute("height", "10");
  svg.setAttribute("aria-hidden", "true"); svg.setAttribute("class", "swatch");
  const l = document.createElementNS("http://www.w3.org/2000/svg", "line");
  for (const [k, v] of Object.entries({ x1: 0, y1: 5, x2: 18, y2: 5, stroke: colour, "stroke-width": 2 })) l.setAttribute(k, v);
  svg.append(l);
  return svg;
}

const srOnly = (text) => el("span", { class: "visually-hidden" }, text);

// Shiny's radioButtons(): label, then one radio per choice.
function radioGroup(id, label, choices, selected, onChange, inline = false) {
  const opts = choices.map(([value, text]) => {
    const input = el("input", { type: "radio", name: id, value, checked: value === selected,
                                onchange: () => onChange(value) });
    const lab = el("label", inline ? { class: "radio-inline" } : {}, input, el("span", {}, text));
    return inline ? lab : el("div", { class: "radio" }, lab);
  });
  const cls = `form-group shiny-input-radiogroup shiny-input-container${inline ? " shiny-input-container-inline" : ""}`;
  return el("div", { id, class: cls, role: "radiogroup", "aria-labelledby": `${id}-label` },
            el("label", { class: "control-label", id: `${id}-label`, for: id }, label),
            el("div", { class: "shiny-options-group" }, opts));
}

// Shiny's checkboxGroupInput().
function checkboxGroup(id, label, values, names, selected, onChange) {
  const opts = values.map((value, k) => el("div", { class: "checkbox" },
    el("label", {}, el("input", { type: "checkbox", name: id, value, checked: selected.includes(value),
                                  onchange: (e) => onChange(value, e.target.checked) }),
       el("span", {}, names[k]))));
  return el("div", { id, class: "form-group shiny-input-checkboxgroup shiny-input-container", role: "group",
                     "aria-labelledby": `${id}-label` },
            el("label", { class: "control-label", id: `${id}-label`, for: id }, label),
            el("div", { class: "shiny-options-group" }, opts));
}

// The scenario list: Show, Bands and Baseline columns (scenario_controls()).
function scenarioControls() {
  const x = listedScenarios();
  const ids = x.map((s) => s.id);
  if (!ids.length) {
    return el("fieldset", { class: "setting" }, el("legend", {}, "Scenarios"),
              el("p", { class: "empty-note" }, "No scenario listed: choose some on the Settings tab."));
  }
  const toggle = (key) => (value, on) => {
    state[key] = on ? [...new Set([...state[key], value])] : state[key].filter((v) => v !== value);
    renderCharts();
  };
  const show = checkboxGroup("show", "Show", ids,
    x.map((s) => el("span", { class: "scenario-name", title: s.description }, swatch(s.colour), s.label)),
    state.show, toggle("show"));
  const bands = isDifferenceMode() ? null : checkboxGroup("bands", "Bands", ids,
    x.map((s) => srOnly(`Bands for ${s.label}`)), state.bands, toggle("bands"));
  const baseline = radioGroup("baseline", "Baseline", x.map((s) => [s.id, srOnly(`${s.label} as baseline`)]),
    state.baseline, (v) => { state.baseline = v; renderCharts(); });
  return el("fieldset", { class: "setting" }, el("legend", {}, "Scenarios"),
            el("div", { class: "scenario-grid" }, show, bands, baseline));
}

// Bands are not drawn in the difference view, so their column is hidden there.
const isDifferenceMode = () => state.mode === "difference";

// The forecast settings header (settings_header()).
function renderForecastHeader() {
  const header = document.getElementById("header-forecast");
  header.replaceChildren(
    el("div", { id: "scenario_controls", class: "shiny-html-output" }, scenarioControls()),
    el("div", { class: "setting" }, radioGroup("mode", "Show",
      [["forecast", "Forecasts"], ["difference", "Difference from baseline (medians)"]], state.mode,
      (v) => { state.mode = v; renderForecastHeader(); renderPage(); })),
    el("div", { class: "setting" }, radioGroup("growth", "Growth",
      [["yoy", "Year over year"], ["qoq_ann", "Quarterly, annualized"], ["natural", "Level"]], state.growth,
      (v) => { state.growth = v; renderPage(); })),
    el("div", { class: "setting setting-history" }, radioGroup("history", "History",
      CORE.settings.history.map((h) => [h.name, h.name]), state.history,
      (v) => { state.history = v; renderCharts(); }, true)),
  );
}

// One tab in the side panel (page_tab()).
function pageTab(value, label, cls) {
  return el("div", { class: `radio ${cls}` },
    el("label", {}, el("input", { type: "radio", name: "page", value, checked: value === state.page,
                                  onchange: () => { state.page = value; renderPage(); } }),
       el("span", {}, label)));
}

// The side panel's tabs (side_panel()).
function renderSidePanel() {
  document.getElementById("page-tabs").replaceChildren(
    el("div", { class: "nav-group" }, "Forecasts"),
    pageTab("overview", "Overview", "sub"),
    ...CORE.settings.topics.map((t) => pageTab(t.id, t.title, "sub")),
    pageTab("irf", "Impulse responses", "top"),
    el("div", { class: "nav-rule" }),
    pageTab("settings", "Settings", "top"),
    pageTab("about", "About the model and data", "top"),
  );
}

// A small table icon, for the chart's data button (table_icon()).
function tableIcon() {
  const span = el("span", {});
  span.innerHTML = '<svg width="16" height="14" viewBox="0 0 16 14" aria-hidden="true" focusable="false">' +
    '<rect x="0.75" y="0.75" width="14.5" height="12.5" rx="1" fill="none" stroke="currentColor" stroke-width="1.5"/>' +
    '<line x1="0.75" y1="4.75" x2="15.25" y2="4.75" stroke="currentColor" stroke-width="1.5"/>' +
    '<line x1="0.75" y1="9" x2="15.25" y2="9" stroke="currentColor" stroke-width="1"/>' +
    '<line x1="6" y1="4.75" x2="6" y2="13.25" stroke="currentColor" stroke-width="1"/></svg>';
  return span.firstChild;
}

// A card holding one chart (chart_card()).
function chartCard(page, variable, unit) {
  const label = varLabel(variable);
  return el("section", { class: "chart-card" },
    el("div", { class: "chart-head" },
       // The space matches the app's HTML, where a line break separates the two.
       el("h2", {}, label, " ", el("span", { class: "unit" }, unit)),
       el("button", { type: "button", class: "btn btn-default action-button chart-data-btn", id: `data_${page}_${variable}`,
                      title: `Data table: ${label}`, "aria-label": `Data table for ${label}`,
                      onclick: () => openDataModal(variable, page === "irf" ? "irf" : "forecast") }, tableIcon())),
    el("div", { class: "chart-plot", id: `c_${page}_${variable}` }));
}

// ---- Pages ----

// The variables a forecast page shows.
function pageVariables(page) {
  if (page === "overview") return state.overviewVars;
  const t = CORE.settings.topics.find((x) => x.id === page);
  return t ? t.variables : null;
}

// Build every tab's container once; renderPage() shows the current one.
function buildPanels() {
  const content = document.getElementById("content");
  const vars = CORE.variables.map((v) => v.variable);
  // Overview: chart selector (selectize, as the app's selectizeInput()) and charts.
  const select = el("select", { id: "overview_vars", multiple: true, "aria-labelledby": "overview_vars-label" });
  const overview = el("div", { class: "page-panel", id: "panel-overview" },
    el("div", { class: "chart-picker" },
       el("div", { class: "form-group shiny-input-container", style: "width: 100%;" },
          el("label", { class: "control-label", id: "overview_vars-label", for: "overview_vars-selectized" }, "Charts"),
          el("div", {}, select)),
       el("a", { href: "#", class: "action-button", id: "overview_reset",
                 onclick: (e) => { e.preventDefault(); setOverview(CORE.settings.overview_default); } },
          "Reset to default")),
    el("div", { id: "grid_overview", class: "shiny-html-output" }),
    dataSection("overview"));
  content.replaceChildren(overview);
  for (const t of CORE.settings.topics) {
    content.append(el("div", { class: "page-panel", id: `panel-${t.id}`, hidden: true },
                      el("div", { id: `grid_${t.id}`, class: "shiny-html-output" }), dataSection(t.id)));
  }
  // The IRF tab is filled when it first opens (loadIrf()).
  content.append(el("div", { class: "page-panel", id: "panel-irf", hidden: true },
                    el("p", { class: "placeholder" }, "Loading the impulse responses…")));
  content.append(el("div", { class: "page-panel", id: "panel-settings", hidden: true }, settingsPanel()));
  content.append(el("div", { class: "page-panel", id: "panel-about", hidden: true }, aboutPanel()));
  renderScenarioTable();
  $(select).selectize({
    plugins: ["remove_button", "selectize-plugin-a11y"], maxItems: CORE.settings.overview_max,
    valueField: "value", labelField: "label", searchField: ["label"],
    options: vars.map((v) => ({ value: v, label: varLabel(v) })), items: state.overviewVars,
    onChange: (value) => { state.overviewVars = value || []; if (state.page === "overview") renderPage(); },
  });
}

function setOverview(vars) {
  const sel = document.getElementById("overview_vars").selectize;
  sel.setValue(vars, true);
  state.overviewVars = [...vars];
  renderPage();
}

// Show the current tab, and on a forecast tab rebuild its cards and charts.
function renderPage() {
  for (const p of document.querySelectorAll(".page-panel")) p.hidden = p.id !== `panel-${state.page}`;
  document.getElementById("header-forecast").hidden = state.page === "irf";
  document.getElementById("header-irf").hidden = state.page !== "irf";
  if (state.page === "irf") {
    loadIrf().then(() => renderIrf(true)).catch((err) => {
      document.getElementById("panel-irf").replaceChildren(
        el("p", { class: "placeholder" }, `The impulse responses could not be loaded (${err.message}).`));
    });
    return;
  }
  const vars = pageVariables(state.page);
  if (!vars) return;
  const difference = state.mode === "difference";
  document.getElementById(`grid_${state.page}`).replaceChildren(
    el("div", { class: "chart-grid" },
       vars.map((v) => chartCard(state.page, v, unitLabel(v, displayTransform(v, state.growth), difference)))));
  renderCharts();
}

// Redraw the charts of the current tab (the cards stay).
function renderCharts() {
  const vars = pageVariables(state.page);
  if (!vars) return;
  for (const v of vars) {
    const { traces, layout } = makeChart(v);
    Plotly.react(`c_${state.page}_${v}`, traces, layout, PLOT_CONFIG);
  }
  refreshTable(state.page);
}

// Footer: data vintage, last observed quarter, run (app_ui()), and a link
// to the source code (static site only).
function renderFooter() {
  const m = CORE.meta;
  document.getElementById("footer").replaceChildren(
    el("span", {}, `Data through ${m.last_observed}`),
    el("span", {}, `FRED-QD vintage ${m.vintage}`),
    el("span", {}, `Run ${m.run_id}`),
    el("span", {}, "Credible bands, not confidence intervals; scenarios are conditional forecasts."),
    el("span", {}, "This product uses the FRED® API but is not endorsed or certified by the Federal Reserve Bank of St. Louis."),
    el("span", {}, "Code and methodology: ",
       el("a", { href: "https://github.com/zayat-alex/large-bvar-usa" }, "github.com/zayat-alex/large-bvar-usa")),
  );
}

// ---- Start ----

async function start() {
  const res = await fetch("data/core.json");
  if (!res.ok) throw new Error(`data/core.json: HTTP ${res.status}`);
  CORE = await res.json();
  // Initial choices, as the app's server: the default scenarios listed and
  // shown, bands and baseline on the first one.
  const cat = CORE.scenarios;
  state.listed = cat.filter((s) => s.default).map((s) => s.id);
  for (const s of cat) { state.colours[s.id] = s.colour; state.names[s.id] = s.label; }
  state.show = [...state.listed];
  state.bands = state.listed.slice(0, 1);
  state.baseline = state.listed[0] ?? null;
  state.known = [...state.listed];
  state.history = CORE.settings.history[0].name;
  state.overviewVars = [...CORE.settings.overview_default];
  renderSidePanel();
  renderForecastHeader();
  renderFooter();
  buildPanels();
  renderPage();
}

start().catch((err) => {
  document.getElementById("content").replaceChildren(
    el("p", { class: "placeholder" }, `The site could not load its data (${err.message}).`));
  throw err;
});
