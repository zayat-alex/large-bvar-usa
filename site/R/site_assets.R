# site_assets.R: the static website's third-party files, taken from the R
# packages the Shiny app uses, so the two sites look the same.
#
# Writes site/vendor/ (rebuilt on every run; nothing in it is edited by
# hand):
#   bootstrap/  Bootstrap 5 compiled with the app's theme by bslib, its
#               fonts and its JavaScript (the table pop-up);
#   shiny/      Shiny's stylesheet, compiled with the same theme;
#   bslib/      bslib's component styles (for example the blurred backdrop
#               of the table pop-up);
#   fontawesome/ the icon font of the download buttons (shiny::icon());
#   selectize/  the multi-select of the chart and response pickers, as
#               Shiny serves it, and jQuery, which it needs;
#   plotly/     plotly.js, the library the R plotly package draws with;
#   app.css     a copy of the app's stylesheet (app/www/app.css);
#   NOTICE.txt  the licences of these files (see site_write_notice()).
# The theme must equal app_ui()'s bs_theme() call (app/R/ui.R).
#
# Licences: Bootstrap, jQuery, Font Awesome and plotly.js carry their notices
# in their own files. The others (selectize.js, Apache 2.0; its
# accessibility plugin, Shiny's and bslib's stylesheets, MIT; the Open Sans
# and Source Code Pro fonts, SIL Open Font License 1.1) are minified or
# compiled without one, so NOTICE.txt gives their copyright lines and
# licence texts. The MIT and Apache texts are R's own copies
# (R.home("share")/licenses); the fonts' licence files, with their
# copyright lines, are in site/licenses/ (from the Google Fonts repository,
# github.com/google/fonts, ofl/opensans and ofl/sourcecodepro).


#' The app's Bootstrap theme (app/R/ui.R, app_ui()).
site_theme <- function() {
  bslib::bs_theme(version = 5, primary = "#0055A4", "font-family-sans-serif" = "'Public Sans', sans-serif")
}

#' Render the theme's dependencies to a folder, as Shiny serves them.
#'
#' A small page with the app's theme and one selectize input is saved with
#' htmltools, which writes every dependency (compiled with the theme) into
#' lib/.
#'
#' @param dir A scratch folder.
#' @return Path of the lib/ folder.
site_render_dependencies <- function(dir) {
  page <- bslib::page(theme = site_theme(),
                      shiny::selectizeInput("x", "x", choices = c("a", "b"), multiple = TRUE,
                                            options = list(plugins = list("remove_button"))),
                      shiny::icon("download"))                       # Font Awesome (download buttons)
  # Shiny's own stylesheet (input widths, label layout), compiled with the
  # theme. Shiny adds it to a running app itself; it has no exported
  # function for it, so the internal one is used (versions pinned by renv).
  page <- htmltools::attachDependencies(page, shiny:::shinyDependencyCSS(site_theme()), append = TRUE)
  htmltools::save_html(page, file.path(dir, "page.html"), libdir = "lib")
  file.path(dir, "lib")
}

#' Copy files, stopping if one is missing.
site_copy <- function(from, to_dir) {
  if (!all(file.exists(from))) stop("missing: ", paste(from[!file.exists(from)], collapse = ", "))
  dir.create(to_dir, showWarnings = FALSE, recursive = TRUE)
  ok <- file.copy(from, to_dir, overwrite = TRUE, recursive = TRUE)
  if (!all(ok)) stop("could not copy: ", paste(from[!ok], collapse = ", "))
}

#' Write site/vendor/.
#'
#' @param vendor The vendor folder (replaced).
#' @param app_css Path of the app's stylesheet.
#' @return The versions written, invisibly (list).
write_site_vendor <- function(vendor, app_css) {
  unlink(vendor, recursive = TRUE)
  scratch <- tempfile("site-deps-")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  lib <- site_render_dependencies(scratch)
  pick <- function(prefix) {
    d <- list.files(lib, pattern = paste0("^", prefix, "-"), full.names = TRUE)
    if (length(d) != 1) stop("expected one ", prefix, " folder in the rendered dependencies")
    d
  }
  bs <- pick("bootstrap")
  site_copy(file.path(bs, c("bootstrap.min.css", "font.css", "fonts", "bootstrap.bundle.min.js")),
            file.path(vendor, "bootstrap"))
  sz <- pick("selectize")
  site_copy(file.path(sz, c("selectize.css", "selectize.min.js", "selectize-plugin-a11y.min.js")),
            file.path(vendor, "selectize"))
  sh <- pick("shiny-sass")
  site_copy(file.path(sh, "shiny-sass.css"), file.path(vendor, "shiny"))
  bc <- pick("bslib-component-css")
  site_copy(file.path(bc, "bslib-component-css.min.css"), file.path(vendor, "bslib"))
  fa <- pick("font-awesome")
  site_copy(file.path(fa, c("css", "webfonts")), file.path(vendor, "fontawesome"))
  jq <- pick("jquery")
  site_copy(file.path(jq, "jquery.min.js"), file.path(vendor, "selectize"))
  pl <- system.file("htmlwidgets", "lib", "plotlyjs", package = "plotly")
  site_copy(file.path(pl, c("plotly-latest.min.js", "LICENSE")), file.path(vendor, "plotly"))
  site_copy(app_css, vendor)
  site_write_notice(vendor, file.path("site", "licenses"))
  invisible(list(bootstrap = basename(bs), shiny = basename(sh), bslib = basename(bc), fontawesome = basename(fa),
                 selectize = basename(sz), jquery = basename(jq),
                 plotly_r = as.character(utils::packageVersion("plotly"))))
}

#' Copyright line from an R package's MIT LICENSE file ("YEAR: ...",
#' "COPYRIGHT HOLDER: ...").
site_package_copyright <- function(pkg) {
  lines <- readLines(system.file("LICENSE", package = pkg, mustWork = TRUE), warn = FALSE)
  field <- function(name) trimws(sub(paste0("^", name, ":"), "", grep(paste0("^", name, ":"), lines, value = TRUE)))
  year <- field("YEAR")
  holder <- field("COPYRIGHT HOLDER")
  if (length(year) != 1 || length(holder) != 1) stop("unexpected LICENSE file in package ", pkg)
  sprintf("Copyright (c) %s %s", year, holder)
}

#' Write vendor/NOTICE.txt: the components without a notice of their own,
#' their copyright lines, and the full licence texts.
#'
#' Versions and copyright lines are read from the installed packages, so
#' the notice follows the versions renv pins. The fonts listed are the ones
#' the compiled theme actually bundles (the families in
#' bootstrap/font.css); the build stops if one has no licence file in
#' `licenses_dir`.
#'
#' @param vendor The vendor folder (already written).
#' @param licenses_dir Folder with one "OFL-<family>.txt" per font.
#' @return Path of the notice, invisibly.
site_write_notice <- function(vendor, licenses_dir) {
  # selectize.js: version and copyright from the header of the unminified
  # file that Shiny ships (the minified copy served here has no header).
  sz_head <- readLines(system.file("www", "shared", "selectize", "js", "selectize.js",
                                   package = "shiny", mustWork = TRUE), n = 20, warn = FALSE)
  sz_version <- sub(".*Selectize [(](v[0-9.]+)[)].*", "\\1", grep("Selectize (v", sz_head, value = TRUE, fixed = TRUE))
  sz_copyright <- trimws(sub("^ [*]", "", grep("Copyright (c)", sz_head, value = TRUE, fixed = TRUE)))
  if (length(sz_version) != 1 || length(sz_copyright) == 0) stop("unexpected selectize.js header")

  css <- readLines(file.path(vendor, "bootstrap", "font.css"), warn = FALSE)
  families <- sort(unique(sub(".*font-family: *'([^']+)'.*", "\\1", grep("font-family:", css, value = TRUE))))
  font_files <- file.path(licenses_dir, paste0("OFL-", families, ".txt"))
  if (!all(file.exists(font_files))) stop("no licence file for font(s): ", paste(families[!file.exists(font_files)], collapse = ", "))

  licence_text <- function(name) readLines(file.path(R.home("share"), "licenses", name), warn = FALSE)
  mit <- licence_text("MIT")
  mit <- mit[grep("^Permission is hereby granted", mit):length(mit)]   # the text after R's template note
  rule <- strrep("=", 72)

  out <- c(
    "Third-party files in this folder",
    "",
    "Bootstrap (bootstrap/), jQuery (selectize/jquery.min.js), Font Awesome",
    "(fontawesome/) and plotly.js (plotly/) carry their copyright and licence",
    "notices in their own files. The files below do not, and are listed here.",
    "",
    sprintf("- selectize.js %s (selectize/selectize.min.js, selectize/selectize.css),", sz_version),
    "  Apache License 2.0 (text below).",
    paste0("  ", sz_copyright),
    "- selectize-plugin-a11y (selectize/selectize-plugin-a11y.min.js), MIT License.",
    "  Copyright 2018-present Salmen Bejaoui",
    sprintf("- Shiny %s stylesheet (shiny/shiny-sass.css), MIT License.", utils::packageVersion("shiny")),
    paste0("  ", site_package_copyright("shiny")),
    sprintf("- bslib %s component styles (bslib/bslib-component-css.min.css), MIT License.", utils::packageVersion("bslib")),
    paste0("  ", site_package_copyright("bslib")),
    sprintf("- The %s fonts (bootstrap/fonts/), SIL Open Font License 1.1;", paste(families, collapse = " and ")),
    "  copyright lines and licence below.",
    "",
    "MIT License text (applies to the MIT-licensed files above, each with its",
    "copyright line):",
    "", mit, "", rule,
    "Apache License 2.0 (selectize.js)", rule, licence_text("Apache-2.0"))
  for (k in seq_along(families)) {
    out <- c(out, "", rule, paste0(families[k], " (SIL Open Font License 1.1)"), rule,
             readLines(font_files[k], warn = FALSE))
  }
  path <- file.path(vendor, "NOTICE.txt")
  writeLines(out, path, useBytes = TRUE)
  invisible(path)
}
