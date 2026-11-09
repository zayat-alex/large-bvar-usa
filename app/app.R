# app.R: the website, a read-only Shiny app.
#
# It reads one file, the export bundle (bundle/bundle.rds, copied here by
# scripts/publish_bundle.R from the project root), and runs no model code,
# so it can be moved to any host that runs Shiny. Shiny sources the files
# in R/ before this one: settings.R (labels, scenarios, tabs), bundle.R
# (loading), chart_data.R and fan_chart.R (charts), ui.R and server.R.
#
# Run from the project root:  shiny::runApp("app")
# The bundle's location can be changed with the environment variable
# BVAR_BUNDLE (used by the tests).

bundle <- load_bundle(Sys.getenv("BVAR_BUNDLE", "bundle/bundle.rds"))

# Every scenario in the bundle can be listed (Settings tab); the ones in
# APP_DEFAULT_SCENARIOS are listed when the site opens.
shiny::shinyApp(ui = app_ui(bundle), server = app_server(bundle))
