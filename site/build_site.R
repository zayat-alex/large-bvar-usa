# build_site.R: write the static website's data and third-party files.
#
# Run from the project root (so .Rprofile activates renv), after
# scripts/publish_bundle.R:
#   Rscript site/build_site.R
#
# Reads app/bundle/bundle.rds, the same file the Shiny app reads, and the
# app's settings (app/R/), and writes site/data/core.json and irf.json
# (site/R/site_data.R) and site/vendor/ (site/R/site_assets.R). The static
# site (site/index.html) then shows the same numbers as the app without a
# server. To view it locally, serve the site/ folder over http (browsers
# do not load the data files from a file:// address), for example:
#   Rscript -e "httpuv::runStaticServer('site')"

app <- new.env()
for (f in list.files(file.path("app", "R"), pattern = "[.]R$", full.names = TRUE)) sys.source(f, app)
source(file.path("site", "R", "site_data.R"))
source(file.path("site", "R", "site_assets.R"))

bundle <- app$load_bundle(file.path("app", "bundle", "bundle.rds"))
paths <- write_site_data(bundle, app, file.path("site", "data"))
for (p in paths) cat(sprintf("%s: %.0f KB\n", p, file.size(p) / 1024))

versions <- write_site_vendor(file.path("site", "vendor"), file.path("app", "www", "app.css"))
cat("site/vendor:", paste(unlist(versions), collapse = ", "), "\n")
