# export.R: run stage 5 of the pipeline (_targets.R), the export, and print a
# summary of the bundle.
#
# Run from the project root (so .Rprofile activates renv):
#   Rscript scripts/export.R
#
# Reads the engine outputs of the published run (config/run.yml, export$run)
# and writes the tables the website reads (R/export.R):
#   <output_root>/<run_id>/bundle/bundle.rds   the tables (one R list)
#   <output_root>/<run_id>/bundle/meta.json    run, data, settings, code
# Run scripts/run_engine.R first if the engine outputs are missing or out of
# date; this script builds only what the export depends on, which is cheap
# when the engine is up to date (the export itself takes under a minute).

library(targets)

tar_make(names = "export_bundle")

paths <- tar_read(export_bundle)
bundle <- readRDS(paths[1])
meta <- bundle$meta
cat("\nBundle:", dirname(paths[1]), "\n")
cat("  run_id:", meta$run_id, "| data through", meta$data$last_observed, "| forecast",
    paste(meta$forecast$quarters, collapse = "-"), "\n")
cat("  size:", round(file.size(paths[1]) / 1e6, 2), "MB\n")
for (tab in setdiff(names(bundle), "meta")) cat(sprintf("  %-13s %7d rows\n", tab, NROW(bundle[[tab]])))

# A few numbers to read at a glance in the fourth quarter of the second
# forecast year (2027Q4 for data ending in 2026Q2): medians over paths of
# 4-quarter growth of real GDP and of core PCE prices (exact percent
# changes), and of the unemployment rate.
f <- bundle$forecasts
q_end <- grep("Q4$", unique(f$quarter), value = TRUE)[2]
med <- function(obj, v, tr) f$q50[f$object == obj & f$variable == v & f$transform == tr & f$quarter == q_end]
cat("\nMedians in", q_end, "\n")
print(data.frame(object = bundle$objects$id, type = bundle$objects$type,
                 rgdp_yoy = vapply(bundle$objects$id, med, 0, "rgdp", "yoy"),
                 corepce_yoy = vapply(bundle$objects$id, med, 0, "corepce", "yoy"),
                 urate = vapply(bundle$objects$id, med, 0, "urate", "level"),
                 ess_share = bundle$objects$ess_share, row.names = NULL), digits = 3)
