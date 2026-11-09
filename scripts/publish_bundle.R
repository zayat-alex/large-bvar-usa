# publish_bundle.R: copy the published run's bundle into the app folder.
#
# Run from the project root (so .Rprofile activates renv), after
# scripts/export.R:
#   Rscript scripts/publish_bundle.R
#
# The app reads only app/bundle/bundle.rds (and meta.json beside it, for a
# reader of the folder). The copy replaces the previous one; app/bundle/ is
# not tracked in git, since it is rebuilt from the pipeline.

library(targets)

paths <- tar_read(export_bundle)                 # outputs/<run_id>/bundle/bundle.rds and meta.json
dest <- file.path("app", "bundle")
dir.create(dest, showWarnings = FALSE, recursive = TRUE)
ok <- file.copy(paths, dest, overwrite = TRUE)
if (!all(ok)) stop("could not copy: ", paste(paths[!ok], collapse = ", "))
cat("Copied", paste(basename(paths), collapse = " and "), "from", dirname(paths[1]), "to", dest, "\n")
