# build_data.R: run the data stage of the pipeline (_targets.R) and print what it built.
#
# Run from the project root (so .Rprofile activates renv):
#   Rscript scripts/build_data.R
#
# Reads config/run.yml, config/series.csv and the raw files they name; never
# goes online. Steps whose inputs have not changed since the last run are
# skipped (cache in _targets/). Writes to cfg$data$output_dir, e.g.
# outputs/data_2026m09/, replacing files of a previous build of that vintage:
#   levels.csv, model_data.csv   the two data tables (docs/data.md, section 1)
#   data.rds                     both tables plus the data and sample config used and session info
#   validation.rds               results of every check (R/data_validate.R)
#   checks/*.png                 plots of the checks
#   exploratory_series.pdf       line plots of all 31 series, six per page
#   data_validation.html         the validation report

# Only the data targets (stage 1) are built; the estimation is run by
# scripts/estimate.R.
# (tar_make() reads `names` in its own environment, so the names are written
# out here rather than taken from a variable.)
targets::tar_make(names = c(data_files, check_plots, validation_file, exploratory_pdf,
                            validation_report))

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
levels <- targets::tar_read(levels)
cfg    <- targets::tar_read(cfg)
cat("\n", sample_description(levels, cfg), "\n", sep = "")
cat("Outputs in", cfg$data$output_dir, "\n")
