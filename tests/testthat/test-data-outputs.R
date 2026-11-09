# test-data-outputs.R: helpers used by the pipeline (R/data_outputs.R).

test_that("raw_file_paths lists every raw file the run config names, once", {
  cfg <- yaml::read_yaml(file.path(project_root(), "config", "run.yml"))
  cfg$data$raw_dir <- file.path(project_root(), cfg$data$raw_dir)
  skip_if_not(file.exists(file.path(cfg$data$raw_dir, cfg$data$qd_vintage_file)),
              "raw data not present (raw data are not in git)")
  p <- raw_file_paths(cfg)
  expect_false(anyDuplicated(p) > 0)
  expect_true(all(file.exists(p)))
  expected <- c(cfg$data$qd_vintage_file, cfg$data$worldbank_file, cfg$data$yahoo_file,
                unlist(cfg$data$alfred_files), cfg$data$backfill$gs2$file)
  expect_true(all(file.path(cfg$data$raw_dir, expected) %in% p))
})

test_that("raw_file_paths stops on a missing file", {
  cfg <- list(data = list(raw_dir = tempdir(), qd_vintage_file = "does_not_exist.csv"))
  expect_error(raw_file_paths(cfg), "raw file not found")
})
