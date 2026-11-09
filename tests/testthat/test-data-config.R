# test-data-config.R: config/series.csv must agree with docs/data.md.
#
# The series table of docs/data.md (section 2) defines the 31 series.
# config/series.csv is the machine-readable copy the pipeline
# uses. These tests parse the markdown table and compare it cell by
# cell, so the two files cannot drift apart silently.

# Read the series table of docs/data.md: the lines that start with "| <number> |".
# Returns a data.frame with the table's columns, backticks and bold removed.
read_doc_series_table <- function(doc_path) {
  lines <- readLines(doc_path, encoding = "UTF-8")
  rows  <- grep("^\\| *[0-9]+ *\\|", lines, value = TRUE)
  cells <- lapply(strsplit(rows, "|", fixed = TRUE), function(x) {
    x <- trimws(x[-1])                  # drop the empty cell before the first "|"
    gsub("`|\\*\\*", "", x)             # remove markdown code and bold marks
  })
  table <- as.data.frame(do.call(rbind, cells), stringsAsFactors = FALSE)
  names(table) <- c("row", "short_name", "variable", "source", "id",
                    "transform", "prior", "notes")
  table$row <- as.integer(table$row)
  table
}

series <- read_series_config(file.path(project_root(), "config", "series.csv"))
doc_table <- read_doc_series_table(file.path(project_root(), "docs", "data.md"))

test_that("series config has the 31 variables with the documented counts", {
  expect_equal(nrow(series), 31)
  expect_equal(sum(series$transform == "log"), 22)   # docs/data.md: "22 log, 9 raw"
  expect_equal(sum(series$transform == "raw"), 9)
  # Only realized volatility has a white-noise prior mean (docs/data.md, section 1).
  expect_equal(series$short_name[series$prior == "WN"], "rvol")
})

test_that("series config matches the docs/data.md series table cell by cell", {
  expect_equal(nrow(doc_table), 31)
  for (col in c("row", "short_name", "variable", "source", "id", "transform", "prior")) {
    expect_identical(series[[col]], doc_table[[col]], label = paste0("series.csv$", col))
  }
})

test_that("API aggregation rules follow the series notes", {
  api <- series[series$source == "API", ]
  expect_setequal(api$id, c("B823RA3Q086SBEA", "HCOMPBS", "GS2", "WTISPLC"))
  # Quarterly series need no aggregation; GS2 and WTISPLC are monthly, averaged to quarterly.
  expect_equal(api$api_aggregation[api$id == "B823RA3Q086SBEA"], "none")
  expect_equal(api$api_aggregation[api$id == "HCOMPBS"], "none")
  expect_equal(api$api_aggregation[api$id == "GS2"], "avg")
  expect_equal(api$api_aggregation[api$id == "WTISPLC"], "avg")
})

test_that("read_series_config rejects broken files", {
  bad <- series
  bad$transform[1] <- "diff"
  path <- tempfile(fileext = ".csv")
  utils::write.csv(bad, path, row.names = FALSE, na = "")
  expect_error(read_series_config(path), "unknown value")

  bad <- series
  bad$api_aggregation[1] <- "avg"          # rgdp is a QD row
  utils::write.csv(bad, path, row.names = FALSE, na = "")
  expect_error(read_series_config(path), "only allowed on API rows")
})
