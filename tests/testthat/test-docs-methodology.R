# test-docs-methodology.R: the settings in docs/methodology.md must agree
# with config/run.yml.
#
# docs/methodology.md states the prior, model and MCMC settings that
# determine the estimates (section 8.1), and the forecast, tilting and
# impulse-response settings of the engine (sections 9.6 and 10) and the
# export settings (section 15), one per
# table row, keyed by their path in config/run.yml. These tests parse those
# tables and compare every row with the configuration, and check that every
# setting of the covered sections appears in them, so the document cannot
# drift from the code's configuration silently.

# Read the settings table: the rows whose first cell is a dotted config path
# in backticks, e.g. "| `prior.lambda.mode` | 0.2 | ... |".
# Returns a data.frame: key (the path), value (the second cell, as text).
read_doc_settings <- function(doc_path) {
  lines <- readLines(doc_path, encoding = "UTF-8")
  rows  <- grep("^\\| *`[a-z_]+(\\.[a-z_]+)+` *\\|", lines, value = TRUE)
  cells <- lapply(strsplit(rows, "|", fixed = TRUE), function(x) trimws(x[2:3]))
  data.frame(key = gsub("`", "", vapply(cells, `[`, "", 1)),
             value = vapply(cells, `[`, "", 2), stringsAsFactors = FALSE)
}

# Value of a dotted path in the configuration list (NULL if absent).
config_value <- function(cfg, key) {
  x <- cfg
  for (part in strsplit(key, ".", fixed = TRUE)[[1]]) {
    if (!is.list(x) || is.null(x[[part]])) return(NULL)
    x <- x[[part]]
  }
  x
}

# Every leaf of a nested list, as dotted paths.
config_leaves <- function(x, prefix) {
  if (!is.list(x)) return(prefix)
  unlist(lapply(names(x), function(nm) config_leaves(x[[nm]], paste(prefix, nm, sep = "."))))
}

cfg <- yaml::read_yaml(file.path(project_root(), "config", "run.yml"))
doc <- read_doc_settings(file.path(project_root(), "docs", "methodology.md"))

test_that("every setting in the methodology table matches config/run.yml", {
  expect_gt(nrow(doc), 30)
  expect_false(any(duplicated(doc$key)))
  for (i in seq_len(nrow(doc))) {
    key <- doc$key[i]
    actual <- config_value(cfg, key)
    expect_false(is.null(actual), label = paste(key, "exists in config/run.yml"))
    if (is.null(actual)) next
    documented <- strsplit(doc$value[i], ", ", fixed = TRUE)[[1]]
    actual <- unlist(actual)
    if (is.numeric(actual)) {
      expect_equal(as.numeric(documented), as.numeric(actual), label = key)
    } else {
      expect_identical(documented, as.character(actual), label = key)
    }
  }
})

test_that("the methodology tables cover every setting of the estimates and the engine", {
  covered <- c(config_leaves(cfg$model, "model"), config_leaves(cfg$prior, "prior"),
               "covid.methods", "covid.cutoff", "mcmc.seed",
               config_leaves(cfg$mcmc$mode_search, "mcmc.mode_search"),
               config_leaves(cfg$mcmc$chain_start, "mcmc.chain_start"),
               config_leaves(cfg$mcmc$tuning, "mcmc.tuning"),
               config_leaves(cfg$mcmc$convergence, "mcmc.convergence"),
               config_leaves(cfg$mcmc$profiles$full, "mcmc.profiles.full"),
               config_leaves(cfg$forecast, "forecast"), config_leaves(cfg$irf, "irf"),
               "export.run", "export.horizon", "export.shown_in_levels")
  missing <- setdiff(covered, doc$key)
  expect_identical(missing, character(0))
})

test_that("the table parser and path lookup work on a known example", {
  path <- tempfile(fileext = ".md")
  writeLines(c("| Setting | Value | Meaning |", "|---|---|---|",
               "| `a.b` | 1, 2 | x |", "| not a key | 3 | y |"), path)
  d <- read_doc_settings(path)
  expect_equal(d$key, "a.b")
  expect_equal(d$value, "1, 2")
  expect_equal(config_value(list(a = list(b = 5)), "a.b"), 5)
  expect_null(config_value(list(a = list(b = 5)), "a.c"))
  expect_setequal(config_leaves(list(x = 1, y = list(z = 2)), "r"), c("r.x", "r.y.z"))
})
