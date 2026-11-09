# helper-project.R: shared setup for all tests (testthat runs helper-*.R files first).
#
# The project is not an R package, so the functions in R/ are sourced here.
# project_root() lets tests find config/ and docs/ whether they are run from the
# project root (scripts/run_tests.R) or from tests/testthat.

project_root <- function() {
  normalizePath(testthat::test_path("..", ".."), winslash = "/", mustWork = TRUE)
}

for (f in list.files(file.path(project_root(), "R"), pattern = "\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}
