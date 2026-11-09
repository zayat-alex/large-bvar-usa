# run_tests.R: run the project's automated tests.
#
# Run from the project root (so .Rprofile activates renv):
#   Rscript scripts/run_tests.R
# The project is not an R package, so tests are run with testthat::test_dir();
# tests/testthat/helper-project.R sources the functions in R/.

testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)
