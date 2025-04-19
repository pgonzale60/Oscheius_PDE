#!/usr/bin/env Rscript
# Test Runner for Oscheius_PDE
# This script makes it easy to run tests with various options

# Load required packages
if (!requireNamespace("testthat", quietly = TRUE)) {
  install.packages("testthat")
}
if (!requireNamespace("mockery", quietly = TRUE)) {
  install.packages("mockery")
}

library(testthat)
library(here)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
verbose <- "--verbose" %in% args
filter_pattern <- NULL

# Check for specific test pattern
for (arg in args) {
  if (grepl("^--filter=", arg)) {
    filter_pattern <- gsub("^--filter=", "", arg)
  }
}

# Print header
cat("\n===== Running Oscheius_PDE Tests =====\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version:", R.version.string, "\n\n")

# Ensure we're in the project root
if (!file.exists(here("tests", "testthat.R"))) {
  stop("Tests not found. Make sure you're running this from the project root.")
}

# Run the tests and capture the results
if (!is.null(filter_pattern)) {
  cat("Running tests matching pattern:", filter_pattern, "\n\n")
  test_results <- testthat::test_dir(
    here("tests", "testthat"),
    filter = filter_pattern,
    reporter = if (verbose) "progress" else "summary",
    stop_on_failure = FALSE
  )
} else {
  cat("Running all tests\n\n")
  test_results <- testthat::test_dir(
    here("tests", "testthat"),
    reporter = if (verbose) "progress" else "summary",
    stop_on_failure = FALSE
  )
}

# Print summary
cat("\n===== Test Results =====\n")
cat("All tests completed.\n")

# Check if there were any failures or errors
if (length(test_results$failed) > 0 || length(test_results$error) > 0) {
  total_issues <- length(test_results$failed) + length(test_results$error)
  cat("\nTest result: FAILED (", total_issues, " issues)\n", sep = "")
  quit(status = 1)  # Exit with non-zero status code to indicate failure
} else {
  cat("\nAll tests passed successfully!\n")
  quit(status = 0)  # Exit with zero status code to indicate success
}
