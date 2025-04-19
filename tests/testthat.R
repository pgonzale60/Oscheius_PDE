library(testthat)
library(here)

# Source the utility functions
source(here("R", "species_metadata.R"))
source(here("R", "file_reading_utils.R"))
source(here("R", "analysis_utils.R"))

# Run all tests
test_dir(here("tests", "testthat"))
