context("Species metadata utilities")

# Setup mock data and environment
library(testthat)
library(here)

# Source the utility functions
source(here("R", "species_metadata.R"))

test_that("get_oscheius_metadata returns correct structure", {
  metadata <- get_oscheius_metadata()
  
  # Check structure
  expect_is(metadata, "data.frame")
  expect_equal(nrow(metadata), 5)
  expect_equal(ncol(metadata), 4)
  
  # Check column names
  expected_cols <- c("assembly_id", "species_name", "strain", "version")
  expect_equal(colnames(metadata), expected_cols)
  
  # Check that each species has required data
  expect_true(all(!is.na(metadata$assembly_id)))
  expect_true(all(!is.na(metadata$species_name)))
  expect_true(all(!is.na(metadata$strain)))
  expect_true(all(!is.na(metadata$version)))
  
  # Check specific values
  expect_true("nxOscDolc1.1" %in% metadata$assembly_id)
  expect_true("Oscheius dolichura" %in% metadata$species_name)
})

test_that("get_species_colors returns correctly formatted colors", {
  colors <- get_species_colors()
  
  # Check structure
  expect_is(colors, "character")
  expect_equal(length(colors), 5)
  
  # Check format of colors (should be hex codes)
  expect_true(all(grepl("^#[0-9A-F]{6}$", toupper(colors))))
  
  # Check that all species have a color
  metadata <- get_oscheius_metadata()
  expect_equal(sort(names(colors)), sort(metadata$assembly_id))
})

test_that("extract_species_id correctly identifies species in filenames", {
  # Test with full paths
  expect_equal(extract_species_id("/path/to/nxOscDolc1.1.assembly.gaps.bed"), "nxOscDolc1.1")
  expect_equal(extract_species_id("/data/assembly/telomeres/nxOscSper1.1.assembly.NCRF_summary.tsv"), "nxOscSper1.1")
  
  # Test with just filenames
  expect_equal(extract_species_id("nxOscOnir1.2.primary.fa.gz.fai"), "nxOscOnir1.2")
  
  # Test with unknown species
  expect_warning(extract_species_id("unknown_species.file"))
})

test_that("get_species_file_paths returns paths for all species", {
  # Mock the list.files function
  mockery::stub(get_species_file_paths, "list.files", function(...) {
    return("mocked_path.bed") 
  })
  
  paths <- get_species_file_paths("assembly_gaps", "*.bed")
  
  # Check structure
  expect_is(paths, "list")
  expect_equal(length(paths), 5)
  
  # Check names
  metadata <- get_oscheius_metadata()
  expect_equal(sort(names(paths)), sort(metadata$assembly_id))
})
