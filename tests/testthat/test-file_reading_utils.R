context("File reading utilities")

library(testthat)
library(here)
library(mockery)

# Source the utility functions
source(here("R", "species_metadata.R"))
source(here("R", "file_reading_utils.R"))

# Create mock data directory for tests
mock_data_dir <- file.path(tempdir(), "mock_oscheius_data")
dir.create(mock_data_dir, recursive = TRUE, showWarnings = FALSE)

# Helper function to create mock data files
create_mock_files <- function() {
  # Create mock directories
  gap_dir <- file.path(mock_data_dir, "assembly", "assembly_gaps")
  tel_dir <- file.path(mock_data_dir, "assembly", "telomeres")
  seq_dir <- file.path(mock_data_dir, "genome_features", "sequence_sizes")
  
  dir.create(gap_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tel_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(seq_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create a mock assembly gaps file
  gap_file <- file.path(gap_dir, "nxOscDolc1.1.assembly.gaps.bed")
  gap_content <- "chr1\t1000\t2000\nchr1\t3000\t3500\nchr2\t500\t700"
  writeLines(gap_content, gap_file)
  
  # Create a mock empty gap file
  empty_gap_file <- file.path(gap_dir, "nxOscTipu1.1.assembly.gaps.bed")
  writeLines("", empty_gap_file)
  
  # Create a mock telomere file
  tel_file <- file.path(tel_dir, "nxOscDolc1.1.assembly.NCRF_summary.tsv")
  tel_content <- "Sequence\tPattern\tSense\tPosition\tLength\tScore\tDirection\nchrom1:1-1000\tTTAGGC\tPlus\t10\t250\t90\tPlus\nchrom2:1-1000\tTTAGGC\tMinus\t5\t180\t85\tMinus"
  writeLines(tel_content, tel_file)
  
  # Create a mock sequence size file
  seq_file <- file.path(seq_dir, "nxOscDolc1.1.primary.fa.gz.fai")
  seq_content <- "chr1\t10000\t0\t80\t81\nchr2\t8000\t10081\t80\t81\nchr3\t12000\t18161\t80\t81"
  writeLines(seq_content, seq_file)
  
  return(list(
    gap_file = gap_file,
    empty_gap_file = empty_gap_file,
    tel_file = tel_file,
    seq_file = seq_file
  ))
}

# Create mock files for tests
mock_files <- create_mock_files()

test_that("read_assembly_gaps correctly parses BED files", {
  # Test with valid file
  gaps <- read_assembly_gaps(mock_files$gap_file)
  
  # Check structure
  expect_is(gaps, "data.frame")
  expect_equal(nrow(gaps), 3)
  expect_equal(ncol(gaps), 5)
  
  # Check column names
  expected_cols <- c("sequence_id", "start", "end", "species_id", "length")
  expect_equal(sort(colnames(gaps)), sort(expected_cols))
  
  # Check values
  expect_equal(gaps$length[1], 1000)  # 2000 - 1000
  expect_equal(gaps$species_id[1], "nxOscDolc1.1")
  
  # Test with empty file
  expect_warning(empty_gaps <- read_assembly_gaps(mock_files$empty_gap_file))
  expect_equal(nrow(empty_gaps), 0)
  expect_equal(ncol(empty_gaps), 5)  # Should still have the right structure
  
  # Test with non-existent file
  expect_error(read_assembly_gaps("non_existent_file.bed"))
})

test_that("read_telomere_summary correctly parses telomere files", {
  # Test with valid file
  telomeres <- read_telomere_summary(mock_files$tel_file)
  
  # Check structure
  expect_is(telomeres, "data.frame")
  expect_equal(nrow(telomeres), 2)
  
  # Check that required columns are present
  required_cols <- c("species_id", "sequence_id", "Position", "Length", "strand")
  expect_true(all(required_cols %in% colnames(telomeres)))
  
  # Check values
  expect_equal(telomeres$species_id[1], "nxOscDolc1.1")
  expect_equal(telomeres$strand[1], "+")
  expect_equal(telomeres$Length[1], 250)
  
  # Test with non-existent file
  expect_error(read_telomere_summary("non_existent_file.tsv"))
})

test_that("read_sequence_sizes correctly parses FAI files", {
  # Test with valid file
  seq_sizes <- read_sequence_sizes(mock_files$seq_file)
  
  # Check structure
  expect_is(seq_sizes, "data.frame")
  expect_equal(nrow(seq_sizes), 3)
  
  # Check column names
  expected_cols <- c("species_id", "sequence_id", "length")
  expect_equal(sort(colnames(seq_sizes)), sort(expected_cols))
  
  # Check values
  expect_equal(seq_sizes$species_id[1], "nxOscDolc1.1")
  expect_equal(seq_sizes$length[1], 10000)
  
  # Test with non-existent file
  expect_error(read_sequence_sizes("non_existent_file.fai"))
})

test_that("process_all_sequence_sizes combines data from all species", {
  # Mock the get_species_file_paths function
  stub(process_all_sequence_sizes, "get_species_file_paths", 
       function(...) list(nxOscDolc1.1 = mock_files$seq_file))
  
  # Run the function
  all_sizes <- process_all_sequence_sizes(mock_data_dir)
  
  # Check structure
  expect_is(all_sizes, "data.frame")
  expect_true(nrow(all_sizes) > 0)
  
  # Check join with metadata
  expect_true("species_name" %in% colnames(all_sizes))
  expect_true("strain" %in% colnames(all_sizes))
})

# Clean up
unlink(mock_data_dir, recursive = TRUE)
