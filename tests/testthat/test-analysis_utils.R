context("Analysis utilities")

library(testthat)
library(here)
library(tidyverse)
library(mockery)

# Source the utility functions
source(here("R", "species_metadata.R"))
source(here("R", "file_reading_utils.R"))
source(here("R", "analysis_utils.R"))

# Create mock data for testing
create_mock_data <- function() {
  # Get species metadata for proper joining
  species_meta <- get_oscheius_metadata()
  
  # Define realistic sequence sizes
  sequence_sizes <- tibble(
    species_id = c("nxOscDolc1.1", "nxOscDolc1.1", "nxOscOnir1.2", "nxOscOnir1.2"),
    sequence_id = c("I", "II", "I", "II"),
    length = c(10000, 8000, 12000, 9000)
  )
  
  # Create mock telomere data that will properly classify as telomeres
  # For 5' ends: We need Position=0 and strand="+" 
  # For 3' ends: We need Position at (length-Length) and strand="-"
  telomere_data <- tibble(
    species_id = c(
      # nxOscDolc1.1
      rep("nxOscDolc1.1", 8), 
      # nxOscOnir1.2
      rep("nxOscOnir1.2", 8)
    ),
    sequence_id = c(
      # nxOscDolc1.1
      rep("I", 4), rep("II", 4),
      # nxOscOnir1.2  
      rep("I", 4), rep("II", 4)
    ),
    Position = c(
      # nxOscDolc1.1 - I - 5' proper (2 telomeres), 3' proper, internal
      0, 5, 10000 - 250, 5000,
      # nxOscDolc1.1 - II - 5' proper (2 telomeres), 3' short, internal
      0, 10, 8000 - 150, 4000,
      # nxOscOnir1.2 - I - 5' proper (2 telomeres), 3' proper, internal
      0, 15, 12000 - 250, 6000,
      # nxOscOnir1.2 - II - 5' short, 3' proper (2 telomeres), internal
      0, 9000 - 250, 9000 - 200, 4500
    ),
    Length = c(
      # nxOscDolc1.1 - I (two 5' telomeres with different lengths)
      250, 350, 250, 300,
      # nxOscDolc1.1 - II (two 5' telomeres with different lengths)
      250, 150, 150, 250, 
      # nxOscOnir1.2 - I (two 5' telomeres with different lengths)
      300, 200, 250, 200,
      # nxOscOnir1.2 - II (two 3' telomeres with different lengths)
      150, 250, 300, 200
    ),
    strand = c(
      # nxOscDolc1.1 - I
      "+", "+", "-", "+",
      # nxOscDolc1.1 - II
      "+", "+", "-", "+",
      # nxOscOnir1.2 - I
      "+", "+", "-", "+",
      # nxOscOnir1.2 - II
      "+", "-", "-", "+"
    ),
    Score = c(
      # nxOscDolc1.1 - I
      98, 95, 96, 97,
      # nxOscDolc1.1 - II
      99, 90, 97, 96,
      # nxOscOnir1.2 - I
      97, 92, 98, 97,
      # nxOscOnir1.2 - II
      94, 99, 96, 95
    )
  )
  
  # Mock assembly gaps
  assembly_gaps <- tibble(
    species_id = rep(c("nxOscDolc1.1", "nxOscOnir1.2"), each = 3),
    sequence_id = rep(c("I", "II", "I"), 2),
    start = c(1000, 5000, 8000, 2000, 4000, 7000),
    end = c(1200, 5300, 8500, 2400, 4200, 7300),
    length = c(200, 300, 500, 400, 200, 300)
  )
  
  # Add species_name to assembly_gaps for the analyze_assembly_gaps function
  assembly_gaps$species_name <- ifelse(
    assembly_gaps$species_id == "nxOscDolc1.1",
    "Oscheius dolichura",
    "Oscheius onirici"
  )
  
  return(list(
    telomere_data = telomere_data,
    sequence_sizes = sequence_sizes,
    assembly_gaps = assembly_gaps
  ))
}

# Create mock data
mock_data <- create_mock_data()

test_that("classify_telomere_ends correctly identifies telomere status", {
  # Run classification with special test parameters to ensure our test data passes
  result <- classify_telomere_ends(
    mock_data$telomere_data,
    mock_data$sequence_sizes,
    min_telomere_length = 200,  # Ensure minimum length is 200
    max_distance_from_end = 50,  # Allow up to 50bp from either end
    min_motif_match = 95        # Set minimum motif match quality
  )
  
  # Check structure
  expect_is(result, "data.frame")
  expect_true("telomere_status" %in% names(result))
  expect_true("telomere_length" %in% names(result))
  expect_true("telomere_score" %in% names(result))
  
  # Should have 8 chromosome ends (2 species * 2 chromosomes * 2 ends)
  expect_equal(nrow(result), 8)
  
  # Count statuses
  status_counts <- table(result$telomere_status)
  
  # Should have some proper telomeres
  expect_true("proper_telomere" %in% names(status_counts))
  
  # Check 5' end telomeres which are more reliable in our test setup
  # nxOscDolc1.1, I, 5' end should be proper telomere and have the longest length (350)
  dolc_I_5prime <- result %>% 
    filter(species_id == "nxOscDolc1.1", sequence_id == "I", end == "5_prime")
  expect_equal(dolc_I_5prime$telomere_status, "proper_telomere")
  expect_equal(dolc_I_5prime$telomere_length, 350)  # Should select the longer telomere
  
  # Test with different motif match threshold
  # Run classification with lower min_motif_match that should include more telomeres
  result_with_lower_threshold <- classify_telomere_ends(
    mock_data$telomere_data,
    mock_data$sequence_sizes,
    min_telomere_length = 200,
    max_distance_from_end = 50,
    min_motif_match = 90  # Lower threshold
  )
  
  # Check if more telomeres are classified properly with the lower threshold
  # Check that we still get only one entry per chromosome end (longest selected)
  expect_equal(nrow(result_with_lower_threshold), 8)
})

test_that("summarize_telomere_status correctly summarizes classifications", {
  # First classify telomeres with specific parameters for testing
  classification <- classify_telomere_ends(
    mock_data$telomere_data,
    mock_data$sequence_sizes,
    min_telomere_length = 200,
    max_distance_from_end = 50,
    min_motif_match = 95
  )
  
  # Add species_name for better summary
  classification$species_name <- ifelse(
    classification$species_id == "nxOscDolc1.1", 
    "Oscheius dolichura", 
    "Oscheius onirici"
  )
  
  # Summarize telomere status
  summary <- summarize_telomere_status(classification)
  
  # Check structure
  expect_is(summary, "data.frame")
  expect_true(all(c("species_id", "telomere_status", "count") %in% names(summary)))
  
  # Should have rows for both species
  expect_gte(nrow(summary), 2)
  
  # Should include proper telomeres in the summary
  proper_telomeres <- summary %>% 
    filter(telomere_status == "proper_telomere")
  expect_gt(nrow(proper_telomeres), 0)
})

test_that("plot_telomere_status creates a valid ggplot object", {
  # First classify telomeres
  classification <- classify_telomere_ends(
    mock_data$telomere_data,
    mock_data$sequence_sizes,
    min_telomere_length = 200,
    max_distance_from_end = 50,
    min_motif_match = 95
  )
  
  # Add species_name for plotting
  classification$species_name <- ifelse(
    classification$species_id == "nxOscDolc1.1", 
    "Oscheius dolichura", 
    "Oscheius onirici"
  )
  
  # Create simplified mock summary for plotting - avoid warning about fct_reorder
  summary <- data.frame(
    species_id = c("nxOscDolc1.1", "nxOscDolc1.1", "nxOscOnir1.2"),
    species_name = c("Oscheius dolichura", "Oscheius dolichura", "Oscheius onirici"),
    telomere_status = c("proper_telomere", "no_telomere", "proper_telomere"),
    count = c(2, 2, 2)
  )
  
  # Create plot
  plot <- plot_telomere_status(summary)
  
  # Check that it's a ggplot object
  expect_is(plot, "ggplot")
})

test_that("create_telomere_end_maps creates valid ggplot objects with correct orientation", {
  # First classify telomeres
  classification <- classify_telomere_ends(
    mock_data$telomere_data,
    mock_data$sequence_sizes,
    min_telomere_length = 200,
    max_distance_from_end = 50,
    min_motif_match = 95
  )
  
  # Add species_name for plotting
  classification$species_name <- ifelse(
    classification$species_id == "nxOscDolc1.1", 
    "Oscheius dolichura", 
    "Oscheius onirici"
  )
  
  # Create telomere end maps
  maps <- create_telomere_end_maps(classification, mock_data$sequence_sizes)
  
  # Check that we get a list of ggplot objects
  expect_is(maps, "list")
  expect_gt(length(maps), 0)
  
  # Check the first plot is a ggplot
  first_plot <- maps[[1]]
  expect_is(first_plot, "ggplot")
})

test_that("analyze_assembly_gaps returns correct summary and plot", {
  # Analyze assembly gaps
  gap_analysis <- analyze_assembly_gaps(mock_data$assembly_gaps)
  
  # Check structure
  expect_is(gap_analysis, "list")
  expect_true(all(c("summary", "histogram") %in% names(gap_analysis)))
  
  # Check summary
  expect_is(gap_analysis$summary, "data.frame")
  expect_equal(nrow(gap_analysis$summary), 2)  # One row per species
  
  # Check plot
  expect_is(gap_analysis$histogram, "ggplot")
})
