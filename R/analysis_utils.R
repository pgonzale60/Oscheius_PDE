#' Analysis Utilities for Genomic Data
#'
#' This file contains utility functions for analyzing and visualizing genomic data
#' for the Oscheius research project.

library(tidyverse)
library(ggplot2)
library(patchwork)

#' Classify telomere presence at chromosome ends
#'
#' @param telomere_data Data frame containing telomere information from read_telomere_summary
#' @param sequence_sizes Data frame containing sequence sizes from read_sequence_sizes
#' @param min_telomere_length Minimum acceptable telomere length in bp
#' @param max_distance_from_end Maximum acceptable distance from sequence end
#' @param min_motif_match Minimum match percentage for telomere motif (0-100)
#' @return A data frame with classification of telomere status for each sequence end
#' @export
classify_telomere_ends <- function(telomere_data, sequence_sizes, 
                                  min_telomere_length = 200, 
                                  max_distance_from_end = 50,
                                  min_motif_match = 95) {
  
  # Verify that necessary columns exist in input data frames
  required_telomere_cols <- c("species_id", "sequence_id", "Position", "Length", "strand", "Score")
  required_seq_cols <- c("species_id", "sequence_id", "length")
  
  missing_tel_cols <- setdiff(required_telomere_cols, names(telomere_data))
  missing_seq_cols <- setdiff(required_seq_cols, names(sequence_sizes))
  
  if (length(missing_tel_cols) > 0) {
    stop("Missing required columns in telomere_data: ", 
         paste(missing_tel_cols, collapse = ", "))
  }
  
  if (length(missing_seq_cols) > 0) {
    stop("Missing required columns in sequence_sizes: ", 
         paste(missing_seq_cols, collapse = ", "))
  }
  
  # Get unique sequence IDs and create a data frame with both ends of each sequence
  sequences <- sequence_sizes %>%
    dplyr::select(species_id, sequence_id, length) %>%
    dplyr::distinct()
  
  # Create a data frame with both ends (5' and 3') of each sequence
  sequence_ends <- sequences %>%
    tidyr::crossing(end = c("5_prime", "3_prime"))
  
  # Process telomere data to identify those near sequence ends
  telomere_at_ends <- tryCatch({
    # Clean up data for NCRF format - ensure there are no name conflicts
    # and columns have consistent types
    processed_telomeres <- telomere_data %>% 
      # Make sure we're working with the right sequence ID column
      dplyr::mutate(
        # Use sequence_id from telomere data, not 'seq' from NCRF format
        sequence_id = as.character(sequence_id),
        strand = as.character(strand),
        Position = as.integer(Position),
        Length = as.integer(Length),
        Score = as.numeric(Score)
      ) %>%
      # Calculate distance from either end of the sequence
      dplyr::mutate(
        # Join with sequence sizes to get sequence length for all calculations
        seq_length = sequences$length[match(paste(species_id, sequence_id), 
                                           paste(sequences$species_id, sequences$sequence_id))],
        # Calculate distance from either end of the sequence
        distance_from_5prime = Position,
        distance_from_3prime = dplyr::case_when(
          # Calculate distance from 3' end consistently for both strands
          !is.na(seq_length) ~ seq_length - (Position + Length),
          TRUE ~ NA_integer_
        )
      )
    
    # Apply multiple filters in sequence:
    # 1. Filter by minimum motif match score
    high_quality_telomeres <- processed_telomeres %>%
      dplyr::filter(Score >= min_motif_match)
    
    # 2. Filter by minimum telomere length
    long_enough_telomeres <- high_quality_telomeres %>%
      dplyr::filter(Length >= min_telomere_length)
    
    # 3. Filter by distance from ends
    near_ends <- long_enough_telomeres %>%
      dplyr::filter(!is.na(distance_from_5prime) & !is.na(distance_from_3prime)) %>%
      dplyr::filter(distance_from_5prime <= max_distance_from_end | 
                   distance_from_3prime <= max_distance_from_end) %>%
      # Determine which end the telomere belongs to
      dplyr::mutate(
        end = dplyr::case_when(
          # 5' end telomeres: near 5' end (low Position) 
          distance_from_5prime <= max_distance_from_end & strand == "+" ~ "5_prime",
          # 3' end telomeres: near 3' end (high Position)
          distance_from_3prime <= max_distance_from_end & strand == "-" ~ "3_prime",
          # Catch other cases that might be close to ends
          distance_from_5prime <= max_distance_from_end ~ "5_prime",
          distance_from_3prime <= max_distance_from_end ~ "3_prime",
          TRUE ~ NA_character_
        )
      ) %>%
      # Filter to keep only valid end classifications
      dplyr::filter(!is.na(end))
    
    # 4. For each chromosome end, keep only the longest telomere
    longest_telomeres <- near_ends %>%
      dplyr::group_by(species_id, sequence_id, end) %>%
      dplyr::slice_max(order_by = Length, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
    
    longest_telomeres
  }, error = function(e) {
    warning("Error processing telomere data: ", e$message)
    # Return empty data frame with correct structure
    data.frame(
      species_id = character(0),
      sequence_id = character(0),
      end = character(0),
      Length = integer(0),
      Position = integer(0),
      Score = numeric(0),
      distance_from_5prime = numeric(0),
      distance_from_3prime = numeric(0),
      stringsAsFactors = FALSE
    )
  })
  
  # Print column names for debugging
  cat("Columns in sequence_ends:", paste(names(sequence_ends), collapse=", "), "\n")
  cat("Columns in telomere_at_ends:", paste(names(telomere_at_ends), collapse=", "), "\n")
  
  # Join with sequence ends to identify which ends have telomeres
  telomere_classification <- tryCatch({
    # Select only the columns we need to avoid conflicts in the join
    minimal_telomere_data <- telomere_at_ends %>%
      dplyr::select(species_id, sequence_id, end, Length, Position, Score,
                   distance_from_5prime, distance_from_3prime)
    
    joined_data <- sequence_ends %>%
      dplyr::left_join(
        minimal_telomere_data,
        by = c("species_id", "sequence_id", "end")
      ) %>%
      # Classify telomere status
      dplyr::mutate(
        telomere_status = dplyr::case_when(
          # No telomere found at this end
          is.na(Length) ~ "no_telomere",
          # Telomere is too short (< min_telomere_length)
          Length < min_telomere_length ~ "too_short",
          # Telomere match quality is too low
          Score < min_motif_match ~ "low_quality_match",
          # Telomere is not at the very end of the sequence
          end == "5_prime" & distance_from_5prime > max_distance_from_end ~ "not_at_end",
          end == "3_prime" & distance_from_3prime > max_distance_from_end ~ "not_at_end",
          # If it passes all checks, it's a proper telomere
          TRUE ~ "proper_telomere"
        ),
        telomere_length = ifelse(is.na(Length), 0, Length),
        telomere_score = ifelse(is.na(Score), NA_real_, Score)
      )
    
    joined_data
  }, error = function(e) {
    warning("Error classifying telomeres: ", e$message)
    
    # Create a fallback dataset with error status
    sequence_ends %>%
      dplyr::mutate(
        telomere_status = "error_in_processing",
        telomere_length = 0,
        telomere_score = NA_real_
      )
  })
  
  # Return the classification result, ensuring column consistency
  result <- telomere_classification %>%
    dplyr::select(species_id, sequence_id, length, end, telomere_status, telomere_length, telomere_score)
  
  return(result)
}

#' Summarize telomere status by species
#'
#' @param telomere_classification Data frame from classify_telomere_ends
#' @return A summary data frame with counts by species and telomere status
#' @export
summarize_telomere_status <- function(telomere_classification) {
  # Check if species_name column exists in the input data
  has_species_name <- "species_name" %in% names(telomere_classification)
  
  # If species_name is missing, try to join with metadata
  if (!has_species_name) {
    metadata <- get_oscheius_metadata()
    telomere_classification <- telomere_classification %>%
      dplyr::left_join(metadata, by = c("species_id" = "assembly_id"))
  }
  
  # Now summarize the data
  summary <- telomere_classification %>%
    dplyr::group_by(species_id, species_name, telomere_status) %>%
    dplyr::summarize(
      count = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::group_by(species_id, species_name) %>%
    dplyr::mutate(
      total = sum(count),
      percentage = round(count / total * 100, 1)
    ) %>%
    dplyr::ungroup()
  
  return(summary)
}

#' Plot telomere status by species
#'
#' @param telomere_summary Data frame from summarize_telomere_status
#' @return A ggplot object
#' @export
plot_telomere_status <- function(telomere_summary) {
  # Reorder factors for better visualization
  telomere_summary <- telomere_summary %>%
    dplyr::mutate(
      telomere_status = factor(telomere_status, 
                             levels = c("proper_telomere", "too_short", "not_at_end", "no_telomere")),
      species_name = forcats::fct_reorder(species_name, species_id)
    )
  
  # Get species colors
  species_colors <- get_species_colors()
  
  # Create the plot
  p <- ggplot(telomere_summary, aes(x = species_name, y = percentage, fill = telomere_status)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_brewer(palette = "Set1", 
                     labels = c("Proper telomere", "Too short", "Not at end", "No telomere")) +
    labs(
      title = "Telomere Status at Chromosome Ends",
      x = "Species",
      y = "Percentage of Chromosome Ends",
      fill = "Telomere Status"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  return(p)
}

#' Analyze and plot assembly gap distribution
#'
#' @param gap_data Data frame containing assembly gap information
#' @param bin_width Width of bins for histogram
#' @return A list containing summary statistics and plots
#' @export
analyze_assembly_gaps <- function(gap_data, bin_width = 100) {
  # Calculate summary statistics
  gap_summary <- gap_data %>%
    dplyr::group_by(species_id, species_name) %>%
    dplyr::summarize(
      total_gaps = dplyr::n(),
      total_gap_length = sum(length),
      mean_gap_length = mean(length),
      median_gap_length = median(length),
      min_gap_length = min(length),
      max_gap_length = max(length),
      .groups = "drop"
    )
  
  # Create histogram of gap lengths
  gap_histogram <- ggplot(gap_data, aes(x = length, fill = species_id)) +
    geom_histogram() +
    facet_wrap(~ species_name, scales = "free_y") +
    scale_fill_manual(values = get_species_colors()) +
    labs(
      title = "Distribution of Assembly Gap Lengths",
      x = "Gap Length (bp)",
      y = "Count"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Return results
  return(list(
    summary = gap_summary,
    histogram = gap_histogram
  ))
}

#' Create telomere end maps for each species
#'
#' @param telomere_classification Data frame from classify_telomere_ends
#' @param sequence_sizes Data frame containing sequence sizes
#' @return A list of plots showing telomere status at chromosome ends
#' @export
create_telomere_end_maps <- function(telomere_classification, sequence_sizes) {
  # Check if species_name is present in telomere_classification
  if (!"species_name" %in% names(telomere_classification)) {
    metadata <- get_oscheius_metadata()
    telomere_classification <- telomere_classification %>%
      dplyr::left_join(metadata, by = c("species_id" = "assembly_id"))
  }
  
  # Get only chromosomes (exclude scaffolds, etc.)
  chromosomes <- sequence_sizes %>%
    # Typically, chromosomes are named differently than scaffolds
    dplyr::filter(grepl("^(chr|chromosome|linkage_group|LG)", sequence_id, ignore.case = TRUE) |
                 sequence_id %in% c(as.character(1:20), as.character(as.roman(1:20))))
  
  # Filter classification data to include only chromosomes
  chrom_telomeres <- telomere_classification %>%
    dplyr::filter(paste(species_id, sequence_id) %in% paste(chromosomes$species_id, chromosomes$sequence_id))
  
  # Create a plot for each species
  species_ids <- unique(chrom_telomeres$species_id)
  
  plots <- list()
  
  for (sp_id in species_ids) {
    species_data <- chrom_telomeres %>%
      dplyr::filter(species_id == sp_id)
    
    # Get sequence length info for this species
    sp_sequence_sizes <- sequence_sizes %>%
      dplyr::filter(species_id == sp_id)
    
    # Join with length data
    species_data <- species_data %>%
      dplyr::left_join(
        dplyr::select(sp_sequence_sizes, species_id, sequence_id, length),
        by = c("species_id", "sequence_id")
      )
    
    # Order the sequences by length (safely, checking lengths match first)
    seq_ids <- unique(species_data$sequence_id)
    if (length(seq_ids) > 0) {
      # Get species name
      sp_name <- unique(species_data$species_name)[1]
      if (is.na(sp_name) || length(sp_name) == 0) {
        sp_name <- sp_id  # Fallback to ID if name not available
      }
      
      # Convert "end" to a factor with specific order to ensure 5' is on left and 3' on right
      species_data$end <- factor(species_data$end, levels = c("5_prime", "3_prime"))
      
      # Create the plot 
      p <- ggplot(species_data, aes(x = end, y = sequence_id, fill = telomere_status)) +
        geom_tile(color = "white") +
        scale_fill_manual(
          values = c(
            "proper_telomere" = "darkgreen", 
            "too_short" = "yellow", 
            "not_at_end" = "orange", 
            "no_telomere" = "red",
            "error_in_processing" = "gray50"
          ),
          name = "Status",
          labels = function(breaks) {
            sapply(breaks, function(x) {
              switch(x,
                "proper_telomere" = "Proper telomere",
                "too_short" = "Too short",
                "not_at_end" = "Not at end",
                "no_telomere" = "No telomere",
                "error_in_processing" = "Error in processing",
                x # Default case returns the original value
              )
            })
          }
        ) +
        scale_x_discrete(
          labels = c("5_prime" = "5' End", "3_prime" = "3' End"),
          position = "bottom"
        ) +
        labs(
          title = paste("Telomere Status for", sp_name),
          x = "Chromosome End",
          y = "Chromosome",
          fill = "Status"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          panel.grid = element_blank()
        )
      
      plots[[sp_id]] <- p
    }
  }
  
  return(plots)
}
