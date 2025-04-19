#' Oscheius Species Metadata
#'
#' This file contains metadata and utility functions for working with the 5 Oscheius species
#' used throughout this research project.

#' Get a data frame containing Oscheius species metadata
#'
#' @return A data frame with columns for assembly_id, species_name, strain, and version
#' @export
get_oscheius_metadata <- function() {
  metadata <- data.frame(
    assembly_id = c(
      "nxOscDolc1.1",
      "nxOscOnir1.2",
      "nxOscSper1.1",
      "nxOscSpeu1.1", 
      "nxOscTipu1.1"
    ),
    species_name = c(
      "Oscheius dolichura",
      "Oscheius onirici",
      "Oscheius sp. DF5120",
      "Oscheius sp. JU1382",
      "Oscheius tipulae"
    ),
    strain = c(
      "PS1017",
      "PS2068",
      "DF5120",
      "JU1382", 
      "CEW1"
    ),
    version = c(
      "1.1",
      "1.2",
      "1.1",
      "1.1",
      "1.1"
    ),
    stringsAsFactors = FALSE
  )
  
  return(metadata)
}

#' Get file paths for all species for a given data type
#'
#' @param data_type Type of data to retrieve (assembly_gaps, telomeres, sequence_sizes)
#' @param file_pattern Pattern to match files
#' @param base_dir Base directory (defaults to here::here("data"))
#' @return A named list of file paths, with assembly_ids as names
#' @export
get_species_file_paths <- function(data_type, file_pattern, base_dir = NULL) {
  if (is.null(base_dir)) {
    base_dir <- here::here("data")
  }
  
  metadata <- get_oscheius_metadata()
  
  # Create full path patterns based on data type
  paths <- switch(data_type,
    "assembly_gaps" = file.path(base_dir, "assembly", "assembly_gaps", file_pattern),
    "telomeres" = file.path(base_dir, "assembly", "telomeres", file_pattern),
    "sequence_sizes" = file.path(base_dir, "genome_features", "sequence_sizes", file_pattern),
    stop("Unknown data type: ", data_type)
  )
  
  # Set up the pattern differently for telomeres
  if (data_type == "telomeres") {
    # Special handling for telomere files
    species_paths <- lapply(metadata$assembly_id, function(id) {
      # Create a pattern that will match any file containing both the species ID and NCRF_summary.tsv
      pattern_dir <- dirname(paths)
      # For telomeres, we want files that contain both the species ID and the NCRF part
      pattern_regex <- paste0("^", id, ".*NCRF_summary\\.tsv$")
      
      # List all files matching the pattern
      files <- list.files(pattern_dir, pattern = pattern_regex, full.names = TRUE)
      
      # Debug output
      if (length(files) == 0) {
        cat("Debug: No files found for species", id, "with pattern", pattern_regex, "in", pattern_dir, "\n")
        # Try listing all files to see what's actually there
        all_files <- list.files(pattern_dir, full.names = FALSE)
        if (length(all_files) > 0) {
          cat("Debug: Files in directory:", paste(all_files, collapse = ", "), "\n")
        } else {
          cat("Debug: Directory is empty\n")
        }
      }
      
      files
    })
  } else {
    # Standard handling for other data types
    species_paths <- lapply(metadata$assembly_id, function(id) {
      pattern <- gsub("\\*", id, paths)
      list.files(path = dirname(pattern), pattern = basename(pattern), full.names = TRUE)
    })
  }
  
  names(species_paths) <- metadata$assembly_id
  
  return(species_paths)
}

#' Get colors for consistent species visualization
#'
#' @return A named vector of colors for each Oscheius species
#' @export
get_species_colors <- function() {
  metadata <- get_oscheius_metadata()
  
  # Define a color palette
  colors <- c(
    "nxOscDolc1.1" = "#E41A1C", # red
    "nxOscOnir1.2" = "#377EB8", # blue
    "nxOscSper1.1" = "#4DAF4A", # green
    "nxOscSpeu1.1" = "#984EA3", # purple
    "nxOscTipu1.1" = "#FF7F00"  # orange
  )
  
  return(colors)
}

#' Get species ID from file name
#'
#' @param file_name Name or path of file
#' @return Character string with the species assembly ID
#' @export
extract_species_id <- function(file_name) {
  # Extract all possible species IDs
  metadata <- get_oscheius_metadata()
  
  # Find which species ID is in the filename
  for (id in metadata$assembly_id) {
    if (grepl(id, basename(file_name))) {
      return(id)
    }
  }
  
  warning("Could not extract species ID from filename: ", file_name)
  return(NA)
}
