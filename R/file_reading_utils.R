#' File Reading Utilities for Genomic Data
#'
#' This file contains utility functions for reading various genomic data file formats
#' used in the Oscheius research project.

#' Read BED file containing assembly gap information
#'
#' @param file_path Path to the BED file
#' @return A data frame containing the parsed BED file with gap information
#' @export
read_assembly_gaps <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }
  
  # Check if this is the special tipulae file which has a different format
  is_tipulae_file <- grepl("nxOscTipu1\\.1\\.assembly\\.gaps\\.bed$", file_path)
  
  if (is_tipulae_file) {
    # Handle the special case for O. tipulae
    tryCatch({
      # Read the raw file content
      file_lines <- readLines(file_path)
      
      if (length(file_lines) == 0) {
        # Empty file
        warning("Empty file: ", file_path)
        return(data.frame(
          sequence_id = character(0),
          start = integer(0),
          end = integer(0),
          species_id = character(0),
          length = integer(0),
          stringsAsFactors = FALSE
        ))
      }
      
      # Parse the content - assuming "V 7444500 8204628" format
      parts <- strsplit(file_lines[1], "\\s+")[[1]]
      
      if (length(parts) >= 3) {
        # Create a data frame with the correct structure
        gaps <- data.frame(
          sequence_id = parts[1],
          start = as.integer(parts[2]),
          end = as.integer(parts[3]),
          stringsAsFactors = FALSE
        )
        
        # Add species identifier
        gaps$species_id <- extract_species_id(file_path)
        
        # Calculate gap length
        gaps$length <- gaps$end - gaps$start
        
        return(gaps)
      } else {
        warning("Invalid format in file: ", file_path)
        return(data.frame(
          sequence_id = character(0),
          start = integer(0),
          end = integer(0),
          species_id = character(0),
          length = integer(0),
          stringsAsFactors = FALSE
        ))
      }
    }, error = function(e) {
      warning("Error parsing tipulae gap file: ", file_path, "\n", e$message)
      return(data.frame(
        sequence_id = character(0),
        start = integer(0),
        end = integer(0),
        species_id = character(0),
        length = integer(0),
        stringsAsFactors = FALSE
      ))
    })
  } else {
    # Standard BED file format: chromosome start end [additional fields]
    tryCatch({
      gaps <- readr::read_tsv(file_path, 
                             col_names = c("sequence_id", "start", "end"), 
                             col_types = readr::cols(
                               sequence_id = readr::col_character(),
                               start = readr::col_integer(),
                               end = readr::col_integer()
                             ))
      
      # Add species identifier
      gaps$species_id <- extract_species_id(file_path)
      
      # Calculate gap length - handle empty data frames safely
      if (nrow(gaps) > 0) {
        gaps$length <- gaps$end - gaps$start
      } else {
        # If there are no rows, create an empty data frame with the right structure
        warning("No gap data found in file: ", file_path)
        gaps$length <- integer(0)
      }
      
      return(gaps)
    }, error = function(e) {
      warning("Error reading assembly gaps file: ", file_path, "\n", e$message)
      # Return empty data frame with correct structure
      empty_df <- data.frame(
        sequence_id = character(0),
        start = integer(0),
        end = integer(0),
        species_id = character(0),
        length = integer(0),
        stringsAsFactors = FALSE
      )
      return(empty_df)
    })
  }
}

#' Read NCRF telomere summary file
#'
#' @param file_path Path to the NCRF summary TSV file
#' @return A data frame containing the parsed telomere information
#' @export
read_telomere_summary <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }
  
  # Read the telomere summary file with tryCatch for robust error handling
  tryCatch({
    # First check the file
    file_lines <- readLines(file_path, n = 5)  # Read first few lines to inspect
    
    # Check if the file has content
    if (length(file_lines) == 0) {
      warning("Empty telomere file: ", file_path)
      return(create_empty_telomere_df_with_species(extract_species_id(file_path)))
    }
    
    # Get the species ID from the filename
    species_id <- extract_species_id(file_path)
    
    # Detect if this is an NCRF format file by checking for specific columns
    is_ncrf_format <- FALSE
    if (length(file_lines) > 0) {
      header_line <- file_lines[1]
      if (grepl("line\\s+motif\\s+seq\\s+start\\s+end\\s+strand", header_line)) {
        is_ncrf_format <- TRUE
      }
    }
    
    # Parse according to the detected format
    if (is_ncrf_format) {
      # Read as NCRF format with proper column names
      telomeres <- readr::read_tsv(file_path, 
                                 col_types = readr::cols(
                                   line = readr::col_integer(),
                                   motif = readr::col_character(),
                                   seq = readr::col_character(),
                                   start = readr::col_integer(),
                                   end = readr::col_integer(),
                                   strand = readr::col_character(),
                                   seqLen = readr::col_integer(),
                                   querybp = readr::col_integer(),
                                   mRatio = readr::col_character(),
                                   m = readr::col_integer(),
                                   mm = readr::col_integer(),
                                   i = readr::col_integer(),
                                   d = readr::col_integer()
                                 ))
      
      # Map NCRF columns to expected format
      mapped_telomeres <- telomeres %>%
        dplyr::mutate(
          species_id = species_id,
          sequence_id = seq,
          Position = start,
          Length = querybp,
          Score = as.numeric(sub("%$", "", mRatio))  # Remove % and convert to numeric
        )
      
      # Make sure all required columns exist
      mapped_telomeres <- mapped_telomeres %>%
        dplyr::select(species_id, sequence_id, Position, Length, strand, Score, motif, 
                     dplyr::everything())
      
      return(mapped_telomeres)
    } else {
      # Try standard format (fallback)
      cat("Attempting to read file as standard format:", file_path, "\n")
      # Use read.table for more robust parsing of potentially inconsistent files
      telomeres <- tryCatch({
        readr::read_tsv(file_path, 
                       col_types = readr::cols(.default = readr::col_character(),
                                              Score = readr::col_double(),
                                              Position = readr::col_integer(),
                                              Length = readr::col_integer()))
      }, error = function(e) {
        # Fallback to more flexible read.table if read_tsv fails
        warning("Using fallback parser for: ", file_path, " - ", e$message)
        df <- read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
        # Ensure proper column types
        if ("Position" %in% names(df)) df$Position <- as.integer(df$Position)
        if ("Length" %in% names(df)) df$Length <- as.integer(df$Length)
        if ("Score" %in% names(df)) df$Score <- as.numeric(df$Score)
        return(df)
      })
      
      # If the data frame is empty, return an empty data frame with correct structure
      if (nrow(telomeres) == 0) {
        warning("No telomere data in file: ", file_path)
        return(create_empty_telomere_df_with_species(species_id))
      }
      
      # Validate and standardize column names - case insensitive matching
      standardize_columns <- function(df) {
        col_map <- list(
          # Standard expected columns - map variations to standard names
          "sequence" = c("sequence", "seq", "chrom", "chromosome", "scaff", "scaffold"),
          "position" = c("position", "pos", "start", "begin"),
          "length" = c("length", "len", "size", "querybp"),
          "direction" = c("direction", "dir", "strand", "orientation")
        )
        
        # Create a new data frame with standardized column names
        new_df <- data.frame(
          species_id = rep(species_id, nrow(df)),
          stringsAsFactors = FALSE
        )
        
        # Find and map columns to standard names
        for (std_name in names(col_map)) {
          variations <- col_map[[std_name]]
          for (col in names(df)) {
            if (tolower(col) %in% variations) {
              new_df[[std_name]] <- df[[col]]
              break
            }
          }
        }
        
        # Handle special case for sequence - try to extract sequence_id if needed
        if ("Sequence" %in% names(df) && !"sequence" %in% names(new_df)) {
          new_df$sequence <- df$Sequence
        }
        
        # Copy any other columns that might be needed
        for (col in c("Position", "Length", "Score", "Direction")) {
          if (col %in% names(df) && !tolower(col) %in% names(new_df)) {
            new_df[[col]] <- df[[col]]
          }
        }
        
        return(new_df)
      }
      
      # Standardize column names
      telomeres <- standardize_columns(telomeres)
      
      # Add species identifier if not already present
      if (!"species_id" %in% names(telomeres)) {
        telomeres$species_id <- species_id
      }
      
      # Parse sequence position information - ensure sequence_id is available
      if ("sequence" %in% names(telomeres)) {
        # Extract sequence_id from sequence field
        telomeres$sequence_id <- sub(":.*", "", telomeres$sequence)
      } else {
        # If no sequence column, use a placeholder
        telomeres$sequence_id <- "unknown"
      }
      
      # Add strand information if not present
      if (!"strand" %in% names(telomeres)) {
        if ("Direction" %in% names(telomeres)) {
          telomeres$strand <- ifelse(grepl("Minus", telomeres$Direction), "-", "+")
        } else if ("direction" %in% names(telomeres)) {
          telomeres$strand <- ifelse(grepl("Minus|minus|-", telomeres$direction), "-", "+")
        } else {
          # Default to + strand if no direction info
          telomeres$strand <- "+"
        }
      }
      
      # Make sure required columns exist
      required_cols <- c("species_id", "sequence_id", "Position", "Length", "strand")
      for (col in required_cols) {
        if (!col %in% names(telomeres)) {
          # Try to map from standardized column names
          if (col == "Position" && "position" %in% names(telomeres)) {
            telomeres$Position <- telomeres$position
          } else if (col == "Length" && "length" %in% names(telomeres)) {
            telomeres$Length <- telomeres$length
          } else {
            # Create a default value
            telomeres[[col]] <- if (col %in% c("Position", "Length")) 0 else ""
            warning("Missing column '", col, "' in file: ", file_path, 
                   ". Using default values.")
          }
        }
      }
      
      return(telomeres)
    }
  }, error = function(e) {
    warning("Error reading telomere summary file: ", file_path, "\n", e$message)
    return(create_empty_telomere_df_with_species(extract_species_id(file_path)))
  })
}

#' Create an empty telomere data frame with correct structure
#' @return An empty data frame with telomere data structure
create_empty_telomere_df <- function() {
  data.frame(
    Sequence = character(0),
    Position = integer(0),
    Length = integer(0),
    Direction = character(0),
    Score = numeric(0),
    species_id = character(0),
    strand = character(0),
    sequence_id = character(0),
    position_info = character(0),
    stringsAsFactors = FALSE
  )
}

#' Create an empty telomere data frame with species ID
#' @param species_id Species ID to include in the empty data frame
#' @return An empty data frame with telomere data structure and the specified species_id
create_empty_telomere_df_with_species <- function(species_id) {
  df <- create_empty_telomere_df()
  if (!is.null(species_id) && length(species_id) > 0) {
    df$species_id <- species_id
  }
  return(df)
}

#' Read genome sequence size information from FAI file
#'
#' @param file_path Path to the FAI file
#' @return A data frame containing sequence IDs and their lengths
#' @export
read_sequence_sizes <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }
  
  # FAI file format: sequence_id length offset line_length line_length_bytes
  seq_sizes <- readr::read_tsv(file_path, 
                             col_names = c("sequence_id", "length", "offset", "line_length", "line_bytes"),
                             col_types = readr::cols(
                               sequence_id = readr::col_character(),
                               length = readr::col_integer(),
                               offset = readr::col_integer(),
                               line_length = readr::col_integer(),
                               line_bytes = readr::col_integer()
                             ))
  
  # Add species identifier
  seq_sizes$species_id <- extract_species_id(file_path)
  
  # Keep only needed columns
  seq_sizes <- seq_sizes %>%
    dplyr::select(species_id, sequence_id, length)
  
  return(seq_sizes)
}

#' Process all species sequence size data
#'
#' @param base_dir Base directory for data files
#' @return A data frame containing sequence sizes for all species
#' @export
process_all_sequence_sizes <- function(base_dir = NULL) {
  # Get file paths for all species
  file_paths <- get_species_file_paths("sequence_sizes", "*.primary.fa.gz.fai", base_dir)
  
  # Read and combine data from all files
  all_sizes <- purrr::map_df(unlist(file_paths), read_sequence_sizes)
  
  # Join with species metadata
  metadata <- get_oscheius_metadata()
  all_sizes <- all_sizes %>%
    dplyr::left_join(metadata, by = c("species_id" = "assembly_id"))
  
  return(all_sizes)
}

#' Process all species assembly gap data
#'
#' @param base_dir Base directory for data files
#' @return A data frame containing assembly gaps for all species
#' @export
process_all_assembly_gaps <- function(base_dir = NULL) {
  # Get file paths for all species
  file_paths <- get_species_file_paths("assembly_gaps", "*.assembly.gaps.bed", base_dir)
  
  # Read and combine data from all files
  all_gaps <- purrr::map_df(unlist(file_paths), read_assembly_gaps)
  
  # Join with species metadata
  metadata <- get_oscheius_metadata()
  all_gaps <- all_gaps %>%
    dplyr::left_join(metadata, by = c("species_id" = "assembly_id"))
  
  return(all_gaps)
}

#' Process all species telomere data
#'
#' Reads and combines telomere data from all species
#' @param base_dir Base directory (defaults to here::here("data"))
#' @return A data frame containing telomere data for all species
#' @export
process_all_telomere_data <- function(base_dir = NULL) {
  # Get file paths for all species
  file_paths <- get_species_file_paths("telomeres", "*NCRF_summary.tsv", base_dir)
  
  # Read and combine data from all files
  all_telomeres <- list()
  
  # Process each file individually to handle potential errors more gracefully
  for (species_id in names(file_paths)) {
    paths <- file_paths[[species_id]]
    for (path in paths) {
      if (file.exists(path)) {
        tryCatch({
          telomere_data <- read_telomere_summary(path)
          if (!is.null(telomere_data) && nrow(telomere_data) > 0) {
            all_telomeres[[length(all_telomeres) + 1]] <- telomere_data
          }
        }, error = function(e) {
          warning("Error processing telomere file ", path, ": ", e$message)
        })
      }
    }
  }
  
  # Combine all telomere data
  if (length(all_telomeres) > 0) {
    combined_telomeres <- dplyr::bind_rows(all_telomeres)
    
    # Add species metadata
    metadata <- get_oscheius_metadata()
    combined_telomeres <- combined_telomeres %>%
      dplyr::left_join(
        dplyr::select(metadata, assembly_id, species_name, strain, version),
        by = c("species_id" = "assembly_id")
      )
    
    return(combined_telomeres)
  } else {
    warning("No telomere data found for any species.")
    return(create_empty_telomere_df())
  }
}
