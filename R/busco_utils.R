# R/busco_utils.R
# Utilities to load BUSCO full table outputs and standardize chromosome identifiers

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

#' Read BUSCO odb10 full table tsv
#'
#' @param busco_file Path to BUSCO full table TSV
#' @return Tibble containing Busco_id, Sequence, start, end for Complete BUSCOs
read_busco <- function(busco_file) {
  read_tsv(
    busco_file,
    col_names = c("Busco_id", "Status", "Sequence",
                  "start", "end", "strand", "Score", "Length",
                  "OrthoDB_url", "Description"),
    col_types = c("ccciicdicc"),
    comment = "#",
    show_col_types = FALSE
  ) %>%
    filter(Status == "Complete") %>%
    select(Busco_id, Sequence, start, end)
}

#' Standardize Oscheius tipulae chromosome accessions
#'
#' Maps NCBI GenBank accessions CP059028.1-CP059033.1 to Roman numerals I-V, X.
#'
#' @param df Tibble containing 'assembly' and 'Sequence' columns
#' @return Standardized tibble
standardize_chromosomes <- function(df) {
  ncbi_map <- c(
    "CP059028.1" = "I",
    "CP059029.1" = "II",
    "CP059030.1" = "III",
    "CP059031.1" = "IV",
    "CP059032.1" = "V",
    "CP059033.1" = "X"
  )
  
  df %>%
    mutate(
      Sequence = ifelse(assembly == "nxOscTipu1.1" & Sequence %in% names(ncbi_map),
                        ncbi_map[Sequence], Sequence),
      Sequence = sub("^chr_", "", Sequence)
    )
}
