# R/synteny_utils.R
# Utilities for coordinate transformations and multispecies Oxford synteny plots

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

#' Transform chromosome sizes into cumulative multispecies coordinate offsets
sizes_to_multi_sp_index <- function(seq_sizes_tbl) {
  seq_sizes_tbl %>%
    arrange(multispecies_sequence) %>%
    mutate(
      multispecies_sequence = factor(multispecies_sequence),
      multispecies_index = cumsum(size) - size
    ) %>%
    select(multispecies_sequence, multispecies_index, size)
}

#' Shift chromosome-level coordinates to genome-wide cumulative offsets
busco_subset_and_recoord <- function(busco_tbl, multi_sp_index) {
  busco_tbl %>%
    inner_join(multi_sp_index, by = "multispecies_sequence") %>%
    mutate(
      multispecies_start = start + multispecies_index,
      multispecies_end = end + multispecies_index
    )
}
