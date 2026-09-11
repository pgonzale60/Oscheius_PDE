#!/usr/bin/env Rscript
# ==============================================================================
# Figure 6: Synteny comparisons link sites of PDE breakage to chromosomal rearrangement
# Manuscript: "Synteny comparisons link sites of PDE breakage to chromosomal rearrangement between Oscheius species"
# Panel: Oxford plot comparing O. dolichura and O. sp. JU1382 with Nigon units and break sites
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

# ------------------------------------------------------------------------------
# 1. Load Utilities & Species Metadata
# ------------------------------------------------------------------------------
source("R/plot_theme.R")
source("R/busco_utils.R")
source("R/synteny_utils.R")

tol_to_publicIds <- tibble(
  assembly = c("nxOscDolc1.1", "nxOscOnir1.2", "nxOscSper1.1", "nxOscSpeu1.1", "nxOscTipu1.1"),
  public_id = c("O. dolichura", "O. onirici", "O. sp. DF5120", "O. sp. JU1382", "O. tipulae")
)

nigonDict <- read_tsv("data/metadata/gene2Nigon_busco20200927.tsv.gz", col_types = "cc", show_col_types = FALSE)

# ------------------------------------------------------------------------------
# 2. Extract and Format Core Data
# ------------------------------------------------------------------------------
cat("Loading sequence sizes...\n")
fai_files <- list.files("data/sequence_sizes/", full.names = TRUE, pattern = "\\.fai$")
names(fai_files) <- sub(".+/([^/]+)\\.primary\\.fa\\.gz\\.fai$", "\\1", fai_files)

seq_sizes <- map_df(fai_files, read_tsv,
                    col_names = c("Sequence", "size", "L1", "L2", "L3"),
                    col_types = "ciiii",
                    show_col_types = FALSE,
                    .id = "assembly") %>%
  dplyr::select(assembly, Sequence, size) %>%
  dplyr::filter(!grepl("loc|MT", Sequence)) %>%
  dplyr::mutate(
    short_species_id = sub("nxOsc([^0-9]+).+", "\\1", assembly),
    multispecies_sequence = paste0(short_species_id, "_", Sequence)
  )

cat("Loading BUSCOs...\n")
buscofiles <- list.files("data/busco/", full.names = TRUE, pattern = "nematoda_odb10_full_table\\.tsv$")
names(buscofiles) <- sub(".+/([^/]+)\\.primary_nematoda_odb10_full_table\\.tsv$", "\\1", buscofiles)

buscos <- map_df(buscofiles, read_busco, .id = "assembly") %>%
  standardize_chromosomes() %>%
  dplyr::left_join(tol_to_publicIds, by = "assembly") %>%
  dplyr::mutate(
    short_species_id = sub("nxOsc([^0-9]+).+", "\\1", assembly),
    multispecies_sequence = paste0(short_species_id, "_", Sequence)
  )

cat("Loading break sites...\n")
break_files <- list.files("data/curated_GRS/", full.names = TRUE, pattern = "\\.chr_diminutions_sites\\.tsv$")
names(break_files) <- sub(".+/([^/]+)\\.chr_diminutions_sites\\.tsv$", "\\1", break_files)

break_sites <- map_df(break_files, read_tsv, show_col_types = FALSE, .id = "assembly") %>%
  dplyr::rename(Sequence = chr, pos = diminution_pos) %>%
  dplyr::mutate(
    start = pos,
    end = pos,
    short_species_id = sub("nxOsc([^0-9]+).+", "\\1", assembly),
    multispecies_sequence = paste0(short_species_id, "_", Sequence)
  ) %>%
  dplyr::filter(!is.na(pos), feature == "telomere_seq_split", type == "internal") %>%
  dplyr::select(assembly, Sequence, multispecies_sequence, pos, start, end, short_species_id)

# Remove breaksites that are less than 100 Kb apart (consolidation)
break_sites_per_grs <- break_sites %>%
  dplyr::group_by(multispecies_sequence) %>%
  dplyr::arrange(pos) %>%
  dplyr::mutate(block = cumsum(c(1, diff(pos) > 100e3))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(multispecies_sequence, block) %>%
  dplyr::slice_max(pos) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(assembly, Sequence, pos)

# ------------------------------------------------------------------------------
# 3. Multispecies Oxford Plot Function
# ------------------------------------------------------------------------------
multispecies_oxford_plot <- function(busco_tbl, seq_sizes_tbl, break_sites_tbl, sp_x, sp_y, xlab, ylab, nigon_dict) {
  
  multispecies_index_x <- seq_sizes_tbl %>% dplyr::filter(assembly %in% sp_x) %>% sizes_to_multi_sp_index()
  multispecies_index_y <- seq_sizes_tbl %>% dplyr::filter(assembly %in% sp_y) %>% sizes_to_multi_sp_index()
  
  x_breaks <- busco_subset_and_recoord(break_sites_tbl, multispecies_index_x) 
  y_breaks <- busco_subset_and_recoord(break_sites_tbl, multispecies_index_y)
  
  x_busco <- busco_subset_and_recoord(busco_tbl, multispecies_index_x) 
  y_busco <- busco_subset_and_recoord(busco_tbl, multispecies_index_y)
  
  multispecies_buscos <- dplyr::inner_join(x_busco, y_busco, by = "Busco_id") %>%
    dplyr::select(Busco_id, multispecies_sequence.x, multispecies_start.x, multispecies_end.x,
                  multispecies_sequence.y, multispecies_start.y, multispecies_end.y) %>%
    dplyr::inner_join(nigon_dict, by = c("Busco_id" = "Orthogroup")) %>%
    dplyr::filter(!is.na(nigon))
  
  ggplot(multispecies_buscos, aes(x = multispecies_start.x, y = multispecies_start.y, colour = nigon)) +
    geom_point(alpha = 0.3, size = 0.5) +
    scale_colour_manual(values = nigon_cols) +
    scale_x_continuous(breaks = multispecies_index_x$multispecies_index + (multispecies_index_x$size)/2,
                       labels = sub(".+_", "", multispecies_index_x$multispecies_sequence),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = multispecies_index_y$multispecies_index + (multispecies_index_y$size)/2,
                       labels = sub(".+_", "", multispecies_index_y$multispecies_sequence),
                       expand = c(0, 0)) +
    geom_vline(xintercept = multispecies_index_x$multispecies_index, linewidth = 0.2) +
    geom_hline(yintercept = multispecies_index_y$multispecies_index, linewidth = 0.2) +
    geom_vline(xintercept = x_breaks$multispecies_start, linetype = "dotted", alpha = 0.8) +
    geom_hline(yintercept = y_breaks$multispecies_start, linetype = "dotted", alpha = 0.8) +
    xlab(xlab) +
    ylab(ylab) +
    guides(colour = guide_legend("Nigon\nunit", override.aes = list(size = 3, alpha = 1))) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(angle = 90, hjust = 0.5),
      axis.ticks = element_blank(),
      axis.title = element_text(face = "italic"),
      text = element_text(family = "Helvetica")
    )
}

# ------------------------------------------------------------------------------
# 4. Generate the Plot
# ------------------------------------------------------------------------------
cat("Generating Figure 6 (O. dolichura vs O. sp. JU1382)...\n")

label_func <- function(tolid, dic_ids) {
  dic_ids$public_id[match(tolid, dic_ids$assembly)]
}

x_species <- "nxOscDolc1.1"
y_species <- "nxOscSpeu1.1"
x_lab <- label_func(x_species, tol_to_publicIds)
y_lab <- label_func(y_species, tol_to_publicIds)

figure6 <- multispecies_oxford_plot(
  buscos, seq_sizes, break_sites_per_grs,
  x_species, y_species, x_lab, y_lab, nigonDict
)

# ------------------------------------------------------------------------------
# 5. Save Outputs
# ------------------------------------------------------------------------------
output_pdf <- "figures/main/Figure6.pdf"
output_png <- "figures/main/Figure6.png"

cat("Saving Figure 6 to", output_pdf, "...\n")
ggsave(output_pdf, figure6, units = "in", width = 4.5, height = 4)

cat("Saving Figure 6 to", output_png, "...\n")
ggsave(output_png, figure6, units = "in", width = 4.5, height = 4, dpi = 300)

cat("Figure 6 successfully generated!\n")
