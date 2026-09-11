#!/usr/bin/env Rscript
# ==============================================================================
# Figure S4: Oxford plots comparing Oscheius species genomes
# Manuscript: "Oxford plots comparing Oscheius species genomes"
# Panel: A-D Synteny dot plots of Oscheius sp. DF5120 vs other species
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

# ------------------------------------------------------------------------------
# 1. Load Utilities & Prepare Data
# ------------------------------------------------------------------------------
source("R/plot_theme.R")
source("R/busco_utils.R")
source("R/synteny_utils.R")

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

cat("Loading BUSCO tables...\n")
buscofiles <- list.files("data/busco/", full.names = TRUE, pattern = "nematoda_odb10_full_table\\.tsv$")
names(buscofiles) <- sub(".+/([^/]+)\\.primary_nematoda_odb10_full_table\\.tsv$", "\\1", buscofiles)

buscos <- map_df(buscofiles, read_busco, .id = "assembly") %>%
  standardize_chromosomes() %>%
  dplyr::mutate(
    short_species_id = sub("nxOsc([^0-9]+).+", "\\1", assembly),
    multispecies_sequence = paste0(short_species_id, "_", Sequence)
  )

cat("Loading GRS break sites...\n")
break_files <- list.files("data/curated_GRS/", full.names = TRUE, pattern = "\\.chr_diminutions_sites\\.tsv$")
names(break_files) <- sub(".+/([^/]+)\\.chr_diminutions_sites\\.tsv$", "\\1", break_files)

break_sites <- map_df(break_files, read_tsv, show_col_types = FALSE, .id = "assembly") %>%
  dplyr::rename(Sequence = chr, pos = diminution_pos) %>%
  dplyr::mutate(
    start = pos,
    end = pos,
    short_species_id = sub("nxOsc([^0-9]+).+", "\\1", assembly),
    multispecies_sequence = paste0(short_species_id, "_", Sequence)
  )

cat("Loading Nigon dictionary...\n")
nigonDict <- read_tsv("data/metadata/gene2Nigon_busco20200927.tsv.gz", col_types = "cc", show_col_types = FALSE)

multispecies_oxford_plot <- function(busco_tbl, seq_sizes_tbl, break_sites_tbl, sp_x, sp_y, xlab = "", ylab = "", nigon_dict) {
  
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
# 2. Generate Oxford Plots
# ------------------------------------------------------------------------------
cat("Generating Oxford plots against O. sp. DF5120 (nxOscSper1.1)...\n")
y_species <- "nxOscSper1.1"

pDo <- multispecies_oxford_plot(buscos, seq_sizes, break_sites, "nxOscDolc1.1", y_species, "O. dolichura", "O. sp. DF5120", nigonDict) + theme(legend.position = "none")
pJU <- multispecies_oxford_plot(buscos, seq_sizes, break_sites, "nxOscSpeu1.1", y_species, "O. sp. JU1382", "O. sp. DF5120", nigonDict) + theme(legend.position = "none")
pOn <- multispecies_oxford_plot(buscos, seq_sizes, break_sites, "nxOscOnir1.2", y_species, "O. onirici", "O. sp. DF5120", nigonDict) + theme(legend.position = "none")
pTi <- multispecies_oxford_plot(buscos, seq_sizes, break_sites, "nxOscTipu1.1", y_species, "O. tipulae", "O. sp. DF5120", nigonDict) + theme(legend.position = "none")

combinedplot <- cowplot::plot_grid(pDo, pJU, pOn, pTi, labels = c("A", "B", "C", "D"), ncol = 2)

# Extract legend from one of the plots
legend_plot <- multispecies_oxford_plot(buscos, seq_sizes, break_sites, "nxOscDolc1.1", y_species, "", "", nigonDict) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(title = "Nigon Element", nrow = 1))
legend <- cowplot::get_legend(legend_plot)

fig_S4 <- cowplot::plot_grid(combinedplot, legend, ncol = 1, rel_heights = c(20, 1))

# ------------------------------------------------------------------------------
# 3. Save Outputs
# ------------------------------------------------------------------------------
output_pdf <- "figures/supplementary/FigureS4.pdf"
output_png <- "figures/supplementary/FigureS4.png"

cat("Saving Figure S4 to", output_pdf, "...\n")
ggsave(output_pdf, fig_S4, units = "in", width = 8, height = 8.4)

cat("Saving Figure S4 to", output_png, "...\n")
ggsave(output_png, fig_S4, units = "in", width = 8, height = 8.4, dpi = 300)

cat("Figure S4 successfully generated!\n")
