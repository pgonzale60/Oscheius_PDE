#!/usr/bin/env Rscript
# ==============================================================================
# Figure 3: Chromatin conformation capture and the structure of Oscheius chromosomes
# Manuscript: "Impact of PDE on genomic architecture"
# Panels:
#   A: Hi-C contact matrix for Oscheius onirici (128 kb resolution)
#   B: Hi-C contact matrix for Oscheius dolichura (128 kb resolution)
#   C: Insulation score (trans-eigenvector 1) profiles across chromosomes with
#      PDE scission sites marked by vertical red dotted lines
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(ggpubr)
  library(cowplot)
})

# ------------------------------------------------------------------------------
# 1. Load Utilities & Metadata
# ------------------------------------------------------------------------------

source("R/plot_theme.R")

# ------------------------------------------------------------------------------
# 2. Load Hi-C Trans Interactions (Cooltools)
# ------------------------------------------------------------------------------

cat("Loading Hi-C trans-eigenvector data...\n")
trans_files <- list.files("data/hic/", full.names = TRUE, pattern = ".*trans\\.vecs\\.tsv\\.gz$")
names(trans_files) <- sub(".+/([^/]+)_eigs-trans\\.trans\\.vecs\\.tsv\\.gz$", "\\1", trans_files)

trans_inter <- map_df(trans_files, read_tsv, show_col_types = FALSE, .id = "assembly") %>%
  mutate(
    short_species_id = sub("nxOsc([^0-9]+).+", "\\1", assembly),
    multispecies_sequence = paste0(short_species_id, "_", chrom)
  ) %>%
  select(assembly, Sequence = chrom, multispecies_sequence, start, end, E1)

# ------------------------------------------------------------------------------
# 3. Load Scission Break Sites
# ------------------------------------------------------------------------------

cat("Loading scission break sites...\n")
break_files <- list.files("data/curated_GRS/", full.names = TRUE, pattern = "chr_diminutions_sites\\.tsv$")
names(break_files) <- sub(".+/([^/]+)\\.chr_diminutions_sites\\.tsv$", "\\1", break_files)

break_sites <- map_df(break_files, read_tsv, show_col_types = FALSE, .id = "assembly") %>%
  rename(Sequence = chr, pos = diminution_pos) %>%
  mutate(
    short_species_id = sub("nxOsc([^0-9]+).+", "\\1", assembly),
    multispecies_sequence = paste0(short_species_id, "_", Sequence)
  ) %>%
  filter(!is.na(pos), feature == "telomere_seq_split") %>%
  select(assembly, Sequence, multispecies_sequence, pos)

# Remove breaksites that are less than 100 Kb apart
break_sites_per_grs <- break_sites %>%
  group_by(multispecies_sequence) %>%
  arrange(pos) %>%
  mutate(block = cumsum(c(1, diff(pos) > 100e3))) %>%
  ungroup() %>%
  group_by(multispecies_sequence, block) %>%
  slice_max(pos) %>%
  ungroup() %>%
  arrange(assembly, Sequence, pos)

# ------------------------------------------------------------------------------
# 4. Data Wrangling & Alignment
# ------------------------------------------------------------------------------

cat("Wrangling insulation scores and scission coordinates...\n")
b_for_plot <- bind_rows(trans_inter, break_sites_per_grs) %>%
  left_join(tol_to_publicIds, by = "assembly") %>%
  mutate(
    public_id = factor(public_id, levels = rev(c("C. elegans", "O. sp. DF5120", "O. dolichura",
                                                 "O. sp. JU1382", "O. onirici"))),
    trans_eigen1 = E1
  ) %>%
  filter(
    !grepl("MT|unloc|scaff", multispecies_sequence),
    (trans_eigen1 > -1.8 & trans_eigen1 < 1.8) | is.na(start),
    assembly != "nxOscTipu1.1"
  ) %>%
  mutate(
    trans_eigen1 = ifelse(assembly %in% c("nxOscDolc1.1", "nxOscOnir1.2"), trans_eigen1, trans_eigen1 * -1)
  )

# ------------------------------------------------------------------------------
# 5. Panel C: Insulation Plot
# ------------------------------------------------------------------------------

cat("Building Panel C: Insulation score profiles...\n")
insulation_plot <- b_for_plot %>%
  ggplot(aes(x = start, y = trans_eigen1, group = multispecies_sequence)) +
  geom_line() +
  facet_grid(public_id ~ Sequence) +
  geom_vline(aes(xintercept = pos), linetype = "dotted", alpha = 0.7, colour = "red") +
  scale_x_continuous(
    labels = function(x) x / 1e6,
    expand = c(0.005, 1),
    breaks = scales::pretty_breaks(2)
  ) +
  ylab("Insulation score") +
  xlab("Position (Mb)") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text.y = element_text(face = "italic"),
    strip.text.x = element_text(face = "bold")
  )

# ------------------------------------------------------------------------------
# 6. Panels A & B: Hi-C Contact Maps
# ------------------------------------------------------------------------------

cat("Loading Panels A & B (Hi-C contact maps)...\n")
onir_hic_path <- "data/hic/O_onirici_hic_def_128k_maxval05.jpg"
odoli_hic_path <- "data/hic/O_dolichura_hic_ICE_128k_maxval05.jpg"

onir_hic <- ggdraw() + draw_image(onir_hic_path)
odoli_hic <- ggdraw() + draw_image(odoli_hic_path)

# ------------------------------------------------------------------------------
# 7. Assemble Composite Figure 3
# ------------------------------------------------------------------------------

cat("Composing Figure 3...\n")
top_half <- align_plots(onir_hic, odoli_hic, insulation_plot, align = "v", axis = "l")
names(top_half) <- c("fig3a", "fig3b", "fig3c")

first_row <- plot_grid(
  top_half$fig3a, top_half$fig3b,
  labels = c("A", "B"),
  nrow = 1
)

second_row <- plot_grid(
  top_half$fig3c,
  labels = c("C"),
  nrow = 1
)

figure3 <- plot_grid(first_row, second_row, ncol = 1, rel_heights = c(1, 1.7))

# ------------------------------------------------------------------------------
# 8. Save Outputs
# ------------------------------------------------------------------------------

output_pdf <- "figures/main/Figure3.pdf"
output_png <- "figures/main/Figure3.png"

cat("Saving Figure 3 to", output_pdf, "...\n")
ggsave(output_pdf, figure3, units = "in", width = 7, height = 8.5)

cat("Saving Figure 3 to", output_png, "...\n")
ggsave(output_png, figure3, units = "in", width = 7, height = 8.5, dpi = 300)

cat("Figure 3 successfully generated!\n")
