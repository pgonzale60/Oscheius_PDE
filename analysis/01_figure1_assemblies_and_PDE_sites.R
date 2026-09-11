#!/usr/bin/env Rscript
# ==============================================================================
# Figure 1: Chromosomal assemblies and sites of programmed DNA elimination
#           in Oscheius species
# Manuscript: "New high quality reference genomes for four Oscheius species"
# Panels:
#   A: Species phylogeny with PDE scission characteristics (number & type of GRS)
#   B: Chromosomal Nigon element painting with PDE scission sites marked
#   C: O. tipulae telomeric FISH microscopy image
#   D: O. dolichura telomeric FISH microscopy image
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(cowplot)
  library(ape)
  library(ggtree)
  library(treeio)
  library(tidytree)
})

# ------------------------------------------------------------------------------
# 1. Load Utilities & Metadata
# ------------------------------------------------------------------------------

source("R/plot_theme.R")
source("R/busco_utils.R")

# Load Nigon element dictionary
nigonDict <- read_tsv("data/metadata/gene2Nigon_busco20200927.tsv.gz",
                      col_types = "cc", show_col_types = FALSE)

# ------------------------------------------------------------------------------
# 2. Load and Prepare BUSCO Data
# ------------------------------------------------------------------------------

cat("Loading BUSCO tables...\n")
buscofiles <- list.files("data/busco/", full.names = TRUE, pattern = "nematoda_odb10_full_table\\.tsv$")
names(buscofiles) <- sub(".+/([^/]+)\\.primary_nematoda_odb10_full_table\\.tsv$", "\\1", buscofiles)

buscos <- map_df(buscofiles, read_busco, .id = "assembly") %>%
  standardize_chromosomes() %>%
  mutate(
    short_species_id = sub("nxOsc([^0-9]+).+", "\\1", assembly),
    multispecies_sequence = paste0(short_species_id, "_", Sequence)
  )

# ------------------------------------------------------------------------------
# 3. Load Curated GRS Coordinates
# ------------------------------------------------------------------------------

cat("Loading GRS bed files...\n")
grsfiles <- list.files("data/curated_GRS/", full.names = TRUE, pattern = "\\.GRS\\.bed$")
names(grsfiles) <- sub(".+/([^/]+)\\.GRS\\.bed$", "\\1", grsfiles)

grs_coords <- map_df(grsfiles, read_tsv,
                     col_names = c("Sequence", "start", "end", "has_GRS"),
                     col_types = "ciic",
                     .id = "assembly") %>%
  mutate(
    short_species_id = sub("nxOsc([^0-9]+).+", "\\1", assembly),
    multispecies_sequence = paste0(short_species_id, "_", Sequence),
    size = ifelse(grepl("no", has_GRS), NA, end - start + 1),
    GRD = (start + end) / 2
  )

# ------------------------------------------------------------------------------
# 4. Panel B: Chromosomal Nigon Painting with Scission Sites
# ------------------------------------------------------------------------------

cat("Building Panel B: Chromosome painting with Nigon elements & PDE sites...\n")

tm_g <- filter(grs_coords, has_GRS == "yes")

b_for_plot <- bind_rows(tm_g, buscos) %>%
  left_join(tol_to_publicIds, by = "short_species_id") %>%
  mutate(
    public_id = factor(public_id, levels = species_order_rev)
  ) %>%
  filter(!grepl("unloc", Sequence), !is.na(public_id)) %>%
  left_join(nigonDict, by = c("Busco_id" = "Orthogroup"))

nigon_plot <- ggplot(b_for_plot) +
  facet_grid(public_id ~ Sequence, scales = "free") +
  geom_rect(
    aes(xmin = start - 2e4, xmax = start + 2e4, ymax = 0, ymin = 1, fill = nigon),
    alpha = 0.8
  ) +
  scale_fill_manual(values = nigon_cols[1:7]) +
  geom_vline(aes(xintercept = GRD), linetype = "dotted", alpha = 1, linewidth = 0.8) +
  xlab("Position (Mb)") +
  scale_x_continuous(
    labels = function(x) x / 1e6,
    expand = c(0.005, 1),
    breaks = scales::pretty_breaks(2)
  ) +
  scale_y_continuous(breaks = NULL) +
  theme(
    strip.text.y = element_blank(),
    strip.text.x = element_text(margin = margin(0, 0, 0, 0, "cm"), face = "bold"),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend("Nigon element", nrow = 1), color = "none")

# ------------------------------------------------------------------------------
# 5. Panel A: Species Phylogeny with PDE Scission Numbers
# ------------------------------------------------------------------------------

cat("Building Panel A: Species phylogeny with PDE features...\n")

treetext <- "(A._rhodensis:0.0916469,
((O._dolichura:0.0635095[&&NHX:nGRS=29],
((O._onirici:0.0235031[&&NHX:nGRS=9],
O._tipulae:0.0229182)100:0.0317448,
O._sp._JU1382:0.0613351)55:0.00585547)51:0.00672745,
O._sp._DF5120:0.057538)1:0.0916469);"

osch_tree <- read.tree(textConnection(treetext))
osch_tree_o <- drop.tip(osch_tree, "A._rhodensis")

dimi_tbl <- tibble(
  label = c("O._dolichura", "O._onirici", "O._tipulae", "O._sp._JU1382", "O._sp._DF5120"),
  nGRS = c("11[1]+7", "10[2]", "12", "11[1]+10", "11[1]+2"),
  type = c("internal and subtelomeric", "subtelomeric only",
           "subtelomeric only", "internal and subtelomeric",
           "internal and subtelomeric")
)

tree_plot <- ggtree(osch_tree_o) %<+% dimi_tbl +
  geom_tiplab(as_ylab = TRUE) +
  xlim(c(0, 0.082)) +
  geom_label(aes(label = nGRS, fill = type), size = 3) +
  scale_fill_discrete(breaks = c("subtelomeric only", "internal and subtelomeric")) +
  guides(fill = guide_legend(nrow = 2, na.translate = FALSE, title = element_blank())) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(family = "Helvetica", face = "italic", color = "black", size = rel(1.1)),
    text = element_text(family = "Helvetica")
  )

# ------------------------------------------------------------------------------
# 6. Panels C & D: Telomeric FISH Microscopy Images
# ------------------------------------------------------------------------------

cat("Loading Panels C & D (FISH microscopy images)...\n")

otipu_fish_path <- "data/images/Otipulae_telomeric_FISH_cropped.jpg"
odoli_fish_path <- "data/images/Odolichura_telomeric_FISH_cropped.jpg"

otipu_fish <- ggdraw() + draw_image(otipu_fish_path)
odoli_fish <- ggdraw() + draw_image(odoli_fish_path)

# ------------------------------------------------------------------------------
# 7. Assemble Composite Figure 1
# ------------------------------------------------------------------------------

cat("Composing Figure 1...\n")

top_half <- align_plots(tree_plot, nigon_plot, otipu_fish, odoli_fish,
                        align = "v", axis = "l")
names(top_half) <- c("fig1a", "fig1b", "fig1c", "fig1d")

first_row <- plot_grid(
  top_half$fig1a, top_half$fig1b,
  labels = c("A", "B"),
  label_x = c(0, -0.02),
  rel_widths = c(1, 2),
  nrow = 1
)

fig1_bottom_panel <- plot_grid(
  otipu_fish, odoli_fish,
  ncol = 2,
  rel_widths = c(1, 1),
  labels = c("C", "D")
)

figure1 <- plot_grid(
  first_row, fig1_bottom_panel,
  ncol = 1
)

# ------------------------------------------------------------------------------
# 8. Save Outputs
# ------------------------------------------------------------------------------

output_pdf <- "figures/main/Figure1.pdf"
output_png <- "figures/main/Figure1.png"

cat("Saving Figure 1 to", output_pdf, "...\n")
ggsave(output_pdf, plot = figure1, width = 12, height = 7)

cat("Saving Figure 1 to", output_png, "...\n")
ggsave(output_png, plot = figure1, width = 12, height = 7, dpi = 300)

cat("Figure 1 successfully generated!\n")
