#!/usr/bin/env Rscript
# ==============================================================================
# Figure 2: Features of PDE in Oscheius species
# Manuscript: "Programmed DNA elimination is present in all five Oscheius species"
# Panels:
#   A: JBrowse visualization of scission & telomere addition in O. sp. DF5120
#   B: JBrowse visualization of scission & telomere addition in O. dolichura
#   C: Span of Germline-Restricted Sequences (GRS, distal vs internal scissions)
#   D: Telomere length distributions (germline ends vs somatic de novo additions)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(ggpubr)
  library(cowplot)
})

# ------------------------------------------------------------------------------
# 1. Load Utilities & Metadata Dictionaries
# ------------------------------------------------------------------------------

source("R/plot_theme.R")

find_mode <- function(x) {
  u <- unique(x)
  tab <- tabulate(match(x, u))
  tu <- u[tab == max(tab)]
  tu[1]
}

species_levels <- rev(c("O. sp. DF5120", "O. dolichura",
                        "O. sp. JU1382", "O. onirici",
                        "O. tipulae"))

# ------------------------------------------------------------------------------
# 3. Load Data
# ------------------------------------------------------------------------------

cat("Loading GRS summary table...\n")
grs <- read_tsv("data/curated_GRS/all_sp_gsheet_20221212.tsv", show_col_types = FALSE) %>%
  mutate(type = ifelse(Is_germline_end == "Yes", "distal", "internal")) %>%
  left_join(strain_to_publicIds, by = "Strain") %>%
  mutate(public_id = factor(public_id, levels = species_levels))

cat("Loading classified telomere positions...\n")
telocoord_per_pos_clas <- read_tsv(
  "data/curated_GRS/telo_position_grouped_and_classified_20221222.tsv",
  show_col_types = FALSE
) %>%
  left_join(tol_to_publicIds, by = "assembly") %>%
  mutate(
    public_id = factor(public_id, levels = species_levels),
    multispecies_sequence = paste0(assembly, "_", Sequence)
  ) %>%
  filter(lineage %in% c("soma", "germline"))

cat("Loading individual clipped telomere read coordinates...\n")
telocoord_files <- list.files("data/curated_GRS/", full.names = TRUE, pattern = "\\.clippedTeloPos\\.tsv$")
names(telocoord_files) <- sub(".+/([^/]+)\\.clippedTeloPos\\.tsv$", "\\1", telocoord_files)

telocoord <- map_df(telocoord_files, read_tsv,
                    col_names = c("chr", "diminution_pos", "telo_orient", "mapq",
                                  "telomere_pos", "n_telomeric_reps",
                                  "clipped_len", "readLen"),
                    col_types = "ciciiiii",
                    .id = "assembly") %>%
  rename(Sequence = chr) %>%
  mutate(
    multispecies_sequence = paste0(assembly, "_", Sequence),
    mapped_len = readLen - clipped_len
  ) %>%
  filter(
    mapq > 0,
    n_telomeric_reps > 20,
    mapped_len > 1000,
    !grepl("MT", Sequence)
  )

# ------------------------------------------------------------------------------
# 4. Data Wrangling
# ------------------------------------------------------------------------------

cat("Wrangling telomere coordinates and clusters...\n")
telocoord_per_pos <- telocoord %>%
  mutate(
    corrected_pos = diminution_pos + telomere_pos,
    break_id = paste0(multispecies_sequence, "_", diminution_pos)
  ) %>%
  group_by(multispecies_sequence, telo_orient) %>%
  arrange(diminution_pos) %>%
  mutate(
    block = ifelse(
      assembly != "nxOscTipu1.1" | telo_orient == "R",
      cumsum(c(1, diff(diminution_pos) > 10)),
      cumsum(c(1, diff(diminution_pos) > max(clipped_len) / 10))
    )
  ) %>%
  ungroup() %>%
  group_by(multispecies_sequence, telo_orient, block) %>%
  mutate(modal_corrected_pos = find_mode(corrected_pos)) %>%
  ungroup() %>%
  left_join(
    select(telocoord_per_pos_clas, multispecies_sequence, modal_corrected_pos, public_id, lineage),
    by = c("multispecies_sequence", "modal_corrected_pos")
  ) %>%
  filter(
    n_telomeric_reps > 20,
    !is.na(lineage),
    mapq > 30
  )

tb_n <- telocoord_per_pos %>%
  group_by(public_id, lineage) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = lineage, values_from = n, values_fill = 0, id_expand = TRUE) %>%
  pivot_longer(names_to = "lineage", values_to = "n", cols = c("germline", "soma"))

# ------------------------------------------------------------------------------
# 5. Generate Panels
# ------------------------------------------------------------------------------

cat("Building Panel C: GRS size plot...\n")
size_plot <- ggplot(grs, aes(x = public_id, y = Span)) +
  geom_jitter(width = 0.2, height = 0.2, size = 0.4, alpha = 0.7) +
  scale_y_continuous(
    labels = function(x) x / 1e3,
    expand = c(0.005, 8000),
    breaks = scales::pretty_breaks(5)
  ) +
  facet_grid(~ type) +
  ylab("Span (Kb)") +
  xlab("") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, face = "italic", hjust = 1, vjust = 0.5),
    strip.text = element_text(face = "bold")
  )

cat("Building Panel D: Telomere length plot...\n")
telo_len_plot <- ggplot(telocoord_per_pos, aes(x = public_id, y = n_telomeric_reps * 6)) +
  geom_jitter(width = 0.25, alpha = 0.6, size = 0.3) +
  geom_text(
    data = tb_n,
    aes(y = 26000, label = n),
    size = 3,
    angle = 90,
    hjust = 1
  ) +
  scale_y_continuous(
    labels = function(x) x / 1e3,
    breaks = scales::pretty_breaks(5)
  ) +
  facet_grid(~ lineage) +
  ylab("Telomere length (Kb)") +
  xlab("") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, face = "italic", hjust = 1, vjust = 0.5),
    strip.text = element_text(face = "bold")
  )

cat("Loading Panels A & B (JBrowse screenshot images)...\n")
df5120_img_path <- "data/images/DF5120_short_X_break.jpg"
odoli_img_path <- "data/images/Odoli_minimal_IV.jpg"

DF5120_break <- ggdraw() + draw_image(df5120_img_path)
odoli_break <- ggdraw() + draw_image(odoli_img_path)

# ------------------------------------------------------------------------------
# 6. Assemble Composite Figure 2
# ------------------------------------------------------------------------------

cat("Composing Figure 2...\n")
first_row <- plot_grid(
  DF5120_break, size_plot,
  labels = c("A", "C"),
  rel_widths = c(2, 1.3),
  nrow = 1
)

second_row <- plot_grid(
  odoli_break, telo_len_plot,
  labels = c("B", "D"),
  rel_widths = c(2, 1.3),
  nrow = 1
)

figure2 <- plot_grid(
  first_row, second_row,
  ncol = 1
)

# ------------------------------------------------------------------------------
# 7. Save Outputs
# ------------------------------------------------------------------------------

dir.create("figures/main", recursive = TRUE, showWarnings = FALSE)

output_pdf <- "figures/main/Figure2.pdf"
output_png <- "figures/main/Figure2.png"

cat("Saving Figure 2 to", output_pdf, "...\n")
ggsave(output_pdf, figure2, units = "in", width = 10, height = 6.5)

cat("Saving Figure 2 to", output_png, "...\n")
ggsave(output_png, figure2, units = "in", width = 10, height = 6.5, dpi = 300)

cat("Figure 2 successfully generated!\n")
