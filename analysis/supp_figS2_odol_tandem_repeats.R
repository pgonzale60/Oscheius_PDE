#!/usr/bin/env Rscript
# ==============================================================================
# Figure S2: Positions of tandem repeats in the Oscheius dolichura genome
# Manuscript: "Positions of tandem repeats in the Oscheius dolichura genome"
# Panel: Binned coverage plot of tandem repeat density across chromosomes
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(scales)
})

# ------------------------------------------------------------------------------
# 1. Load Utilities & Prepare Data
# ------------------------------------------------------------------------------
source("R/plot_theme.R")

cat("Loading sequence sizes...\n")
fai_files <- list.files("data/sequence_sizes/", full.names = TRUE, pattern = "\\.fai$")
names(fai_files) <- sub(".+/([^/]+)\\.primary\\.fa\\.gz\\.fai$", "\\1", fai_files)

seq_sizes <- map_df(fai_files, read_tsv,
                    col_names = c("Sequence", "size", "L1", "L2", "L3"),
                    col_types = "ciiii", show_col_types = FALSE, .id = "assembly") %>%
  dplyr::filter(assembly == "nxOscDolc1.1", !grepl("MT", Sequence)) %>%
  dplyr::mutate(multispecies_sequence = paste0(assembly, "_", Sequence))

cat("Loading break sites (diminution coordinates)...\n")
break_files <- list.files("data/curated_GRS/", full.names = TRUE, pattern = "nxOscDolc1.1\\.chr_diminutions_sites\\.tsv$")
names(break_files) <- "nxOscDolc1.1"

break_sites <- map_df(break_files, read_tsv, show_col_types = FALSE, .id = "assembly") %>%
  dplyr::rename(Sequence = chr, pos = diminution_pos) %>%
  dplyr::filter(!is.na(pos), feature == "telomere_seq_split", type == "internal") %>%
  dplyr::mutate(multispecies_sequence = paste0(assembly, "_", Sequence)) %>%
  dplyr::select(multispecies_sequence, diminution_pos = pos)

cat("Loading TRF output...\n")
trf <- read_tsv("data/repeats/nxOscDolc1.1.trf.tsv", 
                col_names = c("chrom", "start", "end", "unit_size", "num_copies", 
                              "perc_ident", "perc_indels", "score", "entropy", "sequence"),
                show_col_types = FALSE) %>%
  dplyr::mutate(
    assembly = "nxOscDolc1.1",
    id = 1:n(),
    chrom_full = paste0(assembly, "_", chrom)
  ) %>%
  dplyr::filter(!grepl("MT", chrom))

# ------------------------------------------------------------------------------
# 2. Build GRanges and Calculate Binned Coverage
# ------------------------------------------------------------------------------
cat("Calculating tandem repeat genomic density (binned coverage)...\n")
gnm_gr <- Seqinfo(seqnames = seq_sizes$multispecies_sequence, seqlengths = seq_sizes$size, isCircular = rep(FALSE, nrow(seq_sizes)), genome = "nxOscDolc1.1")

trf_gr <- trf %>%
  dplyr::filter(abs(end - start) > 0) %>%
  dplyr::rename(seqnames = chrom_full) %>%
  makeGRangesFromDataFrame(seqinfo = gnm_gr, keep.extra.columns = TRUE, seqnames.field = "seqnames")

# Remove overlapping redundancies by reducing ranges, then calculate per-base coverage
gr1_cov <- GenomicRanges::coverage(GenomicRanges::reduce(trf_gr))

window_size <- 50000
osch_bins <- tileGenome(gnm_gr, tilewidth = window_size, cut.last.tile.in.chrom = TRUE)

trf_binned <- GenomicRanges::binnedAverage(osch_bins, gr1_cov, "binned_cov") %>%
  as_tibble() %>%
  dplyr::bind_rows(break_sites %>% dplyr::rename(seqnames = multispecies_sequence)) %>%
  dplyr::mutate(chr = sub("[^_]+_(.+)", "\\1", seqnames))

# ------------------------------------------------------------------------------
# 3. Generate the Plot
# ------------------------------------------------------------------------------
cat("Generating Figure S2 (Tandem Repeat Density)...\n")

fig_S2 <- trf_binned %>%
  ggplot(aes(x = start, y = binned_cov)) + 
  facet_grid(chr ~ .) +
  geom_bar(position = "stack", stat = "identity", width = window_size, fill = "#0072B2") +
  geom_vline(aes(xintercept = diminution_pos), linetype = "dotted", color = "#D55E00", linewidth = 0.8) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(breaks = scales::pretty_breaks(3), limits = c(0, 1)) +
  labs(
    x = "Genomic Position", 
    y = "Tandem Repeat Fraction"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "grey80", fill = NA),
    text = element_text(family = "Helvetica")
  )

# ------------------------------------------------------------------------------
# 4. Save Outputs
# ------------------------------------------------------------------------------
output_pdf <- "figures/supplementary/FigureS2.pdf"
output_png <- "figures/supplementary/FigureS2.png"

cat("Saving Figure S2 to", output_pdf, "...\n")
ggsave(output_pdf, fig_S2, units = "in", width = 8, height = 6)

cat("Saving Figure S2 to", output_png, "...\n")
ggsave(output_png, fig_S2, units = "in", width = 8, height = 6, dpi = 300)

cat("Figure S2 successfully generated!\n")
