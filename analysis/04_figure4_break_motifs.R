#!/usr/bin/env Rscript
# ==============================================================================
# Figure 4: Motifs identified at the break sites of PDE in Oscheius species
# Manuscript: "Sequence motifs identified at PDE breakage sites"
# Panels:
#   A: Primary sequence logos of conserved cleavage motifs across five Oscheius species
#   B: Secondary cleavage motif identified in O. sp. DF5120
#   C: Motif match score distribution at break sites vs GRS vs core genome
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(cowplot)
  library(ggseqlogo)
  library(motifStack)
  library(GenomicRanges)
})

# ------------------------------------------------------------------------------
# 1. Load Utilities & Species Metadata
# ------------------------------------------------------------------------------

source("R/plot_theme.R")

phylo_order <- c("OscTipu1.1", "OscOnir1.2", "OscSpeu1.1", "OscSper1.1", "OscDolc1.1")

tol_to_motif_ids <- tibble(
  short_id = c("OscDolc1.1", "OscOnir1.2", "OscSper1.1", "OscSpeu1.1", "OscTipu1.1"),
  public_id = c("O. dolichura", "O. onirici", "O. sp. DF5120", "O. sp. JU1382", "O. tipulae")
)

# Custom color scheme for nucleotide logos
cs1 <- make_col_scheme(
  chars = c("A", "C", "G", "T"),
  cols = c("#009E73", "#0072B2", "#E69F00", "#D55E00")
)

# ------------------------------------------------------------------------------
# 2. Panels A & B: Motif Extraction & Alignment (MEME)
# ------------------------------------------------------------------------------

cat("Importing and aligning MEME motifs...\n")

species_meme_files <- c(
  "data/motifs/meme/nxOscDolc1.1.meme.txt",
  "data/motifs/meme/nxOscOnir1.2.meme.txt",
  "data/motifs/meme/nxOscSper1.1.meme.txt",
  "data/motifs/meme/nxOscSpeu1.1.meme.txt",
  "data/motifs/meme/nxOscTipu1.1.meme.txt"
)

all_motifs <- importMatrix(species_meme_files, format = "meme")

# Extract secondary motif of DF5120 (4th motif in file list: DF5120 has 2 motifs)
pcm_secondary_df5120 <- all_motifs[[4]]
pcm_secondary_df5120$name <- "O. sp. DF5120"
pcm_secondary_df5120$color <- colorset(colorScheme = "blindnessSafe")

secondaryMotif <- list()
secondaryMotif$`O. sp. DF5120` <- pcm_secondary_df5120@mat

# Primary motifs (excluding secondary motif at index 4)
primary_motifs <- all_motifs[-4]
strains <- sub(".+/nx([a-zA-Z0-9.]+)\\.meme\\.txt$", "\\1", species_meme_files)
names(primary_motifs) <- strains

ord_motifs <- list()
for (strain_num in seq_along(primary_motifs)) {
  strain <- strains[strain_num]
  primary_motifs[[strain]]$name <- tol_to_motif_ids$public_id[match(strain, tol_to_motif_ids$short_id)]
  primary_motifs[[strain]]$color <- colorset(colorScheme = "blindnessSafe")
  ord_motifs[[match(strain, phylo_order)]] <- trimMotif(primary_motifs[[strain]], t = 0.4)
}

# Align PFMs
pfmsAligned <- DNAmotifAlignment(ord_motifs, rcpostfix = "")

motif_matrices <- list()
for (i in 1:5) {
  motif_matrices[[pfmsAligned[[i]]$name]] <- pfmsAligned[[i]]@mat
}

cat("Building Panels A & B: Motif logos...\n")
main_motif_plot <- ggseqlogo(motif_matrices, ncol = 1, col_scheme = cs1) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "italic"),
    text = element_text(family = "Helvetica")
  )

sec_motif_plot <- ggseqlogo(secondaryMotif, col_scheme = cs1) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "italic"),
    text = element_text(family = "Helvetica")
  )

# ------------------------------------------------------------------------------
# 3. Panel C: Motif Specificity (FIMO Score Distribution)
# ------------------------------------------------------------------------------

cat("Loading break sites, sequence sizes, and FIMO genomic occurrences...\n")

break_sites <- read_tsv("data/curated_GRS/all_sp_teloclipped_based_20220919.tsv", show_col_types = FALSE) %>%
  dplyr::filter(for_coord != "no", feature == "telomere_seq_split") %>%
  dplyr::mutate(chr = sub("chr_", "", chr)) %>%
  dplyr::rename(Sequence = chr) %>%
  dplyr::arrange(assembly, Sequence, diminution_pos)

telocoord_files <- list.files("data/curated_GRS/", full.names = TRUE, pattern = "\\.clippedTeloPos\\.tsv$")
names(telocoord_files) <- sub(".+/([^/]+)\\.clippedTeloPos\\.tsv$", "\\1", telocoord_files)

telocoord <- map_df(telocoord_files, read_tsv,
                    col_names = c("chr", "diminution_pos", "telo_orient", "mapq"),
                    col_types = "cici",
                    show_col_types = FALSE,
                    .id = "assembly") %>%
  dplyr::rename(Sequence = chr) %>%
  dplyr::filter(mapq > 30, !grepl("MT", Sequence)) %>%
  dplyr::mutate(multispecies_sequence = paste0(assembly, "_", Sequence))

GRS_files <- list.files("data/curated_GRS/", full.names = TRUE, pattern = "nx.*\\.GRS\\.bed$")
names(GRS_files) <- sub(".+/([^/]+)\\.GRS\\.bed$", "\\1", GRS_files)

GRS <- map_df(GRS_files, read_tsv,
              col_names = c("Sequence", "start", "end", "has_GRS"),
              col_types = "ciic",
              show_col_types = FALSE,
              .id = "assembly") %>%
  dplyr::filter(!grepl("no", has_GRS)) %>%
  dplyr::mutate(multispecies_sequence = paste0(assembly, "_", Sequence))

fai_files <- list.files("data/sequence_sizes/", full.names = TRUE, pattern = "\\.fai$")
names(fai_files) <- sub(".+/([^/]+)\\.primary\\.fa\\.gz\\.fai$", "\\1", fai_files)

seq_sizes <- map_df(fai_files, read_tsv,
                    col_names = c("Sequence", "size", "L1", "L2", "L3"),
                    col_types = "ciiii",
                    show_col_types = FALSE,
                    .id = "assembly") %>%
  dplyr::select(assembly, Sequence, size) %>%
  dplyr::mutate(multispecies_sequence = paste0(assembly, "_", Sequence))

fimofiles <- list.files("data/motifs/fimo/", full.names = TRUE, pattern = "fimo\\.tsv$")
names(fimofiles) <- sub(".+/([^/]+)\\.fimo\\.tsv$", "\\1", fimofiles)

fimo <- map_df(fimofiles, read_tsv,
               comment = "#",
               show_col_types = FALSE,
               .id = "assembly") %>%
  dplyr::filter(!grepl("MT", sequence_name)) %>%
  dplyr::mutate(
    telo_only = grepl("TTAGGCTTAGGCTTAGGCTTAGGCTT", matched_sequence, fixed = TRUE),
    motif_id = ifelse(is.na(motif_id), "N", motif_id),
    motif_alt_id = ifelse(is.na(motif_alt_id), "MEME-1", motif_alt_id)
  )

# GenomicRanges intersection
cat("Intersecting genomic occurrences with scission and GRS boundaries...\n")
gnm_gr <- Seqinfo(
  seqnames = seq_sizes$multispecies_sequence,
  seqlengths = seq_sizes$size,
  isCircular = rep(FALSE, nrow(seq_sizes)),
  genome = "Oscheius"
)

fimo_gr <- fimo %>%
  dplyr::mutate(chrom = paste0(assembly, "_", sequence_name), fimo_ID = paste0(chrom, start)) %>%
  dplyr::select(chrom, start, end = stop, score = score, pval = `p-value`, qval = `q-value`, fimo_ID, motif_alt_id, telo_only) %>%
  makeGRangesFromDataFrame(seqinfo = gnm_gr, keep.extra.columns = TRUE)

GRS_gr <- GRS %>%
  dplyr::rename(chrom = multispecies_sequence) %>%
  makeGRangesFromDataFrame(seqinfo = gnm_gr)

break_sites_gr <- break_sites %>%
  dplyr::mutate(chrom = paste0(assembly, "_", Sequence), end = diminution_pos + 1) %>%
  dplyr::select(chrom, start = diminution_pos, end, precise_support, type) %>%
  makeGRangesFromDataFrame(seqinfo = gnm_gr, keep.extra.columns = TRUE)

# Select highest scoring non-overlapping ranges
fimo_gr <- fimo_gr[order(fimo_gr$score, decreasing = TRUE)]
fimo_m1_gr <- fimo_gr[fimo_gr$motif_alt_id == "MEME-1"]
idx_m1 <- findOverlaps(fimo_m1_gr, fimo_m1_gr) %>%
  as_tibble() %>%
  arrange(queryHits, subjectHits) %>%
  filter(!duplicated(queryHits)) %>%
  pull(subjectHits) %>%
  unique()
highest_m1_fimo_gr <- fimo_m1_gr[idx_m1]

fimo_m2_gr <- fimo_gr[fimo_gr$motif_alt_id == "MEME-2"]
idx_m2 <- findOverlaps(fimo_m2_gr, fimo_m2_gr) %>%
  as_tibble() %>%
  arrange(queryHits, subjectHits) %>%
  filter(!duplicated(queryHits)) %>%
  pull(subjectHits) %>%
  unique()
highest_m2_fimo_gr <- fimo_m2_gr[idx_m2]

highest_m2_fimo_gr_f15 <- flank(highest_m2_fimo_gr, 15)
OscSper_merge <- mergeByOverlaps(highest_m1_fimo_gr, highest_m2_fimo_gr_f15)
highest_fimo_gr <- c(
  highest_m1_fimo_gr,
  highest_m2_fimo_gr[!highest_m2_fimo_gr$fimo_ID %in% OscSper_merge$highest_m2_fimo_gr_f15$fimo_ID]
)

highest_fimo_gr$score <- ifelse(
  !highest_fimo_gr$fimo_ID %in% OscSper_merge$highest_m1_fimo_gr$fimo_ID,
  highest_fimo_gr$score,
  OscSper_merge$highest_m1_fimo_gr$score[match(highest_fimo_gr$fimo_ID, OscSper_merge$highest_m1_fimo_gr$fimo_ID)] +
    OscSper_merge$highest_m2_fimo_gr_f15$score[match(highest_fimo_gr$fimo_ID, OscSper_merge$highest_m1_fimo_gr$fimo_ID)]
)

highest_fimo_gr$atBreakSite <- highest_fimo_gr %over% break_sites_gr
highest_fimo_gr$atBreakSite <- ifelse(
  highest_fimo_gr$atBreakSite,
  "breaksite",
  ifelse(highest_fimo_gr %over% GRS_gr, "GRS", "core")
)

pd_fimo <- as_tibble(highest_fimo_gr) %>%
  mutate(short_id = sub("nxOsc(.{4}).+", "\\1", seqnames)) %>%
  left_join(
    tibble(
      short_id = c("Dolc", "Onir", "Sper", "Speu", "Tipu"),
      public_id = c("O. dolichura", "O. onirici", "O. sp. DF5120", "O. sp. JU1382", "O. tipulae")
    ),
    by = "short_id"
  ) %>%
  mutate(
    public_id = factor(public_id, levels = c("O. dolichura", "O. sp. DF5120", "O. sp. JU1382", "O. onirici", "O. tipulae")),
    atBreakSite = factor(atBreakSite, levels = c("breaksite", "GRS", "core"))
  ) %>%
  filter(score > 0)

tb_n <- pd_fimo %>%
  group_by(public_id, atBreakSite) %>%
  summarise(n = n(), .groups = "drop")

cat("Building Panel C: Motif score specificity across genomic compartments...\n")
motif_specificity <- ggplot(pd_fimo, aes(y = score, x = atBreakSite)) +
  facet_grid(. ~ public_id, scales = "free") +
  geom_jitter(width = 0.25, size = 0.3, alpha = 0.6) +
  geom_text(data = tb_n, aes(y = 48, label = n), size = 3, angle = 90) +
  scale_y_continuous(expand = c(0.005, 5), breaks = scales::pretty_breaks(5)) +
  ylab("Motif match score") +
  xlab("") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.text = element_text(face = "italic"),
    text = element_text(family = "Helvetica")
  )

# ------------------------------------------------------------------------------
# 4. Assemble Composite Figure 4
# ------------------------------------------------------------------------------

cat("Composing Figure 4...\n")
top_half <- align_plots(main_motif_plot, sec_motif_plot, motif_specificity, align = "v", axis = "l")
names(top_half) <- c("fig4a", "fig4b", "fig4c")

first_col <- plot_grid(
  top_half$fig4a,
  labels = c("A"),
  nrow = 1
)

second_col <- plot_grid(
  sec_motif_plot, top_half$fig4c,
  labels = c("B", "C"),
  rel_widths = c(1.3, 2),
  rel_heights = c(1.2, 4),
  nrow = 2
)

figure4 <- plot_grid(
  first_col, second_col,
  rel_widths = c(1.2, 1),
  ncol = 2
)

# ------------------------------------------------------------------------------
# 5. Save Outputs
# ------------------------------------------------------------------------------

output_pdf <- "figures/main/Figure4.pdf"
output_png <- "figures/main/Figure4.png"

cat("Saving Figure 4 to", output_pdf, "...\n")
ggsave(output_pdf, figure4, units = "in", width = 11, height = 6.5)

cat("Saving Figure 4 to", output_png, "...\n")
ggsave(output_png, figure4, units = "in", width = 11, height = 6.5, dpi = 300)

cat("Figure 4 successfully generated!\n")
