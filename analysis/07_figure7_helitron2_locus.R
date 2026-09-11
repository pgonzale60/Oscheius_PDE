#!/usr/bin/env Rscript
# ==============================================================================
# Figure 7: A helitron2-like locus associated with eliminated sub-telomeric DNA in Oscheius
# Manuscript: "A helitron2-like locus associated with eliminated sub-telomeric DNA in Oscheius"
# Panels:
#   A: Phylogenetic tree of Helitron2-like helicases
#   B: RNA-seq expression (TPM) of the helitron2-like locus across developmental stages
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(magick)
  library(GenomicRanges)
})

# ------------------------------------------------------------------------------
# 1. Load Utilities & Prepare Data
# ------------------------------------------------------------------------------
source("R/plot_theme.R")

cat("Loading and processing RNA-seq metadata...\n")
plot_stages <- c("2 cell, embryo",
                 "2-4 cell, embryo",
                 "4-8 cell, embryo",
                 "8-16 cell, embryo",
                 "16-32 cell, embryo",
                 "32-64 cell, embryo",
                 "64-128 cell, embryo",
                 "100-150 cell, embryo",
                 "150-200 cell, embryo",
                 "Mixed embryos, embryo",
                 "L1 Larvae, whole worm",
                 "L2 Larvae, whole worm",
                 "Dauer, whole worm",
                 "Post-dauer L3, whole worm",
                 "L3 Larvae, whole worm",
                 "Post-dauer L4, whole worm",
                 "L4 Larvae, whole worm",
                 "Post-dauer young adult, whole worm",
                 "Young Adult, whole worm",
                 "Post-dauer mature adult, whole worm",
                 "Mature Adult, whole worm",
                 "Male Adult, whole worm",
                 "Mixed worms, whole worm") %>%
  gsub(" ", "_", .)

sample_data <- read_csv("data/metadata/RNA-seq/samplesheet.csv", show_col_types = FALSE) %>%
  dplyr::mutate(stage = sub(", biol rep [0-9]", "", sample_title),
                stage = gsub(" ", "_", stage)) %>%
  dplyr::select(names = sample, stage, sample_title) %>%
  dplyr::filter(stage %in% plot_stages)

# ------------------------------------------------------------------------------
# 2. Extract gene coordinates & Identify Helitron2-like locus
# ------------------------------------------------------------------------------
cat("Loading gene coordinates and GRS boundaries...\n")
all_expression <- read_tsv("data/metadata/RNA-seq/otipulae_hisat2_stringtie_gene_expression.tsv.gz",
                           col_names = c("ID", "Name", "Reference", "Strand", "Start",
                                         "End", "Coverage", "FPKM", "TPM", "SRA"),
                           show_col_types = FALSE)

gene_coords <- dplyr::select(all_expression, ID:End) %>%
  dplyr::filter(!duplicated(ID))

# Sequence sizes for Gnm GRanges
fai_files <- list.files("data/sequence_sizes/", full.names = TRUE, pattern = "\\.fai$")
names(fai_files) <- sub(".+/([^/]+)\\.primary\\.fa\\.gz\\.fai$", "\\1", fai_files)
seq_sizes <- map_df(fai_files, read_tsv, col_names = c("Sequence", "size", "L1", "L2", "L3"), col_types = "ciiii", show_col_types = FALSE, .id = "assembly") %>%
  dplyr::filter(assembly == "nxOscTipu1.1")

gnm_gr <- Seqinfo(seqnames = seq_sizes$Sequence, seqlengths = seq_sizes$size, isCircular = rep(FALSE, nrow(seq_sizes)), genome = "nxOscTipu1.1")

genes_gr <- dplyr::rename(gene_coords, gene_id = ID, chrom = Reference) %>%
  dplyr::select(-Name) %>%
  makeGRangesFromDataFrame(seqinfo = gnm_gr, keep.extra.columns = TRUE)

GRS <- read_tsv("data/curated_GRS/nxOscTipu1.1.GRS.bed", col_names = c("Sequence", "start", "end"), show_col_types = FALSE) %>%
  dplyr::rename(chrom = Sequence)

GRS_gr <- makeGRangesFromDataFrame(GRS, seqinfo = gnm_gr)

dimi_genes <- genes_gr$gene_id[genes_gr %over% GRS_gr]
helx_ids <- dimi_genes[grepl("g200", dimi_genes)]

# ------------------------------------------------------------------------------
# 3. RNA-seq Expression (TPM) Profile
# ------------------------------------------------------------------------------
cat("Processing RNA-seq expression for Helitron2-like locus...\n")
salmon_expression <- read_tsv("data/gene/salmon/salmon.merged.transcript_tpm.tsv", show_col_types = FALSE) %>%
  pivot_longer(!c(tx, gene_id), names_to = "sample", values_to = "tpm") %>%
  dplyr::group_by(gene_id, sample) %>%
  dplyr::summarise(gTPM = sum(tpm), .groups = "drop") %>%
  dplyr::inner_join(sample_data, by = c("sample" = "names")) %>%
  dplyr::mutate(stage = factor(stage, levels = plot_stages))

helx <- salmon_expression %>%
  dplyr::filter(gene_id %in% helx_ids)

helx_plot <- helx %>%
  dplyr::mutate(stage = factor(gsub("_", " ", stage), levels = gsub("_", " ", rev(levels(stage))))) %>%
  ggplot(aes(y = stage, x = gTPM + 0.5, shape = gene_id)) +
  scale_x_log10(labels = scales::comma) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  theme_classic(base_size = 11) +
  ylab("Developmental Stage") +
  xlab("Expression (TPM + 0.5)") +
  theme(
    legend.position = "top",
    text = element_text(family = "Helvetica")
  ) +
  guides(shape = guide_legend(title = "Gene ID", nrow = 2))

# ------------------------------------------------------------------------------
# 4. Load Phylogenetic Tree Image (Panel A)
# ------------------------------------------------------------------------------
cat("Loading Pre-rendered Phylogenetic Tree Schematic (Panel A)...\n")
helx_phylo <- ggdraw() + draw_image("data/images/hel2_like_in_multimb_seqs_tree.png", scale = 0.95)

# ------------------------------------------------------------------------------
# 5. Assemble Composite Figure 7
# ------------------------------------------------------------------------------
cat("Composing Figure 7...\n")
top_half <- align_plots(helx_phylo, helx_plot, align = "h", axis = "b")

figure7 <- plot_grid(
  top_half[[1]], top_half[[2]],
  labels = c("A", "B"),
  rel_widths = c(1, 0.9),
  nrow = 1
)

# ------------------------------------------------------------------------------
# 6. Save Outputs
# ------------------------------------------------------------------------------
output_pdf <- "figures/main/Figure7.pdf"
output_png <- "figures/main/Figure7.png"

cat("Saving Figure 7 to", output_pdf, "...\n")
ggsave(output_pdf, figure7, units = "in", width = 12, height = 7)

cat("Saving Figure 7 to", output_png, "...\n")
ggsave(output_png, figure7, units = "in", width = 12, height = 7, dpi = 300)

cat("Figure 7 successfully generated!\n")
