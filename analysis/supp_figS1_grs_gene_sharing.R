#!/usr/bin/env Rscript
# ==============================================================================
# Figure S1: Lack of shared genes in eliminated DNA
# Manuscript: "Lack of shared genes in eliminated DNA"
# Panel: UpSet plot of Orthogroups inside the GRS vs the rest of the genome
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(plyranges)
  library(ggupset)
  library(cowplot)
})

# ------------------------------------------------------------------------------
# 1. Load Utilities & Prepare Data
# ------------------------------------------------------------------------------
source("R/plot_theme.R")

tol_to_publicIds <- tibble(
  assembly = c("nxOscDolc1.1", "nxOscOnir1.2", "nxOscSper1.1", "nxOscSpeu1.1", "nxOscTipu1.1", "caenorhabditis_elegans.PRJNA13758"),
  public_id = c("O. dolichura", "O. onirici", "O. sp. DF5120", "O. sp. JU1382", "O. tipulae", "C. elegans")
)

# Load OrthoFinder groups
cat("Loading OrthoFinder groups...\n")
orthoGenes <- read_tsv("data/gene/orthofinder/Orthogroups.tsv", show_col_types = FALSE)
names(orthoGenes)[-1] <- sub("\\.primary\\.fa\\.gz$", "", names(orthoGenes)[-1])

orthoGenelist <- tibble()
for (sp in c("nxOscDolc1.1", "nxOscOnir1.2", "nxOscSper1.1", "nxOscSpeu1.1", "nxOscTipu1.1", "caenorhabditis_elegans.PRJNA13758")) {
  if (!sp %in% names(orthoGenes)) next
  
  og <- orthoGenes %>% dplyr::filter(!is.na(.data[[sp]]))
  orthoGs <- purrr::map2(dplyr::pull(og, .data[[sp]]), og$Orthogroup, function(x, y) {
    spt_gs <- as.character(unlist(stringr::str_split(x, ", ")))
    data.frame(transcript_id = spt_gs, OG_id = rep(y, length(spt_gs)), stringsAsFactors = FALSE)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(assembly = sp)
  orthoGenelist <- dplyr::bind_rows(orthoGenelist, orthoGs)
}

# ------------------------------------------------------------------------------
# 2. Extract gene coordinates & Identify GRS intersections
# ------------------------------------------------------------------------------
cat("Loading gene coordinates and GRS boundaries...\n")
genecoord_files <- list.files("data/gene/GTFs/", full.names = TRUE, pattern = "\\.gff\\.gz$")
names(genecoord_files) <- sub(".+/combined_(.+)\\.gff\\.gz$", "\\1", genecoord_files)

genecoord <- purrr::map_df(genecoord_files, read_tsv,
                           col_names = c("Sequence", "source", "type", "start", "end", "score", "strand", "phase", "ID"),
                           comment = "#", show_col_types = FALSE, .id = "assembly") %>%
  dplyr::filter(type == "mRNA") %>%
  dplyr::filter(!grepl("MT", Sequence)) %>%
  dplyr::mutate(multispecies_sequence = paste0(assembly, "_", Sequence))

# Load GRS files
GRS_files <- list.files("data/curated_GRS/", full.names = TRUE, pattern = "nx.*\\.GRS\\.bed$")
names(GRS_files) <- sub(".+/([^/]+)\\.GRS\\.bed$", "\\1", GRS_files)

GRS <- purrr::map_df(GRS_files, read_tsv,
                     col_names = c("Sequence", "start", "end"),
                     show_col_types = FALSE, .id = "assembly") %>%
  dplyr::mutate(multispecies_sequence = paste0(assembly, "_", Sequence))

# Sequence sizes
fai_files <- list.files("data/sequence_sizes/", full.names = TRUE, pattern = "\\.fai$")
names(fai_files) <- sub(".+/([^/]+)\\.primary\\.fa\\.gz\\.fai$", "\\1", fai_files)
seq_sizes <- purrr::map_df(fai_files, read_tsv,
                           col_names = c("Sequence", "size", "L1", "L2", "L3"),
                           col_types = "ciiii", show_col_types = FALSE, .id = "assembly") %>%
  dplyr::select(assembly, Sequence, size) %>%
  dplyr::mutate(multispecies_sequence = paste0(assembly, "_", Sequence))

# Build GRanges
gnm_gr <- Seqinfo(seqnames = seq_sizes$multispecies_sequence, seqlengths = seq_sizes$size, isCircular = rep(FALSE, nrow(seq_sizes)), genome = "Oscheius")

genes_gr <- genecoord %>%
  dplyr::rename(chrom = multispecies_sequence) %>%
  dplyr::mutate(transcript_id = paste0(assembly, "_", sub(".*ID=([^;]+).*", "\\1", ID))) %>%
  dplyr::select(-ID, -assembly, -Sequence) %>%
  makeGRangesFromDataFrame(seqinfo = gnm_gr, keep.extra.columns = TRUE)

GRS_gr <- GRS %>%
  dplyr::rename(chrom = multispecies_sequence) %>%
  dplyr::group_by(chrom) %>%
  dplyr::mutate(ID = 1:n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(dimiSpan = end - start, uID = paste(chrom, ID, sep = "_")) %>%
  dplyr::select(-ID, -assembly, -Sequence) %>%
  makeGRangesFromDataFrame(seqinfo = gnm_gr, keep.extra.columns = TRUE)

cat("Intersecting genes with GRS...\n")
GR_genes <- join_overlap_inner(genes_gr, GRS_gr, suffix = c(".g", ".coord")) %>%
  as_tibble() %>%
  dplyr::select(-source, -type, -score, -phase) %>%
  dplyr::left_join(
    dplyr::mutate(orthoGenelist, transcript_id = paste0(assembly, "_", transcript_id)),
    by = "transcript_id"
  )

notGR_genes <- as_tibble(genes_gr) %>%
  dplyr::filter(!transcript_id %in% GR_genes$transcript_id) %>%
  dplyr::left_join(
    dplyr::mutate(orthoGenelist, transcript_id = paste0(assembly, "_", transcript_id)),
    by = "transcript_id"
  )

# ------------------------------------------------------------------------------
# 3. Wrangle for UpSet Plot
# ------------------------------------------------------------------------------
allInvOGs <- unique(GR_genes$OG_id)

allInvGenes <- as_tibble(genes_gr) %>%
  dplyr::select(-source, -type, -score, -phase) %>%
  dplyr::left_join(
    dplyr::mutate(orthoGenelist, transcript_id = paste0(assembly, "_", transcript_id)),
    by = "transcript_id"
  ) %>%
  dplyr::filter(OG_id %in% allInvOGs) %>%
  dplyr::mutate(in_GRS = ifelse(transcript_id %in% GR_genes$transcript_id, "in_GRS", "in_core"))

nonCoreFamGR <- as_tibble(genes_gr) %>%
  dplyr::select(-source, -type, -score, -phase) %>%
  dplyr::left_join(
    dplyr::mutate(orthoGenelist, transcript_id = paste0(assembly, "_", transcript_id)),
    by = "transcript_id"
  ) %>%
  dplyr::filter(!OG_id %in% notGR_genes$OG_id)

fdupType <- dplyr::filter(allInvGenes, in_GRS == "in_GRS") %>%
  dplyr::group_by(assembly, OG_id) %>%
  dplyr::mutate(fake_OG = paste0(OG_id, "_", 1:n())) %>%
  dplyr::ungroup()

nfamMultiOGType <- nonCoreFamGR %>%
  dplyr::group_by(assembly, OG_id) %>%
  dplyr::mutate(fake_OG = paste0(OG_id, "_", 1:n())) %>%
  dplyr::ungroup()

gForPlot <- fdupType %>%
  dplyr::group_by(fake_OG) %>%
  dplyr::summarise(
    gType = list(tol_to_publicIds$public_id[match(assembly, tol_to_publicIds$assembly)]),
    .groups = "drop"
  )

nfamgForPlot <- nfamMultiOGType %>%
  dplyr::group_by(fake_OG) %>%
  dplyr::summarise(
    gType = list(tol_to_publicIds$public_id[match(assembly, tol_to_publicIds$assembly)]),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 4. Generate the Plot
# ------------------------------------------------------------------------------
cat("Generating Figure S1 (UpSet plot of Orthogroups)...\n")

personal_theme <- function() { 
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = "grey20"), 
    panel.grid.major = element_line(colour = "grey92"),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "italic"),
    text = element_text(family = "Helvetica")
  )
}
 
upsGRgns <- ggplot(gForPlot, aes(x = gType)) +
  geom_bar(fill = "#0072B2") +
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5, family = "Helvetica", size = 3) +
  scale_x_upset() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  xlab("") + ylab("Orthogroups inside GRS") +
  personal_theme()

nCoreupsGRgns <- ggplot(nfamgForPlot, aes(x = gType)) +
  geom_bar(fill = "#D55E00") +
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5, family = "Helvetica", size = 3) +
  scale_x_upset() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  xlab("") + ylab("Exclusive Orthogroups (GRS only)") +
  personal_theme()
  
fig_S1 <- cowplot::plot_grid(upsGRgns, nCoreupsGRgns, ncol = 1, rel_heights = c(1, 0.8), labels = c("A", "B"))

# ------------------------------------------------------------------------------
# 5. Save Outputs
# ------------------------------------------------------------------------------
output_pdf <- "figures/supplementary/FigureS1.pdf"
output_png <- "figures/supplementary/FigureS1.png"

cat("Saving Figure S1 to", output_pdf, "...\n")
ggsave(output_pdf, fig_S1, units = "in", width = 5, height = 7)

cat("Saving Figure S1 to", output_png, "...\n")
ggsave(output_png, fig_S1, units = "in", width = 5, height = 7, dpi = 300)

cat("Figure S1 successfully generated!\n")
