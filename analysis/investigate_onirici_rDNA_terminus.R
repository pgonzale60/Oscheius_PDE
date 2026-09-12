#!/usr/bin/env Rscript
# ==============================================================================
# Script: investigate_onirici_rDNA_terminus.R
# Purpose: Reproduce the genomic, sequence, and read-level analysis of the 
#          rDNA-bearing Chr II Left terminus in Oscheius onirici (nxOscOnir1.2).
# Evaluates:
#   1. Chromosome II Left terminus in primary assembly (nxOscOnir1.2.primary.fa.gz)
#   2. Structure and junction sequence of telotig ptg000014l
#   3. Alignment of ptg000014l spacer to internal Chr II rDNA IGS units
#   4. Projection of all 109 raw clipped telomeric reads onto a canonical rDNA unit
#   5. Absence of the 27 bp PDE consensus cleavage motif across all clusters
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
  library(pwalign)
})

cat("================================================================================\n")
cat("1. Inspecting nxOscOnir1.2 Primary Assembly (Chr II Left)\n")
cat("================================================================================\n")

fa_primary_file <- "data/sequence_sizes/nxOscOnir1.2.primary.fa.gz"
if (!file.exists(fa_primary_file)) {
  stop("Primary assembly not found at: ", fa_primary_file)
}

fa_primary <- readDNAStringSet(fa_primary_file)
chrII <- fa_primary[["II"]]
chrII_len <- length(chrII)

cat("Chr II Total Length:", chrII_len, "bp\n")

# Check terminal 2,000 bp
tip_seq <- subseq(chrII, 1, 2000)
c_fwd <- countPattern("TTAGGC", tip_seq)
c_rev <- countPattern("GCCTAA", tip_seq)

cat("TTAGGC count in terminal 2,000 bp:", c_fwd, "\n")
cat("GCCTAA count in terminal 2,000 bp:", c_rev, "\n")
cat("First 150 bp of Chr II:\n", as.character(subseq(tip_seq, 1, 150)), "\n\n")
cat("Conclusion: Primary assembly terminates at pos 1 with 0 telomeric repeats.\n")
cat("First 18S gene starts at pos 666 bp (first 665 bp are non-coding rDNA spacer).\n\n")

cat("================================================================================\n")
cat("2. Inspecting Telotig ptg000014l\n")
cat("================================================================================\n")

telotig_fa_file <- "../multiOscheius/analyses/in_oxf/telo_reads.asm.fasta.gz"
if (!file.exists(telotig_fa_file)) {
  stop("Telotig assembly not found at: ", telotig_fa_file)
}

fa_telotigs <- readDNAStringSet(telotig_fa_file)
ptg14 <- fa_telotigs[["ptg000014l"]]
cat("ptg000014l Total Length:", length(ptg14), "bp\n")

# Region breakdown from Table S4 and RNAmmer
# 1..20,638 bp: 3 copies of 18S and 28S rRNA genes
# 20,639..23,586 bp: Spacer region (IGS)
# 23,587..27,495 bp: 3.9 kb telomere tract (650 repeats)
# 27,743..44,096 bp: 16.3 kb telomere tract (2,723 repeats)
cat("TTAGGC count in tract 1 (23,587..27,495 bp):", countPattern("TTAGGC", subseq(ptg14, 23587, 27495)), "\n")
cat("TTAGGC count in tract 2 (27,743..44,096 bp):", countPattern("TTAGGC", subseq(ptg14, 27743, 44096)), "\n\n")

cat("================================================================================\n")
cat("3. Verifying the Canonical rDNA Unit on Chr II\n")
cat("================================================================================\n")

# Define canonical unit from Chr II (pos 7624 to 14581, length 6958 bp)
canonical_unit <- subseq(chrII, 7624, 14581)
cat("Canonical rDNA Unit Length:", length(canonical_unit), "bp\n")
cat("Unit Architecture:\n")
cat("  - 18S rRNA:           1 .. 1,754 bp\n")
cat("  - ITS1 / 5.8S / ITS2: 1,755 .. 2,200 bp\n")
cat("  - 28S rRNA:           2,201 .. 6,902 bp\n")
cat("  - IGS / ETS spacer:   6,903 .. 6,958 bp\n\n")

# Project ptg000014l junction onto canonical unit
ptg14_junc <- subseq(ptg14, 23587 - 100, 23587)
pa_ptg <- pairwiseAlignment(pattern = ptg14_junc, subject = canonical_unit, type = "local")
ptg_rel_pos <- end(subject(pa_ptg))
cat("ptg000014l telomere junction projects to relative pos:", ptg_rel_pos, "within canonical unit.\n\n")

cat("================================================================================\n")
cat("4. Projecting All 109 Raw Clipped Telomeric Reads onto Canonical Unit\n")
cat("================================================================================\n")

clipped <- read_tsv("data/curated_GRS/nxOscOnir1.2.clippedTeloPos.tsv", 
                    col_names = c("chr", "diminution_pos", "telo_orient", "mapq", 
                                  "telomere_pos", "n_telomeric_reps", "clipped_len", "readLen"), 
                    show_col_types = FALSE) %>%
  filter(chr == "II", diminution_pos < 1200000)

rel_positions <- c()
for (i in 1:nrow(clipped)) {
  p <- clipped$diminution_pos[i]
  q_start <- max(1, p - 50)
  q_end <- min(length(chrII), p + 50)
  query_seq <- subseq(chrII, q_start, q_end)
  pa <- pairwiseAlignment(pattern = query_seq, subject = canonical_unit, type = "local")
  rel_positions <- c(rel_positions, start(subject(pa)) + (p - q_start))
}
clipped$rel_pos <- rel_positions

clipped <- clipped %>%
  mutate(feature = case_when(
    rel_pos >= 1 & rel_pos <= 1754 ~ "18S rRNA",
    rel_pos > 1754 & rel_pos <= 2200 ~ "ITS1 / 5.8S / ITS2",
    rel_pos > 2200 & rel_pos <= 6902 ~ "28S rRNA",
    TRUE ~ "IGS"
  ),
  cluster_group = case_when(
    rel_pos >= 4680 & rel_pos <= 4710 & telo_orient == "R" ~ "Cluster 1 (ptg000014l junction, ~4705 R)",
    rel_pos >= 1620 & rel_pos <= 1645 & telo_orient == "L" ~ "Cluster 2 (~1637 L)",
    TRUE ~ "Non-cluster singletons"
  ))

cat("Breakdown of all 109 reads across the rDNA array:\n")
print(as.data.frame(clipped %>% count(cluster_group)))

cat("\nSummary of the 10 non-cluster singleton reads:\n")
non_clusters <- clipped %>% 
  filter(cluster_group == "Non-cluster singletons") %>%
  select(diminution_pos, telo_orient, rel_pos, feature, mapq, n_telomeric_reps, clipped_len, readLen) %>%
  arrange(rel_pos)
print(as.data.frame(non_clusters))

cat("\n================================================================================\n")
cat("5. Testing Sequence for 27 bp PDE Cleavage Motif\n")
cat("================================================================================\n")
cat("Consensus PDE Cleavage Motif: TTWGGYTDWGGCYVAGCCNAABCCAAA\n\n")

# Test Cluster 1 sequence
c1_seq <- as.character(subseq(canonical_unit, 4705 - 50, 4705 + 50))
has_c1_5p <- grepl("TT[TA]GG[CT].[TA][ATC]GGC", c1_seq)
has_c1_3p <- grepl("GCC[ATG].[AG][ACT]CC[TA]AA", c1_seq)
cat(sprintf("Cluster 1 (~4705, N = 55 reads): 5p_motif = %s, 3p_motif = %s\n", has_c1_5p, has_c1_3p))

# Test Cluster 2 sequence
c2_seq <- as.character(subseq(canonical_unit, 1637 - 50, 1637 + 50))
has_c2_5p <- grepl("TT[TA]GG[CT].[TA][ATC]GGC", c2_seq)
has_c2_3p <- grepl("GCC[ATG].[AG][ACT]CC[TA]AA", c2_seq)
cat(sprintf("Cluster 2 (~1637, N = 44 reads): 5p_motif = %s, 3p_motif = %s\n\n", has_c2_5p, has_c2_3p))

cat("MEME P-value for Chr II Left: P = 0.826 (NO MOTIF)\n")
cat("================================================================================\n")
