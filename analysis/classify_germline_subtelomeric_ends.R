#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
})

# 1. Species list: 5 species in total
# 4 species from Table S4 (dolichura, onirici, DF5120, JU1382) with 1 rDNA end excluded (11 * 4 = 44 ends)
# O. tipulae with all 12 ends evaluated (12 ends)
species_list <- c("nxOscTipu1.1", "nxOscOnir1.2", "nxOscSper1.1", "nxOscSpeu1.1", "nxOscDolc1.1")

# Strain names for readability
strain_map <- c(
  "nxOscTipu1.1" = "O. tipulae (CEW1)",
  "nxOscOnir1.2" = "O. onirici (PS2068)",
  "nxOscSper1.1" = "O. sp. DF5120 (DF5120)",
  "nxOscSpeu1.1" = "O. sp. JU1382 (JU1382)",
  "nxOscDolc1.1" = "O. dolichura (PS1017)"
)

# 2. Known rDNA-bearing ends (from Table S4 & RNAmmer)
# - O. dolichura: Chr IV Right (ptg000001l, coords 12453593-12464839)
# - O. onirici: Chr II Left (ptg000014l, coords 666-1197170)
# - O. sp. DF5120: Chr II Right (ptg000003l, coords 8139838-8176750)
# - O. sp. JU1382: Chr II Right (ptg000020l, coords 9071364-9093655)
# - O. tipulae: none excluded (all 12 evaluated)
rDNA_flags <- tribble(
  ~species, ~chrom, ~end, ~rDNA_notes,
  "nxOscDolc1.1", "IV", "R", "Table S4 (ptg000001l, coords 12453593-12464839)",
  "nxOscOnir1.2", "II", "L", "Table S4 (ptg000014l, coords 666-1197170)",
  "nxOscSper1.1", "II", "R", "Table S4 (ptg000003l, coords 8139838-8176750)",
  "nxOscSpeu1.1", "II", "R", "Table S4 (ptg000020l, coords 9071364-9093655)"
)

# 3. File paths
rnammer_files <- list.files("data/rnammer", pattern = "rnammer\\.gff3$", full.names = TRUE)
names(rnammer_files) <- sub(".+/([^/]+)\\.rnammer\\.gff3", "\\1", rnammer_files)

clipped_files <- list.files("data/curated_GRS", pattern = "clippedTeloPos\\.tsv$", full.names = TRUE)
names(clipped_files) <- sub(".+/([^/]+)\\.clippedTeloPos\\.tsv", "\\1", clipped_files)

results <- list()

for (sp in species_list) {
  message(paste("\n================== Processing", sp, "(", strain_map[sp], ") =================="))
  
  # Load FASTA
  fa_file <- file.path("data/sequence_sizes", paste0(sp, ".primary.fa.gz"))
  if (!file.exists(fa_file)) {
    warning("FASTA not found: ", fa_file, ". Skipping...")
    next
  }
  fa <- readDNAStringSet(fa_file)
  
  # Load clipped reads
  clipped <- read_tsv(clipped_files[[sp]], 
                      col_names = c("chr", "diminution_pos", "telo_orient", "mapq", 
                                    "telomere_pos", "n_telomeric_reps", "clipped_len", "readLen"), 
                      show_col_types = FALSE)
  
  chrs <- c("I", "II", "III", "IV", "V", "X")
  
  for (chr_name in chrs) {
    fa_idx <- which(names(fa) == chr_name)
    if (length(fa_idx) == 0) {
      fa_idx <- grep(paste0("^", chr_name, " |^", chr_name, "$"), names(fa))
    }
    if (length(fa_idx) == 0) {
      warning("Chromosome ", chr_name, " not found in FASTA for ", sp)
      next
    }
    
    seq_chr <- fa[[fa_idx[1]]]
    chrom_len <- length(seq_chr)
    
    # -------------------------------------------------------------
    # Left End (L)
    # -------------------------------------------------------------
    is_rDNA_L <- nrow(rDNA_flags %>% filter(species == sp, chrom == chr_name, end == "L")) > 0
    rDNA_note_L <- ifelse(is_rDNA_L, (rDNA_flags %>% filter(species == sp, chrom == chr_name, end == "L"))$rDNA_notes[1], NA_character_)
    
    somatic_pos_L <- NA_real_
    somatic_supp_L <- NA_integer_
    germline_pos_L <- NA_real_
    germline_supp_L <- NA_integer_
    class_L <- NA_character_
    telo_len_L <- 0
    non_telo_grs_L <- NA_real_
    
    if (is_rDNA_L) {
      class_L <- "rDNA-Bearing"
    } else {
      # Left reads
      l_reads <- clipped %>% filter(chr == chr_name, telo_orient == "L", diminution_pos < chrom_len / 2)
      
      # Somatic scission: OUTWARD from center (chrom_len/2 -> 1)
      l_soma_clusters <- l_reads %>%
        arrange(desc(diminution_pos)) %>%
        mutate(cluster = cumsum(c(1, abs(diff(diminution_pos)) > 1000))) %>%
        group_by(cluster) %>%
        summarize(pos = round(median(diminution_pos)), n = n(), .groups = "drop") %>%
        filter(n >= 15) %>%
        arrange(desc(pos))
      
      if (nrow(l_soma_clusters) > 0) {
        somatic_pos_L <- l_soma_clusters$pos[1]
        somatic_supp_L <- l_soma_clusters$n[1]
      }
      
      # Germline reads: INWARD from 1
      max_search_L <- ifelse(!is.na(somatic_pos_L), somatic_pos_L - 1000, chrom_len / 4)
      l_germ_clusters <- l_reads %>%
        filter(diminution_pos < max_search_L) %>%
        arrange(diminution_pos) %>%
        mutate(cluster = cumsum(c(1, diff(diminution_pos) > 1000))) %>%
        group_by(cluster) %>%
        summarize(pos = round(median(diminution_pos)), n = n(), .groups = "drop") %>%
        filter(n >= 2) %>%
        arrange(pos)
      
      if (nrow(l_germ_clusters) > 0) {
        germline_pos_L <- l_germ_clusters$pos[1]
        germline_supp_L <- l_germ_clusters$n[1]
      }
      
      # Inspect FASTA terminus (first 2,000 bp)
      term_seq_L <- subseq(seq_chr, 1, min(2000, chrom_len))
      c_fwd_L <- countPattern("TTAGGC", term_seq_L)
      c_rev_L <- countPattern("GCCTAA", term_seq_L)
      fasta_complete_L <- (c_fwd_L + c_rev_L) >= 3
      
      if (fasta_complete_L) {
        class_L <- "Complete"
        telo_len_L <- (c_fwd_L + c_rev_L) * 6
        if (!is.na(somatic_pos_L)) {
          non_telo_grs_L <- max(0, somatic_pos_L - 1 - telo_len_L)
        }
      } else if (!is.na(germline_pos_L)) {
        class_L <- "Bounded Truncation"
        if (!is.na(somatic_pos_L)) {
          non_telo_grs_L <- max(0, somatic_pos_L - germline_pos_L)
        }
      } else {
        class_L <- "Unbounded Truncation"
      }
    }
    
    results[[length(results) + 1]] <- tibble(
      species = sp, strain = strain_map[sp], chrom = chr_name, end = "L", chrom_length = chrom_len,
      is_rDNA = is_rDNA_L, rDNA_notes = rDNA_note_L,
      somatic_scission_pos = somatic_pos_L, somatic_read_support = somatic_supp_L,
      germline_telo_pos = germline_pos_L, germline_read_support = germline_supp_L,
      telo_length_bp = telo_len_L, non_telo_grs_bp = non_telo_grs_L, classification = class_L
    )
    
    # -------------------------------------------------------------
    # Right End (R)
    # -------------------------------------------------------------
    is_rDNA_R <- nrow(rDNA_flags %>% filter(species == sp, chrom == chr_name, end == "R")) > 0
    rDNA_note_R <- ifelse(is_rDNA_R, (rDNA_flags %>% filter(species == sp, chrom == chr_name, end == "R"))$rDNA_notes[1], NA_character_)
    
    somatic_pos_R <- NA_real_
    somatic_supp_R <- NA_integer_
    germline_pos_R <- NA_real_
    germline_supp_R <- NA_integer_
    class_R <- NA_character_
    telo_len_R <- 0
    non_telo_grs_R <- NA_real_
    
    if (is_rDNA_R) {
      class_R <- "rDNA-Bearing"
    } else {
      # Right reads
      r_reads <- clipped %>% filter(chr == chr_name, telo_orient == "R", diminution_pos > chrom_len / 2)
      
      # Somatic scission: OUTWARD from center (chrom_len/2 -> chrom_len)
      r_soma_clusters <- r_reads %>%
        arrange(diminution_pos) %>%
        mutate(cluster = cumsum(c(1, diff(diminution_pos) > 1000))) %>%
        group_by(cluster) %>%
        summarize(pos = round(median(diminution_pos)), n = n(), .groups = "drop") %>%
        filter(n >= 15) %>%
        arrange(pos)
      
      if (nrow(r_soma_clusters) > 0) {
        somatic_pos_R <- r_soma_clusters$pos[1]
        somatic_supp_R <- r_soma_clusters$n[1]
      }
      
      # Germline reads: INWARD from chrom_len
      min_search_R <- ifelse(!is.na(somatic_pos_R), somatic_pos_R + 1000, chrom_len * 3 / 4)
      r_germ_clusters <- r_reads %>%
        filter(diminution_pos > min_search_R) %>%
        arrange(desc(diminution_pos)) %>%
        mutate(cluster = cumsum(c(1, abs(diff(diminution_pos)) > 1000))) %>%
        group_by(cluster) %>%
        summarize(pos = round(median(diminution_pos)), n = n(), .groups = "drop") %>%
        filter(n >= 2) %>%
        arrange(desc(pos))
      
      if (nrow(r_germ_clusters) > 0) {
        germline_pos_R <- r_germ_clusters$pos[1]
        germline_supp_R <- r_germ_clusters$n[1]
      }
      
      # Inspect FASTA terminus (last 2,000 bp)
      term_seq_R <- subseq(seq_chr, max(1, chrom_len - 2000), chrom_len)
      c_fwd_R <- countPattern("TTAGGC", term_seq_R)
      c_rev_R <- countPattern("GCCTAA", term_seq_R)
      fasta_complete_R <- (c_fwd_R + c_rev_R) >= 3
      
      if (fasta_complete_R) {
        class_R <- "Complete"
        telo_len_R <- (c_fwd_R + c_rev_R) * 6
        if (!is.na(somatic_pos_R)) {
          non_telo_grs_R <- max(0, chrom_len - somatic_pos_R - telo_len_R)
        }
      } else if (!is.na(germline_pos_R)) {
        class_R <- "Bounded Truncation"
        if (!is.na(somatic_pos_R)) {
          non_telo_grs_R <- max(0, germline_pos_R - somatic_pos_R)
        }
      } else {
        class_R <- "Unbounded Truncation"
      }
    }
    
    results[[length(results) + 1]] <- tibble(
      species = sp, strain = strain_map[sp], chrom = chr_name, end = "R", chrom_length = chrom_len,
      is_rDNA = is_rDNA_R, rDNA_notes = rDNA_note_R,
      somatic_scission_pos = somatic_pos_R, somatic_read_support = somatic_supp_R,
      germline_telo_pos = germline_pos_R, germline_read_support = germline_supp_R,
      telo_length_bp = telo_len_R, non_telo_grs_bp = non_telo_grs_R, classification = class_R
    )
  }
}

df <- bind_rows(results) %>%
  mutate(total_assembled_grs_bp = ifelse(is.na(non_telo_grs_bp), NA_real_, telo_length_bp + non_telo_grs_bp))

# Save full detailed table
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write_tsv(df, "results/tables/Table_germline_subtelomeric_ends_detailed.tsv")

# Summary for the 4 species (44 ends)
summary_4sp <- df %>%
  filter(!is_rDNA, species %in% c("nxOscDolc1.1", "nxOscOnir1.2", "nxOscSper1.1", "nxOscSpeu1.1")) %>%
  group_by(strain, classification) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = classification, values_from = n, values_fill = 0) %>%
  mutate(Total_Non_rDNA_Ends = rowSums(select(., -strain)))

write_tsv(summary_4sp, "results/tables/Table_germline_subtelomeric_ends_summary_4sp.tsv")

# Summary for all 5 species (56 ends)
summary_5sp <- df %>%
  filter(!is_rDNA) %>%
  group_by(strain, classification) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = classification, values_from = n, values_fill = 0) %>%
  mutate(Total_Non_rDNA_Ends = rowSums(select(., -strain)))

write_tsv(summary_5sp, "results/tables/Table_germline_subtelomeric_ends_summary_all5sp.tsv")

message("\n================== SUMMARY: 4 SPECIES (44 ENDS) ==================")
print(as.data.frame(summary_4sp))

message("\n================== SUMMARY: ALL 5 SPECIES (56 ENDS) ==================")
print(as.data.frame(summary_5sp))

