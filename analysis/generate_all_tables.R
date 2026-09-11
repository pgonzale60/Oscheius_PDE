#!/usr/bin/env Rscript
# ==============================================================================
# Generate Tables
# Consolidates and formats tables into results/tables/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

cat("Generating Tables...\n")
dir.create("results/tables/", recursive = TRUE, showWarnings = FALSE)

# Table 1: Assembly statistics
cat("Generating Table 1: Assembly Statistics\n")
assm_stats <- read_delim("data/tables/oscheius_assmstats.tsv", 
                         delim = " ", col_names = c("Strain", "Size_Mb", "Contigs", "N50_Mb", "Type"), 
                         show_col_types = FALSE)
write_tsv(assm_stats, "results/tables/Table1_Assembly_Stats.tsv")

# Table 2: Sequencing read statistics
cat("Generating Table 2: Sequencing Read Statistics\n")
read_stats <- read_tsv("data/tables/bothBatch_filtered_stats.tsv", show_col_types = FALSE)
write_tsv(read_stats, "results/tables/Table2_Sequencing_Reads.tsv")

# Table S4: Telomere positions and coverage
cat("Generating Table S4: Telomere position and coverage\n")
table_s4 <- read_tsv("data/tables/supp_table_s4.tsv", show_col_types = FALSE)
write_tsv(table_s4, "results/tables/TableS4_Telomere_Coverage.tsv")

# Table S8: Orthogroups and Gene Content
cat("Generating Table S8: Orthogroups in Eliminated Regions\n")
table_s8 <- read_tsv("data/tables/GRS_genes_and_their_interpro_and_celeg_hit_20220904.tsv", show_col_types = FALSE)
write_tsv(table_s8, "results/tables/TableS8_Orthogroups_or_accessions.tsv")

cat("All tables generated successfully in results/tables/!\n")
