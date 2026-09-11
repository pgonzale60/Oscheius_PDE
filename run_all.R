#!/usr/bin/env Rscript
# ==============================================================================
# Master Execution Script for Oscheius_PDE
# Reproduces all Figures and Tables for the manuscript:
# "Programmed DNA elimination in Oscheius nematodes"
# ==============================================================================

cat("=================================================================\n")
cat("Starting Oscheius PDE Pipeline...\n")
cat("=================================================================\n\n")

# Ensure directories exist
dir.create("figures/main", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/supplementary", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Function to run a script and check its status
run_script <- function(script_path, description) {
  cat(sprintf("\n---> Running %s: %s\n", description, script_path))
  
  start_time <- Sys.time()
  status <- system2("Rscript", args = script_path)
  end_time <- Sys.time()
  
  if (status != 0) {
    stop(sprintf("\n[ERROR] Script %s failed with exit code %d.", script_path, status))
  } else {
    cat(sprintf("\n[SUCCESS] %s completed in %.1f seconds.\n", description, as.numeric(difftime(end_time, start_time, units = "secs"))))
  }
}

# Run Main Figures
run_script("analysis/01_figure1_assemblies_and_PDE_sites.R", "Figure 1")
run_script("analysis/02_figure2_features_of_PDE.R", "Figure 2")
run_script("analysis/03_figure3_hic_chromatin_structure.R", "Figure 3")
run_script("analysis/04_figure4_break_motifs.R", "Figure 4")
run_script("analysis/05_figure5_odol_variable_elimination.R", "Figure 5")
run_script("analysis/06_figure6_synteny_and_rearrangements.R", "Figure 6")
run_script("analysis/07_figure7_helitron2_locus.R", "Figure 7")

# Run Supplementary Figures
run_script("analysis/supp_figS1_grs_gene_sharing.R", "Figure S1")
run_script("analysis/supp_figS2_odol_tandem_repeats.R", "Figure S2")
run_script("analysis/supp_figS4_oxford_plots_multispecies.R", "Figure S4")
run_script("analysis/supp_figS5_chrX_rearrangements.R", "Figure S5")

# Run Tables
run_script("analysis/generate_all_tables.R", "Tables 1-3 & S1-S8")

cat("\n=================================================================\n")
cat("Pipeline completed successfully!\n")
cat("All figures are available in: figures/\n")
cat("All tables are available in: results/tables/\n")
cat("=================================================================\n")
