# R/plot_theme.R
# Shared color palettes, factor levels, and theme helpers for Oscheius_PDE

suppressPackageStartupMessages({
  library(ggplot2)
  library(tibble)
})

# Nigon element colors
nigon_cols <- c(
  "A" = "#af0e2b",
  "B" = "#e4501e",
  "C" = "#4caae5",
  "D" = "#f3ac2c",
  "E" = "#57b741",
  "N" = "#8880be",
  "X" = "#81008b",
  "-" = "#aaaaaa"
)

# Strain to public species name dictionary
strain_to_publicIds <- tibble(
  Strain = c("PS1017", "PS2068", "DF5120", "JU1382", "CEW1", "caenorhabditis_elegans"),
  public_id = c("O. dolichura", "O. onirici", "O. sp. DF5120",
                "O. sp. JU1382", "O. tipulae", "C. elegans")
)

# Tree of Life assembly ID to public species name dictionary
tol_to_publicIds <- tibble(
  assembly = c("nxOscDolc1.1", "nxOscOnir1.2", "nxOscSper1.1", 
              "nxOscSpeu1.1", "nxOscTipu1.1", "caenorhabditis_elegans"),
  short_species_id = c("Dolc", "Onir", "Sper", "Speu", "Tipu", "Cele"),
  short_id = c("OscDolc1.1", "OscOnir1.2", "OscSper1.1", "OscSpeu1.1", "OscTipu1.1", "Cele"),
  public_id = c("O. dolichura", "O. onirici", "O. sp. DF5120",
                "O. sp. JU1382", "O. tipulae", "C. elegans")
)

# Standard factor levels for consistent species ordering
species_order <- c("O. dolichura", "O. sp. DF5120", "O. sp. JU1382", "O. onirici", "O. tipulae")
species_order_rev <- rev(species_order)

# Minimal publication theme
theme_pde <- function(base_size = 11, base_family = "Helvetica") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}
