#!/usr/bin/env Rscript
# ==============================================================================
# Figure 5: An X chromosome region is variably eliminated in Oscheius dolichura
# Manuscript: "An X chromosome region is variably eliminated in Oscheius dolichura"
# Panels:
#   A: Diagram (loaded as image) showing coverage vectors on X chromosome
#   B: Correlation and distribution plots of the alpha/beta/gamma fragments
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(cowplot)
  library(magick)
})

# ------------------------------------------------------------------------------
# 1. Load Utilities
# ------------------------------------------------------------------------------
source("R/plot_theme.R")

# ------------------------------------------------------------------------------
# 2. Extract statistics / proportions of alpha, beta, gamma fragments
# ------------------------------------------------------------------------------
cat("Processing Table of Coverage Proportions for alpha, beta, gamma...\n")

# Proportions given in the original script
# X = alpha, Y = beta, Z = gamma
X <- c(45.83, 36.44, 31.36, 77.71, 33.58, 63.16)
Y <- c(18.5, 4.0, 1.5, 31.2, 1.4, 20.4)
Z <- c(33.3, 40.5, 29.6, 76.5, 30.7, 61.4)
dset_names <- c("HiFi", "ONT", "Hi-C", "Hi-C adults", "Hi-C larvae", "WGS")

X <- X/100
Y <- Y/100
Z <- Z/100

df <- data.frame(X, Y, Z, sequencing = dset_names)

cat("Generating Scatterplot...\n")
p_scatter <- ggplot(df, aes(x=Y*100, y=Z*100)) + 
  geom_point(color = "black", size = 2) + 
  geom_smooth(data = df %>% filter(sequencing != "HiFi"), method="lm", se=TRUE, color = "#0072B2") +
  geom_text_repel(aes(label=sequencing), box.padding = 0.5, point.padding = 0.5, family = "Helvetica", size = 3) +
  labs(x="% beta", y="% gamma") +
  theme_minimal(base_size = 11) +
  theme(
    text = element_text(family = "Helvetica")
  )

# ------------------------------------------------------------------------------
# 3. Process coverage vectors for the alpha/beta/gamma regions in ONT
# ------------------------------------------------------------------------------
cat("Processing Coverage distributions...\n")

dolcov <- read_tsv("data/coverage/nxOscDolc1.1.multiple.regions.bed.gz",
                   col_names = c("chr", "start", "end", "mdn_cov", "type", "sample"),
                   show_col_types = FALSE)

p_boxplot <- filter(dolcov, type != "nonrep.GRS", mdn_cov < 3000, sample == "ont") %>% 
  mutate(sample = "ONT") %>%
  ggplot(aes(x = type, y = mdn_cov, fill = type)) +
  facet_grid( ~ sample) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  xlab("Region Type") +
  ylab("Median Coverage") +
  theme_minimal(base_size = 11) +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# ------------------------------------------------------------------------------
# 4. Load JBrowse/Diagram Vector (Panel A)
# ------------------------------------------------------------------------------
cat("Loading Pre-rendered JBrowse Schematic (Panel A)...\n")
fig5a_img <- ggdraw() + draw_image("data/images/fig5/figure_odol_variable.jpg", scale = 1)

# ------------------------------------------------------------------------------
# 5. Assemble Composite Figure 5
# ------------------------------------------------------------------------------
cat("Composing Figure 5...\n")

bottom_row <- plot_grid(p_boxplot, p_scatter, labels = c("B", "C"), ncol = 2, rel_widths = c(1, 1.2))

figure5 <- plot_grid(fig5a_img, bottom_row, labels = c("A", ""), ncol = 1, rel_heights = c(1.5, 1))

# ------------------------------------------------------------------------------
# 6. Save Outputs
# ------------------------------------------------------------------------------
output_pdf <- "figures/main/Figure5.pdf"
output_png <- "figures/main/Figure5.png"

cat("Saving Figure 5 to", output_pdf, "...\n")
ggsave(output_pdf, figure5, units = "in", width = 8, height = 8)

cat("Saving Figure 5 to", output_png, "...\n")
ggsave(output_png, figure5, units = "in", width = 8, height = 8, dpi = 300)

cat("Figure 5 successfully generated!\n")
