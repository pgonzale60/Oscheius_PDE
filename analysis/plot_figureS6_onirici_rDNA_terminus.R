#!/usr/bin/env Rscript
# ==============================================================================
# Figure S6: Resolution of the Chromosome II rDNA Terminus in Oscheius onirici
# 
# Panel A: Chromosome II Architecture, Assembly Inversion, & Depth
# Panel B: Canonical rDNA Unit, Junction to ptg000014l, & Read Support
# Panel C: Direct Sequence Contrast: Chr II Left (No PDE) vs Chr II Right (PDE)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(Biostrings)
  library(pwalign)
})

# ------------------------------------------------------------------------------
# 1. Source Themes & Paths
# ------------------------------------------------------------------------------
source("R/plot_theme.R")

output_pdf <- "figures/supplementary/FigureS6.pdf"
output_png <- "figures/supplementary/FigureS6.png"

# ------------------------------------------------------------------------------
# 2. Panel A: Chromosome II Architecture & Assembly Conflict
# ------------------------------------------------------------------------------
df_chr <- tibble(
  xmin = c(0, 1.185509, 10.439112),
  xmax = c(1.185509, 10.439112, 10.440939),
  ymin = 0.7,
  ymax = 1.3,
  region = factor(c("rDNA Array (~1.2 Mb)", "Euchromatic Core", "GRS (Eliminated)"),
                  levels = c("rDNA Array (~1.2 Mb)", "Euchromatic Core", "GRS (Eliminated)"))
)

p_a <- ggplot() +
  # Chromosome body
  geom_rect(data = df_chr, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = region),
            color = "black", linewidth = 0.4) +
  
  # Coordinate 1 marker (Left assembly truncation)
  annotate("segment", x = 0, xend = 0, y = 1.45, yend = 1.7, arrow = arrow(length = unit(0.2, "cm"), ends = "first"), color = "#d95f02", linewidth = 0.8) +
  annotate("text", x = 0, y = 1.88, label = "Coord 1 bp: Truncated in assembly\n(0 telomeric repeats in 2 kb)", 
           hjust = 0, size = 3.3, fontface = "bold", color = "#d95f02", lineheight = 0.95) +
  
  # Coordinate 1.185 Mb marker (rDNA-core junction where reads actually map)
  annotate("segment", x = 1.185509, xend = 1.185509, y = 0.55, yend = 0.3, arrow = arrow(length = unit(0.2, "cm"), ends = "first"), color = "#7570b3", linewidth = 0.8) +
  annotate("text", x = 1.185509, y = 0.15, label = "Coord 1.185 Mb: True rDNA-telomere junction\n(Inverted cassette stitched into euchromatin)", 
           hjust = 0, size = 3.3, fontface = "bold", color = "#7570b3", lineheight = 0.95) +
  
  # Right end PDE marker
  annotate("segment", x = 10.439112, xend = 10.439112, y = 1.45, yend = 1.7, arrow = arrow(length = unit(0.2, "cm"), ends = "first"), color = "#1b9e77", linewidth = 0.8) +
  annotate("text", x = 10.439112, y = 1.88, label = "Coord 10.439 Mb: Bona Fide PDE End\n(141 clipped reads, P = 3.2e-12)", 
           hjust = 1, size = 3.3, fontface = "bold", color = "#1b9e77", lineheight = 0.95) +
  
  scale_fill_manual(values = c("rDNA Array (~1.2 Mb)" = "#fc8d62", 
                               "Euchromatic Core" = "#8da0cb", 
                               "GRS (Eliminated)" = "#e78ac3")) +
  scale_x_continuous(breaks = seq(0, 10, by = 2), limits = c(-0.5, 11), labels = function(x) paste0(x, " Mb")) +
  ylim(-0.1, 2.2) +
  labs(title = "A   Chromosome II Architecture & Assembly Inversion in Oscheius onirici",
       x = "Chromosome II Assembly Coordinate (Mb)",
       fill = "Genomic Region") +
  theme_pde() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.box.margin = margin(t = -5, b = 5)
  )

# ------------------------------------------------------------------------------
# 3. Panel B: Canonical rDNA Unit & Read Support
# ------------------------------------------------------------------------------
df_unit <- tibble(
  xmin = c(1, 1755, 2201, 6903),
  xmax = c(1754, 2200, 6902, 6958),
  ymin = 0.7,
  ymax = 1.3,
  feature = factor(c("18S rRNA", "ITS1 / 5.8S / ITS2", "28S rRNA", "IGS / Spacer"),
                   levels = c("18S rRNA", "ITS1 / 5.8S / ITS2", "28S rRNA", "IGS / Spacer"))
)

clipped_data <- tibble(
  rel_pos = c(637, 1624, 1637, 2314, 2559, 2724, 3956, 4705, 5601, 5784, 6375, 6586),
  n_reads = c(1, 1, 44, 1, 1, 1, 1, 55, 1, 1, 1, 1),
  group = factor(c(rep("Non-cluster Singletons (n = 1)", 2), 
                   "Cluster 2: Rev-strand entry (L)", 
                   rep("Non-cluster Singletons (n = 1)", 4), 
                   "Cluster 1: Fwd-strand / ptg000014l (R)", 
                   rep("Non-cluster Singletons (n = 1)", 4)),
                 levels = c("Cluster 1: Fwd-strand / ptg000014l (R)", 
                            "Cluster 2: Rev-strand entry (L)", 
                            "Non-cluster Singletons (n = 1)"))
)

p_b <- ggplot() +
  # Unit gene model
  geom_rect(data = df_unit, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = feature),
            color = "black", linewidth = 0.3) +
  geom_text(aes(x = (1 + 1754)/2, y = 1.0, label = "18S rRNA"), fontface = "bold", size = 3.2) +
  geom_text(aes(x = (2201 + 6902)/2, y = 1.0, label = "28S rRNA"), fontface = "bold", size = 3.2) +
  
  # Stems and points
  geom_segment(data = clipped_data, aes(x = rel_pos, xend = rel_pos, y = 1.3, yend = 1.3 + n_reads * 0.055),
               color = "grey50", linewidth = 0.5) +
  geom_point(data = clipped_data, aes(x = rel_pos, y = 1.3 + n_reads * 0.055, color = group, size = group)) +
  
  # Labels
  annotate("text", x = 4705, y = 4.7, 
           label = "Cluster 1 (N = 55 reads)\nptg000014l junction (pos 4,703)\nFwd strand (R) -> (TTAGGC)n\nMean >1,200 reps",
           size = 3.1, fontface = "bold", color = "#1b9e77", hjust = 0, lineheight = 0.95) +
  
  annotate("text", x = 1637, y = 4.1, 
           label = "Cluster 2 (N = 44 reads)\nRev strand (L) <- (GCCTAA)n\nDistal 18S gene boundary\nMean >1,200 reps",
           size = 3.1, fontface = "bold", color = "#7570b3", hjust = 1, lineheight = 0.95) +
  
  annotate("text", x = 3300, y = 1.8, 
           label = "10 Singletons (n = 1 each, median 62 reps, library chimeras)",
           size = 2.9, fontface = "italic", color = "grey35") +
  
  scale_fill_manual(values = c("18S rRNA" = "#66c2a5", 
                               "ITS1 / 5.8S / ITS2" = "#ffd92f", 
                               "28S rRNA" = "#8da0cb", 
                               "IGS / Spacer" = "#e78ac3")) +
  scale_color_manual(values = c("Cluster 1: Fwd-strand / ptg000014l (R)" = "#1b9e77",
                                "Cluster 2: Rev-strand entry (L)" = "#7570b3",
                                "Non-cluster Singletons (n = 1)" = "#e41a1c")) +
  scale_size_manual(values = c("Cluster 1: Fwd-strand / ptg000014l (R)" = 4.5,
                               "Cluster 2: Rev-strand entry (L)" = 4.0,
                               "Non-cluster Singletons (n = 1)" = 2.0)) +
  scale_x_continuous(limits = c(-100, 7500), breaks = seq(0, 7000, 1000), labels = function(x) paste0(x, " bp")) +
  ylim(0.4, 5.2) +
  labs(title = "B   Canonical 6,958 bp rDNA Unit Architecture & Telomere-Clipped Read Support",
       x = "Position along Canonical rDNA Unit (bp)",
       y = "Read Support",
       fill = "rDNA Component",
       color = "Telomeric Read Evidence",
       size = "Telomeric Read Evidence") +
  theme_pde() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = -5, b = 5)
  )

# ------------------------------------------------------------------------------
# 4. Panel C: Direct Sequence Contrast: Chr II Left vs Chr II Right
# ------------------------------------------------------------------------------
p_c <- ggplot() +
  # Background cards
  annotate("rect", xmin = 0.4, xmax = 4.85, ymin = 0.4, ymax = 4.85, fill = "#f8f9fa", color = "#cbd5e1", linewidth = 0.6) +
  annotate("rect", xmin = 5.15, xmax = 9.6, ymin = 0.4, ymax = 4.85, fill = "#f8f9fa", color = "#cbd5e1", linewidth = 0.6) +
  
  # Headers
  annotate("text", x = 2.62, y = 4.5, label = "Chromosome II Left (rDNA Terminus)", fontface = "bold", size = 4.0, color = "#2b83ba") +
  annotate("text", x = 7.38, y = 4.5, label = "Chromosome II Right (Bona Fide PDE End)", fontface = "bold", size = 4.0, color = "#d7191c") +
  
  # Left details
  annotate("text", x = 0.7, y = 3.9, label = "Biological Status:", fontface = "bold", hjust = 0, size = 3.3) +
  annotate("text", x = 2.3, y = 3.9, label = "NO DIMINUTION (Retained in Soma)", hjust = 0, size = 3.3, color = "#2b83ba", fontface = "bold") +
  
  annotate("text", x = 0.7, y = 3.4, label = "Sequencing Depth:", fontface = "bold", hjust = 0, size = 3.3) +
  annotate("text", x = 2.3, y = 3.4, label = "Uniform ~72x depth across 1.2 Mb array", hjust = 0, size = 3.3) +
  
  annotate("text", x = 0.7, y = 2.9, label = "Clipped Read Support:", fontface = "bold", hjust = 0, size = 3.3) +
  annotate("text", x = 2.3, y = 2.9, label = "99 reads confirm ptg000014l junction", hjust = 0, size = 3.3) +
  
  annotate("text", x = 0.7, y = 2.4, label = "MEME Motif Match:", fontface = "bold", hjust = 0, size = 3.3) +
  annotate("text", x = 2.3, y = 2.4, label = "P = 0.826 (MOTIF ABSENT)", hjust = 0, size = 3.3, fontface = "bold", color = "#d7191c") +
  
  annotate("text", x = 0.7, y = 1.8, label = "Sequence at Terminus:", fontface = "bold", hjust = 0, size = 3.3) +
  annotate("text", x = 0.7, y = 1.25, 
           label = "rDNA:  ... T T G C T T T C G G G C T T A G T A T C T ...\nTelo:   | ( T T A G G C ) n   [16.3 kb pure tract]", 
           family = "Courier", hjust = 0, size = 2.9, color = "black", lineheight = 1.1) +

  # Right details
  annotate("text", x = 5.45, y = 3.9, label = "Biological Status:", fontface = "bold", hjust = 0, size = 3.3) +
  annotate("text", x = 7.15, y = 3.9, label = "PDE SCISSION (1.8 kb GRS Eliminated)", hjust = 0, size = 3.3, color = "#d7191c", fontface = "bold") +
  
  annotate("text", x = 5.45, y = 3.4, label = "Sequencing Depth:", fontface = "bold", hjust = 0, size = 3.3) +
  annotate("text", x = 7.15, y = 3.4, label = "Drops 94x (Core) -> 11x (GRS residual)", hjust = 0, size = 3.3) +
  
  annotate("text", x = 5.45, y = 2.9, label = "Clipped Read Support:", fontface = "bold", hjust = 0, size = 3.3) +
  annotate("text", x = 7.15, y = 2.9, label = "141 reads at exact coord 10,439,112 bp", hjust = 0, size = 3.3) +
  
  annotate("text", x = 5.45, y = 2.4, label = "MEME Motif Match:", fontface = "bold", hjust = 0, size = 3.3) +
  annotate("text", x = 7.15, y = 2.4, label = "P = 3.2 x 10^-12 (HIGHLY SIGNIFICANT)", hjust = 0, size = 3.3, fontface = "bold", color = "#1b9e77") +
  
  annotate("text", x = 5.45, y = 1.8, label = "Cleavage Scar & Motif:", fontface = "bold", hjust = 0, size = 3.3) +
  annotate("text", x = 5.45, y = 1.25, 
           label = "Motif:   T T G G A T G A G G C C G A G C C G A A A C C T A A T\nScar:  ... C T A A T A T T T T T T T G G A T G A G G C C G A | (TTAGGC)n", 
           family = "Courier", hjust = 0, size = 2.9, color = "black", lineheight = 1.1) +

  xlim(0, 10) +
  ylim(0.3, 5.0) +
  labs(title = "C   Direct Within-Chromosome Contrast: Chromosome II Left vs. Chromosome II Right") +
  theme_void() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0, family = "Helvetica", margin = margin(t = 10, b = 10))
  )

# ------------------------------------------------------------------------------
# 5. Assemble Multi-Panel Figure
# ------------------------------------------------------------------------------
final_fig <- plot_grid(
  p_a,
  p_b,
  p_c,
  ncol = 1,
  rel_heights = c(1.0, 1.35, 1.2)
)

cat("Saving figure to:\n")
cat(" - PDF:", output_pdf, "\n")
cat(" - PNG:", output_png, "\n")

ggsave(output_pdf, final_fig, width = 11, height = 11, device = "pdf")
ggsave(output_png, final_fig, width = 11, height = 11, dpi = 300)

cat("Figure generation complete!\n")
