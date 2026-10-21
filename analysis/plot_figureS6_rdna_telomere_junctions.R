#!/usr/bin/env Rscript
# ==============================================================================
# Figure S6: Telomere-bearing reads at the rDNA array of three Oscheius species
#
# Question: does the terminal rDNA array undergo somatic breakage? If it did, telomere
# junctions would be recovered at many positions (or at supra-stoichiometric depth). Every
# long read that carries a telomeric tract and maps to the array contributes one junction. The
# non-telomeric flank of the tract is mapped to ONE rDNA monomer per species, so position and
# orientation do not depend on how the assembly arranged (or collapsed / inverted) the copies.
#
#   above the unit model : telomere lies downstream of the junction on the monomer + strand
#   below the unit model : telomere lies upstream (monomer - strand runs into it)
#   (reads from either DNA strand are pooled: a junction is one molecular site)
#   height               : reads / somatic autosomal coverage  (1.0 = one junction per genome)
#   shape                : species;  filled = >= 2 reads at that position, open = singleton
#
# Inputs (all small tables in the repo; no BAMs):
#   data/rdna_telomere/<asm>.monomer_junctions.tsv
#       (scripts/rdna_telomere_junctions/: run_extraction.sh -> build_monomers.py ->
#        map_junctions_to_monomer.py)
#   data/rnammer/<asm>.rnammer.gff3       (unit model and phase anchors)
#   data/coverage/<asm>.mosdepth.summary.txt (somatic baseline)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})
source("R/plot_theme.R")

output_pdf <- "figures/supplementary/FigureS6.pdf"
output_png <- "figures/supplementary/FigureS6.png"
dir.create("figures/supplementary", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Species table
# ------------------------------------------------------------------------------
sp <- tribble(
  ~asm,           ~strain,  ~species,           ~rdna_chr, ~lo,       ~hi,        ~shape,
  "nxOscOnir1.2", "PS2068", "O. onirici",       "II",      1,         1260000,    21,
  "nxOscSper1.1", "DF5120", "O. sp. DF5120",    "II",      8130000,   8180000,    24,
  "nxOscDolc1.1", "PS1017", "O. dolichura",     "IV",      12450000,  12464839,   22
)

# ------------------------------------------------------------------------------
# 2. Somatic baseline: median of mean coverage over the assembled autosomes I-V, excluding the
#    rDNA chromosome (unplaced scaffolds, "_region" rows and X are not used)
# ------------------------------------------------------------------------------
baseline <- map_dfr(seq_len(nrow(sp)), function(i) {
  read_tsv(sprintf("data/coverage/%s.mosdepth.summary.txt", sp$asm[i]), show_col_types = FALSE) %>%
    filter(chrom %in% c("I", "II", "III", "IV", "V"), chrom != sp$rdna_chr[i]) %>%
    summarise(asm = sp$asm[i], baseline = median(mean))
})

# ------------------------------------------------------------------------------
# 3. rDNA unit anchors (18S start) and per-species unit model
# ------------------------------------------------------------------------------
read_rrna <- function(i) {
  read_tsv(sprintf("data/rnammer/%s.rnammer.gff3", sp$asm[i]), comment = "#",
           col_names = c("chr", "src", "type", "start", "end", "score", "strand", "phase", "attr"),
           show_col_types = FALSE) %>%
    mutate(gene = sub("^ID=([0-9]+)s_.*", "\\1", attr)) %>%
    filter(chr == sp$rdna_chr[i], start >= sp$lo[i], start <= sp$hi[i], gene %in% c("18", "28")) %>%
    arrange(start)
}
rrna <- map(seq_len(nrow(sp)), read_rrna) %>% set_names(sp$asm)

unit_len <- map_dbl(rrna, function(r) {
  s18 <- r$start[r$gene == "18"]
  if (length(s18) >= 2) median(diff(s18)) else {
    s28 <- r$start[r$gene == "28"]            # PS1017 collapsed array: 28S start to 28S start
    median(diff(s28))
  }
})
anchors <- map(rrna, ~ .x$start[.x$gene == "18"])

models <- map_dfr(seq_len(nrow(sp)), function(i) {
  r <- rrna[[i]]; a <- anchors[[i]][1]; ul <- unit_len[[i]]
  g18 <- r %>% filter(gene == "18") %>% slice(1)
  g28 <- r %>% filter(gene == "28", start > g18$start) %>% slice(1)
  # 28S length: median over copies when the array has several; longest copy when it has few
  # (the collapsed PS1017 array ends in a 28S truncated by the scaffold end)
  len28 <- (r %>% filter(gene == "28") %>% mutate(len = end - start))$len
  len28 <- if (length(len28) >= 5) median(len28) else max(len28)
  tibble(asm = sp$asm[i],
         `18S` = c((g18$start - a) / ul, (g18$end - a) / ul),
         `28S` = c((g28$start - a) / ul, (g28$start - a + len28) / ul)) %>%
    mutate(edge = c("xmin", "xmax")) %>%
    pivot_longer(c(`18S`, `28S`), names_to = "gene", values_to = "x") %>%
    pivot_wider(names_from = edge, values_from = x)
})
cat("Unit lengths (bp):\n"); print(unit_len)
cat("Somatic baselines:\n"); print(baseline)

# ------------------------------------------------------------------------------
# 4. Telomere junctions placed on a single monomer per species
#    (flank sequence mapped to the monomer: scripts/rdna_telomere_junctions/)
# ------------------------------------------------------------------------------
jn <- map_dfr(seq_len(nrow(sp)), function(i) {
  read_tsv(sprintf("data/rdna_telomere/%s.monomer_junctions.tsv", sp$asm[i]),
           show_col_types = FALSE, col_types = cols(.default = "c")) %>%
    mutate(asm = sp$asm[i])
}) %>%
  select(-strain) %>%
  left_join(sp, by = "asm") %>% left_join(baseline, by = "asm")

cat("\nFlanks placed on the monomer:\n")
print(jn %>% group_by(strain) %>% summarise(flanks = n(), mapped = sum(mapped == "1")))

jn <- jn %>%
  filter(mapped == "1") %>%
  mutate(mono_phase = as.numeric(mono_phase),
         ysign = if_else(telomere_side == "downstream", 1, -1),
         bin = round(mono_phase / 0.01) * 0.01)

pts <- jn %>%
  group_by(strain, species, shape, baseline, telomere_side, ysign, bin) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(frac = n / baseline, singleton = n == 1)

tot <- jn %>% count(strain, baseline) %>% mutate(ratio = n / baseline)
cat("\nTotals (junctions, genome-equivalents):\n"); print(tot)
cat("\nDominant position per species:\n")
print(pts %>% group_by(strain) %>% slice_max(n, n = 2, with_ties = FALSE) %>% select(strain, telomere_side, bin, n))

lab <- tot %>% mutate(label = sprintf("%s: %d reads (%.2f x)", strain, n, ratio))
pts <- pts %>% left_join(select(lab, strain, label), by = "strain")

# ------------------------------------------------------------------------------
# 5. Plot
# ------------------------------------------------------------------------------
band <- 0.22                      # half-height of the central unit-model band
track_y <- c(nxOscOnir1.2 = 0.12, nxOscSper1.1 = 0.00, nxOscDolc1.1 = -0.12)
mod <- models %>%
  mutate(y = track_y[asm], gene_fill = gene) %>%
  left_join(select(sp, asm, strain, shape), by = "asm")

pts <- pts %>% mutate(y = ysign * (band + frac))

ymax <- ceiling((max(pts$frac) + 0.02) / 0.1) * 0.1
brk <- seq(0, ymax, by = max(0.1, round(ymax / 4, 1)))

shape_vals <- setNames(sp$shape, sp$strain)
lab_levels <- lab$label[match(sp$strain, lab$strain)]
shape_by_label <- setNames(sp$shape, lab_levels)
pts$label <- factor(pts$label, levels = lab_levels)

p <- ggplot() +
  geom_hline(yintercept = c(-band, band), colour = "grey80", linewidth = 0.3) +
  geom_rect(data = mod, aes(xmin = xmin, xmax = xmax, ymin = y - 0.045, ymax = y + 0.045, fill = gene),
            colour = "black", linewidth = 0.25) +
  geom_text(data = mod %>% filter(gene == "18S"),
            aes(x = -0.012, y = y, label = strain), hjust = 1, size = 2.8) +
  geom_segment(data = pts, aes(x = bin, xend = bin, y = ysign * band, yend = y),
               colour = "grey55", linewidth = 0.35) +
  geom_point(data = pts, aes(x = bin, y = y, shape = label, fill = if_else(singleton, "open", "filled")),
             size = 2.6, colour = "black", stroke = 0.5) +
  annotate("text", x = 1.0, y = band + ymax + 0.14, label = "telomere downstream of junction (monomer + strand)", hjust = 1, vjust = 1, size = 3, fontface = "bold") +
  annotate("text", x = 1.0, y = -(band + ymax) - 0.14, label = "telomere upstream of junction (monomer - strand)", hjust = 1, vjust = 0, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("18S" = "#66c2a5", "28S" = "#8da0cb", "filled" = "black", "open" = "white"),
                    breaks = c("filled", "open"),
                    labels = c(filled = ">= 2 reads at position", open = "singleton"),
                    name = NULL) +
  scale_shape_manual(values = shape_by_label, name = NULL) +
  scale_x_continuous(limits = c(-0.09, 1.0), breaks = seq(0, 1, 0.25),
                     labels = function(x) sprintf("%.2f", x),
                     name = "Telomere junction position on the rDNA monomer (fraction of each species' monomer; 18S start = 0)",
                     sec.axis = sec_axis(~ . * median(unit_len), breaks = seq(0, 7000, 1000),
                                         labels = function(x) paste0(x / 1000, " kb"),
                                         name = sprintf("bp on the median unit (%s bp)", format(round(median(unit_len)), big.mark = ",")))) +
  scale_y_continuous(breaks = c(-(band + brk), band + brk), labels = c(brk, brk),
                     limits = c(-(band + ymax) - 0.2, band + ymax + 0.2),
                     name = "Telomere-bearing reads / somatic coverage") +
  guides(shape = guide_legend(order = 1, nrow = 1, override.aes = list(fill = "grey70")),
         fill = guide_legend(order = 2, override.aes = list(shape = 21, size = 2.6))) +
  labs(title = "Telomere junctions along the rDNA unit in three Oscheius species",
       subtitle = "18S (green) and 28S (blue) drawn on each species' own repeat unit; height 1.0 = one junction per somatic genome") +
  theme_pde() +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.box = "vertical",
        plot.title = element_text(face = "bold", size = 12))

ggsave(output_pdf, p, width = 9, height = 6.5, device = "pdf")
ggsave(output_png, p, width = 9, height = 6.5, dpi = 300)

# plotted table for the supplement / reuse
write_tsv(pts %>% select(strain, telomere_side, bin, n, frac, singleton),
          "results/tables/FigureS6_rdna_telomere_positions.tsv")
cat("Figure S6 written.\n")
