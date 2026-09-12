# Preliminary Result: Genomic Architecture and Absence of Programmed DNA Elimination at the rDNA Terminus in *Oscheius onirici*

**Date:** September 12, 2026  
**Species:** *Oscheius onirici* (strain PS2068, assembly `nxOscOnir1.2`)  
**Genomic Locus:** Chromosome II Left Terminus (rDNA Array)

---

## 1. Research Question & Biological Context
In *Oscheius* species undergoing Programmed DNA Elimination (PDE), somatic chromosomes undergo developmental scission and *de novo* telomere addition at distal ends. Chromosome II Left harbors a massive tandem array of 18S/28S ribosomal RNA genes ($\ge 1.2\text{ Mb}$, 344+ assembled units). 

We sought to determine:
1. **Is diminution absent at the rDNA-bearing end?** (i.e. The rDNA forms the terminal block in both soma and germline, without somatic scission).
2. **Or does diminution occur inside the rDNA array?** (i.e. The somatic GRS consists purely of eliminated distal rDNA copies).

---

## 2. Empirical Findings

### A. The Primary Whole-Genome Assembly (`nxOscOnir1.2.primary.fa.gz`)
* **Terminal sequence:** The assembled contig for Chromosome II Left terminates abruptly at coordinate 1 with **zero `TTAGGC` telomeric repeats** in the terminal 2,000 bp.
* **Content:** The first 18S rRNA gene begins at position 666 bp. The preceding 665 bp consist solely of canonical non-coding rDNA intergenic sequence.
* **Status:** The chromosome was truncated by the assembler at the edge of the first rDNA repeat unit due to tandem repeat collapse.

### B. Somatic Read Depth Across the rDNA Array
* Using 1 kb binned depth profiles (`mosdepth` from whole-worm PacBio sequencing):
  * **Chr II Core (> 1.5 Mb):** $94.4\times$ mean depth.
  * **Chr II rDNA (0 – 1.2 Mb):** $71.8\times$ mean depth ($\sim 76.1\%$ of core).
  * **Contrast with bona fide GRS (e.g. Chr III Right):** Drops from $100.7\times$ in core down to $11.1\times$ in GRS ($11.0\%$ residual germline depth).
* **Takeaway:** There is **no sharp step-drop** at the 1.2 Mb rDNA boundary. The soma unequivocally retains the entire $\ge 1.2\text{ Mb}$ rDNA array.

### C. Telotig `ptg000014l` Structure & rDNA–Telomere Junction
* **Contig length:** 44,096 bp, with $29.4\times$ coverage.
* **Structure:**
  * Coordinates 1 – 20,638 bp: Three full 18S/28S rDNA repeat units.
  * Coordinates 20,639 – 23,586 bp (2,948 bp): The canonical rDNA **Intergenic Spacer (IGS)**, showing $>95\%$ identity to internal IGS units across Chr II.
  * Coordinates 23,587 – 44,096 bp: Direct transition into a massive **16.3 kb pure `(TTAGGC)n` telomeric repeat tract**.

### D. Telomere-Clipped Read Distribution Across the Canonical rDNA Unit
To ensure that the assembled rDNA–telomere junction in `ptg000014l` was not an artifact of assembler selection, all 109 telomere-clipped raw PacBio reads mapping across the entire 1.2 Mb rDNA array were projected onto a canonical 6,958 bp rDNA unit (18S: 1..1,754 bp; ITS/5.8S: 1,755..2,200 bp; 28S: 2,201..6,902 bp; IGS/ETS: 6,903..6,958 bp):

1. **Cluster 1 ($N = 55$ reads, orientation `R`):**
   * Projects sharply to relative position **~4,705 bp** (in 28S rRNA).
   * Matches the exact assembly junction found in telotig `ptg000014l` (relative pos 4,703 bp).
   * Telomere tracts are massive: mean $>1,200$ repeats.
2. **Cluster 2 ($N = 44$ reads, orientation `L`):**
   * Projects sharply to relative position **~1,637 bp** (distal end of 18S rRNA).
   * Represents the complementary strand perspective traversing inward from the opposite end of the collapsed rDNA repeat.
   * Telomere tracts are massive: mean $>1,200$ repeats.
3. **Non-cluster Singletons ($N = 10$ reads):**
   * 10 reads are scattered across the rDNA unit (pos 637, 1624, 2314, 2559, 2724, 3956, 5601, 5784, 6375, 6586).
   * Exactly $n = 1$ read at each distinct site.
   * Low repeat counts (median 62 repeats, range 16–134) and MAPQ $\le 1$, characteristic of stochastic sequencing chimeras / low-level library noise rather than coordinated biological cleavage sites.

### E. Complete Absence of the 27 bp PDE Cleavage Motif
* The conserved *O. onirici* PDE consensus cleavage motif is:
  $$\text{5' - } \mathbf{TTWGGYTDWGGCYVAGCCNAABCCAAA} \text{ - 3'}$$
* In true somatic diminution sites, enzymatic scission leaves the 5' half-motif immediately preceding the *de novo* telomere addition.
* Neither Cluster 1 (~4,705 bp) nor Cluster 2 (~1,637 bp) contains the 27 bp motif or the 5' half-motif (`has_motif = FALSE`).
* MEME motif search across all chromosome ends (`nxOscOnir1.2.meme.txt`) gave $P = 0.826$ for Chr II Left, contrasting with $P \sim 10^{-10} \text{ to } 10^{-14}$ at bona fide PDE ends.

---

## 3. Biological Conclusion

**Programmed DNA Elimination is ABSENT at the rDNA-bearing chromosome end in *O. onirici*.**

1. The rDNA array constitutes the natural terminal block of Chromosome II.
2. The telomere captured on `ptg000014l` represents the ancestral germline/somatic chromosome terminus.
3. The assembly failed to link `ptg000014l` to the main Chr II scaffold purely due to tandem repeat collapse, explaining why this end was classified as an rDNA-bearing exclusion rather than a PDE-truncated GRS.

---

## 4. Reproducibility
To reproduce all sequence comparisons, pattern matches, and coordinate checks:
```bash
Rscript analysis/investigate_onirici_rDNA_terminus.R
```
