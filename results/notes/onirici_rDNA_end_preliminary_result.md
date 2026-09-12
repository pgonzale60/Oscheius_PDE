# Preliminary Result: Genomic Architecture, Assembly Inversion, and Motif-Independent Telomerization of the rDNA Locus in *Oscheius onirici*

**Date:** September 12, 2026  
**Species:** *Oscheius onirici* (strain PS2068, assembly `nxOscOnir1.2`)  
**Genomic Locus:** Chromosome II Left Terminus (rDNA Array)

---

## 1. Research Question & Biological Context
In *Oscheius* species undergoing Programmed DNA Elimination (PDE), somatic chromosomes undergo developmental scission and *de novo* telomere addition at distal ends guided by a conserved 27 bp cleavage motif. Chromosome II Left harbors a massive tandem array of 18S/28S ribosomal RNA genes ($\ge 1.2\text{ Mb}$, 344+ assembled units).

We sought to determine:
1. **Is diminution absent at the rDNA-bearing end, or does diminution occur within the rDNA array?**
2. **What is the true physical relationship between the rDNA array, the euchromatic core, and the terminal telomere?**
3. **Do the telomere-clipped reads and telotigs show evidence of the standard PDE cleavage motif?**

---

## 2. Empirical Findings

### A. Chromosome II Architecture & Assembly Misassembly
Examination of the primary sequence of Chromosome II (`nxOscOnir1.2.primary.fa.gz`) across the 0 to 1.25 Mb region uncovered a striking structural misassembly:

1. **Leftmost Edge (Coordinate 1 bp):**
   * Terminates abruptly in the middle of an rDNA repeat unit ($18S$ begins at pos 666 bp).
   * Possesses **0 telomeric hexamer repeats** (`TTAGGC`/`GCCTAA`) in the terminal 2,000 bp.
2. **The Internal Telomere Junction (Coordinate 1,185,509 bp):**
   * The tandem array runs from 1 bp to 1,185,509 bp.
   * At **pos 1,185,509 bp**, the sequence contains the **exact rDNA–telomere junction** captured in telotig `ptg000014l` ($100\%$ identity across 87 bp).
   * 3 raw telomere-clipped PacBio reads map directly to this coordinate.
3. **The Inversion & Stitch into Euchromatin (1,185,510 – 1,190,877 bp):**
   * Instead of continuing into a telomere, the assembler stitched a **4,703 bp inverted rDNA segment** (the reverse complement of 18S to 28S pos 4,703) directly downstream of the junction.
   * Coordinates 1,190,212 – 1,190,877 bp represent the reverse complement of the 666 bp non-coding spacer seen at pos 1 – 666 bp.
   * The true **euchromatic chromosome core** (containing canonical protein-coding genes like `gene-8173`) begins cleanly at **coordinate 1,190,878 bp**.
4. **Assembly Pathology Conclusion:**
   The assembler inverted the terminal rDNA cassette, stitching the true telomeric junction internally into the boundary of the euchromatic core while leaving coordinate 1 bp dangling without a telomere. The physical connection between the rDNA array and the chromosome core at 1.191 Mb is genuine.

---

### B. Somatic Read Depth Across the rDNA Array
* Whole-worm PacBio sequencing depth (`mosdepth` 1 kb windows):
  * **Chr II Core (> 1.5 Mb):** $94.4\times$ mean depth.
  * **Chr II rDNA (0 – 1.2 Mb):** $71.8\times$ mean depth ($\sim 76.1\%$ of core).
  * **Contrast with bona fide GRS (Chr II Right GRS):** Drops from $94\times$ in core to $11.7\times$ in GRS ($12\%$ residual depth).
* **Takeaway:** The soma unequivocally retains the massive $\ge 1.2\text{ Mb}$ rDNA array; there is no whole-locus elimination.

---

### C. Telomere-Clipped Read Clustering Reveals Two Inward-Facing Junctions
Projection of all 109 telomere-clipped raw PacBio reads mapping across the rDNA locus onto a canonical 6,958 bp rDNA repeat unit (18S: 1..1,754 bp; ITS/5.8S: 1,755..2,200 bp; 28S: 2,201..6,902 bp; IGS: 6,903..6,958 bp) revealed **two distinct, prominent clusters**:

1. **Cluster 1 ($N = 55$ reads, Orientation `R`):**
   * Maps sharply to **position ~4,705 bp** (inside the **28S rRNA gene**).
   * Represents a forward-strand junction where rDNA transitions rightward into telomere repeats: `[...rDNA...28S(4705)] ── (TTAGGC)n`.
   * Matches the exact assembly junction found in telotig `ptg000014l` (pos 4,703 bp).
   * Possesses massive telomere tracts (mean $>1,200$ repeats).
2. **Cluster 2 ($N = 44$ reads, Orientation `L`):**
   * Maps sharply to **position ~1,637 bp** (distal end of the **18S rRNA gene**).
   * Represents a reverse-strand junction where telomere repeats transition into rDNA: `(GCCTAA)n ── [18S(1637)...28S...]`.
   * Possesses massive telomere tracts (mean $>1,200$ repeats).
3. **Non-cluster Singletons ($N = 10$ reads):**
   * Isolated singletons ($n = 1$ read per site) scattered across the repeat unit.
   * Low repeat counts (median 62 repeats, some $\le 34$) and MAPQ $\le 1$, consistent with stochastic PCR-free library chimeras / sequencing noise.
4. **Summary:** 99 of 109 reads ($90.8\%$) support exactly two specific, inward-facing telomeric cleavage junctions.

---

### D. Complete Absence of the 27 bp PDE Cleavage Motif
* The conserved *O. onirici* PDE consensus cleavage motif is:
  $$\text{5' - } \mathbf{TTWGGYTDWGGCYVAGCCNAABCCAAA} \text{ - 3'}$$
* In bona fide chromosomal PDE ends (such as Chromosome II Right at 10,439,112 bp):
  * The scission motif matches with extraordinary statistical significance (**$P = 3.16 \times 10^{-12}$**).
  * The enzymatic cut cleanly leaves the 5' half-motif directly preceding the *de novo* telomere addition (`...CTAATATTTTTTTGGATGAGGCCGA | (TTAGGC)n`).
* In stark contrast, sequence examination of both rDNA telomeric junctions:
  * **Cluster 1 (28S, pos ~4,705 bp):** `...TTGCTTTCGGGCTTAGTATCT | (TTAGGC)n`
  * **Cluster 2 (18S, pos ~1,637 bp):** `(GCCTAA)n | TTTCGAGAAGACTGAG...`
  * Neither junction displays sequence homology to the 27 bp motif, nor contains a 5' half-motif (MEME $P = 0.826$).
* **Conclusion:** Telomerization at the rDNA locus operates via a **motif-independent mechanism**, distinct from standard PDE scission.

---

## 3. Biological Synthesis: Extrachromosomal rDNA Minichromosomes vs. Array Capping

The discovery of two distinct, opposing telomeric junctions inside coding ribosomal genes (28S at ~4,705 bp and 18S at ~1,637 bp) points to two compelling biological models:

1. **Extrachromosomal Linear rDNA Minichromosomes:**
   * In organisms undergoing developmental chromatin diminution (such as the nematode *Ascaris suum* or ciliates like *Tetrahymena*), the chromosomal rDNA array is excised from the germline chromosomes during development and processed into independent, linear extrachromosomal elements capped by telomeres on both ends:
     $$\text{[Telomere]} \longleftarrow \mathbf{18S}_{(1637..1754)} \text{—[ITS/28S/IGS]—[Full Units]}_n\text{—}\mathbf{28S}_{(2201..4705)} \longrightarrow \text{[Telomere]}$$
   * Cluster 1 ($N=55$) and Cluster 2 ($N=44$) represent the two terminal caps of these linear minichromosomes.
2. **Specialized Array Healing / Truncation:**
   * Alternatively, if the array remains contiguous with the chromosome core at 1.19 Mb, somatic cells undergo a specialized, motif-independent double-strand break that truncates the distal units of the array and caps it with telomeres.

Under either model, **the rDNA terminus is not governed by the canonical 27 bp PDE machinery**.

---

## 4. Reproducibility
* Standalone validation and sequence check script:
  ```bash
  Rscript analysis/investigate_onirici_rDNA_terminus.R
  ```
* Multi-panel figure generation script:
  ```bash
  Rscript analysis/plot_figureS6_onirici_rDNA_terminus.R
  ```
