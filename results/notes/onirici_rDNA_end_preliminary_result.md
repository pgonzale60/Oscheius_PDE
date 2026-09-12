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

### D. Complete Absence of the 27 bp PDE Cleavage Motif
* The conserved *O. onirici* PDE consensus cleavage motif is:
  $$\text{5' - } \mathbf{TTWGGYTDWGGCYVAGCCNAABCCAAA} \text{ - 3'}$$
* In true somatic diminution sites, enzymatic scission leaves the 5' half-motif immediately preceding the *de novo* telomere addition.
* In `ptg000014l`, the sequence immediately upstream of pos 23,587 (`...TTGCTTTCGGGCTTAGTATCT - [TTAGGC]n`) **completely lacks the 27 bp motif and contains no 5' half-motif**.
* This matches the MEME motif search (`nxOscOnir1.2.meme.txt`), where 10 of the 12 chromosome ends had significant motif matches ($P \sim 10^{-10}\text{ to }10^{-14}$), while Chr II Left scored **$P = 0.826$ (no motif detected)**.

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
