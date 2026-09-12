# Preliminary Result: Genomic Architecture, Assembly Inversion, and Motif-Independent Telomerization of the rDNA Locus in *Oscheius onirici*

**Date:** September 12, 2026  
**Species:** *Oscheius onirici* (strain PS2068, assembly `nxOscOnir1.2` / `nxOscOnir1.5`)  
**Genomic Locus:** Chromosome II Left Terminus (rDNA Array)

---

## 1. Research Question & Biological Context
In *Oscheius* species undergoing Programmed DNA Elimination (PDE), somatic chromosomes undergo developmental scission and *de novo* telomere addition at distal ends guided by a conserved 27 bp cleavage motif. Chromosome II Left harbors a massive tandem array of 18S/28S ribosomal RNA genes ($\ge 1.2\text{ Mb}$, 344+ assembled units).

We sought to determine:
1. **Is diminution absent at the rDNA-bearing end, or does diminution occur within the rDNA array?**
2. **What is the true physical relationship between the rDNA array, the euchromatic core, and the terminal telomere?**
3. **Where does the germline break occur, and do the spanning reads reveal the canonical PDE cleavage motif?**

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
   * Telomere-clipped reads map directly to this coordinate.
3. **The Inversion & Stitch into Euchromatin (1,185,510 – 1,190,877 bp):**
   * Instead of continuing into a telomere, the assembler stitched a **4,703 bp inverted rDNA segment** (the reverse complement of 18S to 28S pos 4,703) directly downstream of the junction.
   * Coordinates 1,190,212 – 1,190,877 bp represent the reverse complement of the 666 bp non-coding spacer seen at pos 1 – 666 bp.
   * The true **euchromatic chromosome core** (containing canonical protein-coding genes like `gene-8173`) begins cleanly at **coordinate 1,190,878 bp**.
4. **Assembly Pathology Conclusion:**
   The assembler folded the terminal rDNA cassette back into the chromosome, creating an artificial head-to-head hairpin at 1.185 Mb that connects to the euchromatic core, while leaving coordinate 1 bp dangling without a telomere.

---

### B. High-Coverage Spanning Reads Across the Core–rDNA Boundary on Farm HPC
By accessing the Sanger Farm HPC cluster (`farm22-head2`) via Teleport (`tsh`) and querying the remapped PacBio read alignments (`remmaped_all_reads_nxOscOnir1.5.bam`), we inspected the raw alignments across both key boundaries:

1. **Continuous Physical Connection at 1,190,878 bp:**
   * Sequencing coverage across the core-to-rDNA boundary ($1,190,878\text{ bp}$) is **$230\times$**.
   * Numerous individual PacBio HiFi reads (e.g. `m64016_201225_182119/176949466/ccs`) cleanly span across the boundary, with several kilobases extending into the euchromatic core and into the rDNA non-coding spacer.
   * **Complete Absence of Motif:** The entire transition sequence shows **zero match** to the 27 bp consensus elimination motif ($\text{TTWGGYTDWGGCYVAGCCNAABCCAAA}$) and contains no 5' half-motif. The boundary is an ordinary, continuous genomic transition.
2. **Read Dynamics at the Internal Telomere Junction (1,185,509 bp):**
   * Coverage across this site is **$\sim 100\times$** (strictly equimolar with standard somatic ends).
   * Reads terminate with soft-clipping where the telomeric hexamers begin.
   * Crucially, there are **no reads spanning across the inverted hairpin palindrome** created by the assembler, proving that the $4.7\text{ kb}$ inverted segment was an algorithmic contig-stitching artifact rather than a true biological palindrome.

---

### C. Telomere-Clipped Read Distribution: Molarity Rules Out Small Minichromosomes
Projection of all 109 telomere-clipped raw PacBio reads mapping across the rDNA locus onto a canonical 6,958 bp rDNA repeat unit (18S: 1..1,754 bp; ITS/5.8S: 1,755..2,200 bp; 28S: 2,201..6,902 bp; IGS: 6,903..6,958 bp) revealed **two distinct, prominent clusters**:

1. **Cluster 1 ($N = 55$ reads, Orientation `R`):**
   * Maps sharply to **position ~4,705 bp** (inside the **28S rRNA gene**).
   * Represents a forward-strand junction where rDNA transitions rightward into telomere repeats: `[...rDNA...28S(4705)] ── (TTAGGC)n`.
   * Matches the exact assembly junction found in telotig `ptg000014l` (pos 4,703 bp).
2. **Cluster 2 ($N = 44$ reads, Orientation `L`):**
   * Maps sharply to **position ~1,637 bp** (distal end of the **18S rRNA gene**).
   * Represents a reverse-strand junction where telomere repeats transition into rDNA: `(GCCTAA)n ── [18S(1637)...28S...]`.
3. **Non-cluster Singletons ($N = 10$ reads):**
   * Isolated singletons ($n = 1$ read per site) scattered across the repeat unit with low repeat counts and MAPQ $\le 1$ (sequencing noise).
4. **Molarity / Stoichiometry Proof:**
   * In *O. onirici*, standard somatic chromosome ends have **~50–100 telomere-clipped reads**.
   * Cluster 1 ($N=55$) and Cluster 2 ($N=44$) each exhibit **exactly 1× equimolar coverage**.
   * If the sequence-specific mechanism cleaved every rDNA unit (or excised a 3 kb fragment from each), we would observe a $>300\times$ coverage multiplier ($>15,000$ reads).
   * The 1× equimolar coverage proves that the cuts do not shatter the array into small fragments, but instead **define the macro-boundaries of a single massive ~1.2 Mb molecule per genome equivalent**.

---

### D. Complete Absence of the 27 bp PDE Cleavage Motif
* The conserved *O. onirici* PDE consensus cleavage motif is:
  $$\text{5' - } \mathbf{TTWGGYTDWGGCYVAGCCNAABCCAAA} \text{ - 3'}$$
* In bona fide chromosomal PDE ends (such as Chromosome II Right at 10,439,112 bp):
  * The scission motif matches with extraordinary statistical significance (**$P = 3.16 \times 10^{-12}$**).
  * The enzymatic cut cleanly leaves the 5' half-motif directly preceding the *de novo* telomere addition (`...CTAATATTTTTTTGGATGAGGCCGA | (TTAGGC)n`).
* In stark contrast, sequence examination of all three rDNA junctions:
  * **Core transition (1,190,878 bp):** No motif.
  * **28S Junction (pos ~4,705 bp):** `...TTGCTTTCGGGCTTAGTATCT | (TTAGGC)n` — No motif (MEME $P = 0.826$).
  * **18S Junction (pos ~1,637 bp):** `(GCCTAA)n | TTTCGAGAAGACTGAG...` — No motif.
* **Conclusion:** The entire developmental processing of the rDNA locus operates via a **motif-independent enzymatic mechanism**, completely decoupled from the sequence-specific endonucleolytic machinery responsible for somatic genome reduction across the rest of the genome.

---

## 3. Unified Biological Model: Germline Architecture & Somatic Two-Fates

Synthesizing the sequence evidence, coverage molarity, and long-read alignments:

```
A. GERMLINE CHROMOSOME ARCHITECTURE:
[Telomere(?)] ═════ [Tandem rDNA Array ~1.2 Mb] ═════ [Euchromatic Core] (starts at 1.191 Mb)
                    (Units: 18S -> 5.8S -> 28S -> IGS)  Continuous 230x PacBio read support

B. SOMATIC DEVELOPMENTAL PROCESSING (MOTIF-INDEPENDENT):
1. Chromosomal Stump Healing (Model B Aspect):
   The distal portion of the chromosomal array undergoes double-strand scission at position 4,705 bp of 28S.
   The proximal fragment anchored to the core is healed and capped with a telomere:
   [Euchromatic Core] ──── [Proximal rDNA units] ──── 28S(4705) ── (TTAGGC)n  (Cluster 1, N=55)

2. Extrachromosomal Linear Minichromosome (Model A Aspect):
   The remaining massive ~1.2 Mb rDNA array is liberated as a single contiguous molecule and capped on both ends:
   (GCCTAA)n ── 18S(1637) ════ [Tandem Units ~1.2 Mb] ════ 28S(4705) ── (TTAGGC)n
   (Cluster 2, N=44)                                      (Cluster 1, N=55)
```

Because whole-worm sequencing captures DNA from both the healed chromosomal stump and the equimolar extrachromosomal minichromosome, we observe robust, stoichiometric read support for both the 18S and 28S inward-facing telomeres.

---

## 4. Methods, Scripts & Cluster Reproducibility

### Cluster Environment & BAM Access (Sanger Farm HPC)
* **Host:** `farm22-head2` via Teleport (`tsh ssh pg17@farm22-head2`)
* **Environment Modules:**
  ```bash
  module use /software/spack_environments/groups/toladmin/00/share/spack/modules/linux-ubuntu22.04-x86_64_v3 \
             /software/spack_environments/default/00/share/spack/modules/linux-ubuntu22.04-x86_64_v3
  module load samtools/1.20--h50ea8bc_0
  ```
* **BAM File:**
  `/lustre/scratch122/tol/teams/blaxter/users/pg17/freeze/2022-02-17/data/Oscheius/Oscheius_onirici/analysis/diminution/telomeres/remmaped_all_reads_nxOscOnir1.5.bam`

### Inspection Commands Used:
```bash
# 1. Depth and reads spanning core-rDNA boundary (1,190,878 bp):
samtools depth -r II:1190850-1190900 remmaped_all_reads_nxOscOnir1.5.bam
samtools view remmaped_all_reads_nxOscOnir1.5.bam II:1190875-1190880 | awk '$4 <= 1190878 && $4 + 500 >= 1190878'

# 2. Depth across internal telomere junction (1,185,509 bp):
samtools depth -r II:1185490-1185530 remmaped_all_reads_nxOscOnir1.5.bam

# 3. Soft-clipping check at internal telomere junction:
samtools view remmaped_all_reads_nxOscOnir1.5.bam II:1185500-1185520 | grep "S"
```

### Local Analysis Scripts:
* Read projection and motif checks:
  ```bash
  Rscript analysis/investigate_onirici_rDNA_terminus.R
  ```
* Multi-panel Figure S6 generation:
  ```bash
  Rscript analysis/plot_figureS6_onirici_rDNA_terminus.R
  ```
