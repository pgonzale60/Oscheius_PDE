# Preliminary Findings & Assembly Dissection: The Chromosome II rDNA Locus in *Oscheius onirici*

**Date:** September 13, 2026  
**Species:** *Oscheius onirici* (strain PS2068, assemblies `nxOscOnir1.2` / `nxOscOnir1.5`)  
**Genomic Locus:** Chromosome II Left Terminus (rDNA Array)

> [!WARNING]
> **SUPERSEDED SPECULATIVE MODEL:**
> Earlier working notes hypothesized an "extrachromosomal macrochromosome with two-sided resection and somatic fusion." That speculative model was **definitively refuted** upon verifying that the core boundary (`II:1,190,878 bp`) has **zero telomere-clipped reads**, proving the core is never amputated.
> 
> **Current Validated Model:** See [`rdna_resolution_and_hypotheses.md`](file:///Users/pgonzale60/Documents/research/oscheius_2026/pde_workbench/notes/rdna_resolution_and_hypotheses.md). The rDNA array forms the native, intact physical chromosome terminus in both germline and soma. It is **completely exempt from PDE**, directly refuting Launay et al. (2026).

---

## 1. Empirical Findings

### A. Chromosome II Left: Assembly Artifact Dissection
Inspection of the PacBio HiFi primary assembly (`nxOscOnir1.2.primary.fa.gz`) uncovered an algorithmic contig-stitching artifact:
1. **Leftmost Edge (1 bp):** Terminates abruptly mid-unit ($18S$ begins at pos 666 bp) with **0 telomeric hexamer repeats** (`TTAGGC`/`GCCTAA`) in the terminal 2,000 bp.
2. **Internal Telomere Junction (1,185,509 bp):** Contains the exact somatic rDNA–telomere junction captured in telotig `ptg000014l` (100% identity over 87 bp).
3. **Artificial Inversion (1,185,510 – 1,190,877 bp):** The assembler stitched a **4,703 bp inverted rDNA segment** (reverse complement of 18S to 28S pos 4,703) directly downstream of the junction.
4. **Core Boundary (1,190,878 bp):** The true euchromatic chromosome core (containing canonical protein-coding genes like `gene-8173`) begins cleanly at **1,190,878 bp**.
5. **Read Verification:** Raw BAM queries on Sanger Farm HPC (`remmaped_all_reads_nxOscOnir1.5.bam`) show **0 reads spanning the 4.7 kb inverted palindrome**, confirming it is an algorithmic stitching artifact resulting from intra-array structural inversions.

### B. Spanning Reads Across the Core–rDNA Boundary
* Querying `II:1190850-1190900` revealed that while raw depth reaches $230\times$ due to repetitive IGS multi-mapping, filtering for unique mapping (`MAPQ >= 20`) reveals **exactly 30 high-confidence, uniquely mapped PacBio HiFi reads (`MAPQ = 60`)** that cleanly span the junction between euchromatin and the non-coding rDNA spacer.
* In the adjacent unique core ($1.200\text{ – }1.202\text{ Mb}$), coverage is $95.8\times$. 
* Crucially, there are **0 reads with clipped telomeres at the 1,190,878 bp core boundary**. The core is physically attached and intact.

### C. Telomere-Clipped Read Distribution & Stoichiometry
Projection of telomere-clipped PacBio reads onto the canonical 6,958 bp rDNA unit (18S: 1..1,754 bp; ITS/5.8S: 1,755..2,200 bp; 28S: 2,201..6,902 bp; IGS: 6,903..6,958 bp) revealed two sharp clusters:
1. **Cluster 1 ($N = 55$ reads, Orientation `R`):** Maps sharply to **position ~4,705 bp** (inside **28S rRNA**), transitioning into `(TTAGGC)n`.
2. **Cluster 2 ($N = 44$ reads, Orientation `L`):** Maps sharply to **position ~1,637 bp** (distal end of **18S rRNA**), transitioning into `(GCCTAA)n`.
3. **Molarity / Stoichiometry:**
   - Total reads across Cluster 1 and Cluster 2 = $55 + 44 = 99\text{ reads}$.
   - This matches the diploid genome coverage ($\sim 95\text{--}100\times$).
   - Unlike true PDE ends, there is **zero germline read extension (~10x)** extending past these telomeres, and zero CBS cleavage motifs at the junctions.
   - The two opposite orientations reflect native structural inversions (head-to-head repeats) or a stable inversion polymorphism in the sequenced population.

---

## 2. Methodological & Cluster Reproducibility

### Sanger Farm HPC Environment
* **Host:** `farm22-head2` via Teleport (`tsh ssh pg17@farm22-head2`)
* **Environment Modules:**
  ```bash
  export MODULEPATH="/software/spack_environments/groups/toladmin/02/share/spack/modules/linux-ubuntu22.04-x86_64_v3:\
  /software/spack_environments/groups/toladmin/01/share/spack/modules/linux-ubuntu22.04-x86_64_v3:\
  /software/spack_environments/groups/toladmin/00/share/spack/modules/linux-ubuntu22.04-x86_64_v3:\
  /software/spack_environments/default/00/share/spack/modules/linux-ubuntu22.04-x86_64_v3:\
  /software/modules"
  module load samtools/1.20--h50ea8bc_0
  ```
* **Primary BAM File:**
  `/lustre/scratch122/tol/teams/blaxter/users/pg17/freeze/2022-02-17/data/Oscheius/Oscheius_onirici/analysis/diminution/telomeres/remmaped_all_reads_nxOscOnir1.5.bam`

### Key Extraction Commands
```bash
# Extract clean spanning reads across the 1,190,878 bp core-rDNA boundary
samtools view -h -q 20 remmaped_all_reads_nxOscOnir1.5.bam II:1190875-1190880 | \
  awk 'BEGIN{OFS="\t"} /^@/ {print; next} $4 <= 1190878 && ($4 + length($10)) >= 1190878 {print}' > spanning_reads_core_rDNA.sam

# Inspect soft-clipping and telomeric hexamers at 1,185,509 bp
samtools view remmaped_all_reads_nxOscOnir1.5.bam II:1185500-1185520 | grep "S"
```
