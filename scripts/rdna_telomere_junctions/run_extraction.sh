#!/usr/bin/env bash
# Reproduce data/rdna_telomere/*.tracts.tsv and *.region_reads.txt for Figure S6.
#
# Runs on a workstation with a Teleport (tsh) session to farm22. Only small region queries
# are issued on the head node; BAMs never leave the cluster (only telomere-bearing read
# records are streamed back and parsed locally by extract_junctions.py).
# Cluster conventions: see pde-curation-context/farm22/farm22_execution_rules.md
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${HERE}/../../data/rdna_telomere"
mkdir -p "${OUT}"

FREEZE=/lustre/scratch122/tol/teams/blaxter/users/pg17/freeze/2022-02-17/data/Oscheius
SIF=/software/tola/images/htstools-1.17.sif
TELO_RE='(TTAGGC){4}|(GCCTAA){4}'

# strain | assembly | species dir | rDNA array region on the assembly
while IFS='|' read -r strain asm dir region; do
  bam="${FREEZE}/${dir}/analysis/coverage/minimap/pb/${asm}.bam"
  echo "[${strain}] ${asm} ${region}" >&2
  # primary mapped reads overlapping the array (denominator for read-support normalisation)
  tsh ssh pg17@farm22-head1 "sexec ${SIF} samtools view -c -F 2308 ${bam} ${region}" </dev/null \
    > "${OUT}/${asm}.region_reads.txt"
  # every telomere-bearing read, parsed block-by-block
  tsh ssh pg17@farm22-head1 </dev/null \
    "sexec ${SIF} samtools view -F 2308 ${bam} ${region} | grep -E '${TELO_RE}'" \
    | python3 "${HERE}/extract_junctions.py" --species "${strain}" --gap-report "${OUT}/${asm}.block_gaps.tsv" \
        --flank-fasta "${OUT}/${asm}.flanks.fa" \
    > "${OUT}/${asm}.tracts.tsv"
done <<'EOF'
PS2068|nxOscOnir1.2|Oscheius_onirici|II:1-1260000
DF5120|nxOscSper1.1|Oscheius_sp._DF5120|II:8130000-8180000
PS1017|nxOscDolc1.1|Oscheius_dolichura|IV:12450000-12464839
EOF
