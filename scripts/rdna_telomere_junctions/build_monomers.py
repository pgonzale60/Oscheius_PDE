#!/usr/bin/env python3
"""
Build one rDNA monomer per species from the assembled array, with 18S start as position 0.

The monomer is a reference for placing telomere junctions on a single repeat unit,
independent of how many copies (or in which orientation) each assembly carries.

  * array with >= 2 annotated 18S: the 18S-start to 18S-start unit whose length is closest to
    the species median spacing.
  * collapsed array with a single 18S (PS1017): 28S-start to next 28S-start unit, rotated so that
    the 18S start is position 0 (the monomer is treated as circular downstream).
"""
import gzip
import re
import statistics as st
import sys

SPECIES = [  # strain, assembly, rDNA chromosome, array window (lo, hi)
    ("PS2068", "nxOscOnir1.2", "II", 1, 1260000),
    ("DF5120", "nxOscSper1.1", "II", 8130000, 8180000),
    ("PS1017", "nxOscDolc1.1", "IV", 12450000, 12464839),
]


def load_contig(path, name):
    seq, on = [], False
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                if on:
                    break
                on = line[1:].split()[0] == name
            elif on:
                seq.append(line.strip())
    return "".join(seq).upper()


def rrna(gff, chrom, lo, hi):
    out = []
    for line in open(gff):
        if line[0] == "#":
            continue
        p = line.rstrip("\n").split("\t")
        m = re.match(r"ID=(\d+)s_", p[8])
        if p[0] == chrom and lo <= int(p[3]) <= hi and m and m.group(1) in ("18", "28"):
            out.append((int(p[3]), int(p[4]), m.group(1)))
    return sorted(out)


def main(repo, out_fa, out_tsv):
    fa = open(out_fa, "w")
    meta = open(out_tsv, "w")
    meta.write("strain\tassembly\tcontig\tstart\tend\tlength\trotation\tn18S\tmedian_spacing\n")
    for strain, asm, chrom, lo, hi in SPECIES:
        r = rrna(f"{repo}/data/rnammer/{asm}.rnammer.gff3", chrom, lo, hi)
        s18 = [x[0] for x in r if x[2] == "18"]
        s28 = [x[0] for x in r if x[2] == "28"]
        contig = load_contig(f"{repo}/data/sequence_sizes/{asm}.primary.fa.gz", chrom)
        if len(s18) >= 2:
            sp = [b - a for a, b in zip(s18, s18[1:])]
            med = st.median(sp)
            i = min(range(len(sp)), key=lambda k: abs(sp[k] - med))
            start, end, rot = s18[i], s18[i + 1] - 1, 0
            mono = contig[start - 1:end]
        else:
            med = "NA"
            start, end = s28[0], s28[1] - 1                # 28S start to next 28S start
            unit = contig[start - 1:end]
            rot = s18[0] - start                           # 18S start inside that unit
            mono = unit[rot:] + unit[:rot]
        fa.write(f">{strain}|{asm}|{chrom}:{start}-{end}|rot{rot}\n{mono}\n")
        meta.write(f"{strain}\t{asm}\t{chrom}\t{start}\t{end}\t{len(mono)}\t{rot}\t{len(s18)}\t{med}\n")


if __name__ == "__main__":
    main(*sys.argv[1:4])
