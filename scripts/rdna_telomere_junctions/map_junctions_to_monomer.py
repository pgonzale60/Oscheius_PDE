#!/usr/bin/env python3
"""
Place telomere junctions on a single rDNA monomer using the flanking read sequence.

For every clipped telomere tract, extract_junctions.py writes the array-side flank (up to 1 kb)
adjacent to the tract. Each flank is aligned (minimap2, map-hifi) to its species' monomer,
duplicated head-to-tail so that junctions near the unit boundary are not cut. The junction is the
flank end that touched the telomere:

    clip_right tract : flank 3' end    clip_left tract : flank 5' end

Output per junction: monomer position (0-based, 18S start = 0, wrapped to the monomer length)
and the side of the junction the telomere lies on, in monomer-plus-strand terms:
    downstream : telomere continues toward higher monomer coordinates (monomer + strand runs into it)
    upstream   : telomere lies at lower coordinates (monomer - strand runs into it)
This is independent of assembly orientation, so inverted or collapsed copies in the assembly do
not matter.
"""
import subprocess
import sys
import tempfile
from collections import OrderedDict


def read_fa(path):
    d, name = OrderedDict(), None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:]
            d[name] = []
        else:
            d[name].append(line)
    return OrderedDict((k, "".join(v)) for k, v in d.items())


def main(monomers_fa, flanks_fa, strain, out_tsv, max_edge_clip=50):
    mono = {k.split("|")[0]: v for k, v in read_fa(monomers_fa).items()}[strain]
    L = len(mono)
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as t:
        t.write(f">{strain}_x2\n{mono + mono}\n")
        ref = t.name
    flanks = read_fa(flanks_fa)
    paf = subprocess.run(
        ["minimap2", "-c", "-x", "map-hifi", "--secondary=no", "-t", "4", ref, flanks_fa],
        capture_output=True, text=True, check=True).stdout
    best = {}
    for line in paf.splitlines():
        f = line.split("\t")
        q, ql, qs, qe, strand = f[0], int(f[1]), int(f[2]), int(f[3]), f[4]
        ts, te, matches, blk, mapq = int(f[7]), int(f[8]), int(f[9]), int(f[10]), int(f[11])
        if q not in best or blk > best[q][-1]:
            best[q] = (ql, qs, qe, strand, ts, te, matches, blk, mapq, blk)
    with open(out_tsv, "w") as o:
        o.write("\t".join(["strain", "read", "tract_idx", "location", "motif", "flank_len", "mapped",
                           "strand", "edge_clip", "identity", "mapq", "mono_pos", "mono_phase",
                           "telomere_side"]) + "\n")
        for name, seq in flanks.items():
            sp_, rid, idx, loc, dom = name.split("|")
            ql = len(seq)
            if name not in best:
                o.write("\t".join([sp_, rid, idx, loc, dom, str(ql), "0"] + ["NA"] * 7) + "\n")
                continue
            _, qs, qe, strand, ts, te, mt, blk, mapq, _ = best[name]
            right = loc == "clip_right"               # junction at flank 3' end
            edge_clip = (ql - qe) if right else qs
            if right:
                pos = te + edge_clip if strand == "+" else ts - edge_clip
            else:
                pos = ts - edge_clip if strand == "+" else te + edge_clip
            side = "downstream" if ((right and strand == "+") or (not right and strand == "-")) else "upstream"
            pos = pos % L
            o.write("\t".join(map(str, [sp_, rid, idx, loc, dom, ql, 1, strand, edge_clip,
                                        round(mt / blk, 4), mapq, pos, round(pos / L, 5), side])) + "\n")


if __name__ == "__main__":
    main(*sys.argv[1:5])
