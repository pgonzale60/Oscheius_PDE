#!/usr/bin/env python3
"""
Extract telomeric tracts carried by long reads mapped to an rDNA array (SAM text on stdin).

Input  : `samtools view <bam> <region>` records on stdin, pre-filtered with
         `grep -E '(TTAGGC){4}|(GCCTAA){4}'` so only telomere-bearing reads arrive.
Output : one TSV row per telomeric TRACT. No read names are written (short hash only).

A HiFi telomere is not a perfect (TTAGGC)n string: indels/substitutions interrupt it, so a
single physical tract appears as many short exact blocks. Blocks of >= MIN_REPS hexamers that
are separated by <= MAX_GAP bp are therefore merged into one tract. Gaps between blocks are
reported before merging (see --gap-report) so the choice of MAX_GAP can be audited.

Each tract is placed relative to the read alignment:
  clip_right / clip_left : the tract touches the soft-clipped tail, i.e. a junction between
                           the array and a telomere. `ref_pos` is the alignment end/start.
  internal               : the tract lies wholly inside the aligned part (interstitial tract).
"""
import argparse
import hashlib
import re
import sys

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def walk(ops, aln_start):
    segs, rpos, qpos = [], aln_start, 0
    for n, op in ops:
        if op in "M=X":
            segs.append((qpos, qpos + n, rpos))
            qpos += n
            rpos += n
        elif op in "IS":
            qpos += n
        elif op in "DN":
            rpos += n
    return segs, rpos - 1          # last aligned reference base (1-based)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--min-reps", type=int, default=4)
    ap.add_argument("--max-gap", type=int, default=250)
    ap.add_argument("--gap-report", help="write raw inter-block gaps here for auditing")
    ap.add_argument("--flank-fasta", help="write the non-telomeric flank of each clipped tract here")
    ap.add_argument("--flank-len", type=int, default=1000)
    ap.add_argument("--min-flank", type=int, default=150)
    a = ap.parse_args()
    flank_out = open(a.flank_fasta, "w") if a.flank_fasta else None
    pat = re.compile(r"(?:TTAGGC){%d,}|(?:GCCTAA){%d,}" % (a.min_reps, a.min_reps))
    gaps_out = open(a.gap_report, "w") if a.gap_report else None
    if gaps_out:
        gaps_out.write("species\tgap\n")

    print("\t".join([
        "species", "read", "mapq", "read_len", "aln_start", "aln_end", "clip_left", "clip_right",
        "tract_idx", "n_tracts_read", "motif_fwd_reps", "motif_rev_reps", "tract_len",
        "n_blocks_merged", "location", "ref_pos", "offset_from_clip", "gap_to_next_tract"]))

    for line in sys.stdin:
        f = line.rstrip("\n").split("\t")
        flag, mapq, cigar, seq = int(f[1]), int(f[4]), f[5], f[9]
        if flag & 2304 or cigar == "*" or seq == "*":
            continue
        ops = [(int(n), op) for n, op in CIGAR.findall(cigar)]
        aln_start = int(f[3])
        segs, aln_end = walk(ops, aln_start)
        clip_l = ops[0][0] if ops[0][1] == "S" else 0
        clip_r = ops[-1][0] if ops[-1][1] == "S" else 0
        rl = len(seq)
        blocks = [(m.start(), m.end(), m.group(0)[:6]) for m in pat.finditer(seq)]
        if not blocks:
            continue
        # merge blocks into tracts
        tracts = []
        for s, e, mot in blocks:
            if tracts and s - tracts[-1]["e"] <= a.max_gap:
                t = tracts[-1]
                if gaps_out:
                    gaps_out.write(f"{a.species}\t{s - t['e']}\n")
                t["e"] = e
                t["n"] += 1
            else:
                tracts.append({"s": s, "e": e, "n": 1, "fwd": 0, "rev": 0})
                t = tracts[-1]
            t["fwd" if mot == "TTAGGC" else "rev"] += (e - s) // 6
        rid = hashlib.md5(f[0].encode()).hexdigest()[:8]
        for i, t in enumerate(tracts):
            s, e = t["s"], t["e"]
            if e >= rl - clip_r and clip_r > 0:
                loc, ref_pos, off = "clip_right", aln_end, s - (rl - clip_r)
            elif s <= clip_l and clip_l > 0:
                loc, ref_pos, off = "clip_left", aln_start, clip_l - e
            else:
                loc, off, ref_pos = "internal", "NA", "NA"
                for qs, qe, rs in segs:
                    if qs <= s < qe:
                        ref_pos = rs + (s - qs)
                        break
            nxt = tracts[i + 1]["s"] - e if i + 1 < len(tracts) else "NA"
            if flank_out and loc in ("clip_right", "clip_left"):
                # flank = array-side sequence adjacent to the tract, in read (reference) orientation
                prev_e = tracts[i - 1]["e"] if i > 0 else 0
                nxt_s = tracts[i + 1]["s"] if i + 1 < len(tracts) else rl
                fl = seq[max(prev_e, s - a.flank_len):s] if loc == "clip_right" else seq[e:min(nxt_s, e + a.flank_len)]
                if len(fl) >= a.min_flank:
                    dom = "F" if t["fwd"] > t["rev"] else "R"
                    flank_out.write(f">{a.species}|{rid}|{i + 1}|{loc}|{dom}\n{fl}\n")
            print("\t".join(map(str, [
                a.species, rid, mapq, rl, aln_start, aln_end, clip_l, clip_r, i + 1, len(tracts),
                t["fwd"], t["rev"], e - s, t["n"], loc, ref_pos, off, nxt])))


if __name__ == "__main__":
    main()
