#!/usr/bin/env python3
"""Filter or lift an ancestry-frequency panel and validate target REF alleles.

Input columns: chrom, pos, ref, alt, followed by population ALT frequencies.
"""

from __future__ import annotations

import argparse
import csv

import pysam
from pyliftover import LiftOver

BASE_COLUMNS = ("chrom", "pos", "ref", "alt")
AUTOSOMES = {f"chr{i}" for i in range(1, 23)}
COMP = str.maketrans("ACGT", "TGCA")


def chrom(value: str) -> str:
    return value if value.startswith("chr") else "chr" + value


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--reference", required=True, help="target-assembly FASTA")
    p.add_argument("--chain", help="source-to-target chain; omit when input already uses the target assembly")
    p.add_argument("--reverse-chain", help="target-to-source chain for reciprocal validation")
    p.add_argument("--min-af", type=float, default=0.005)
    p.add_argument("--max-af", type=float, default=0.995)
    p.add_argument("--min-spread", type=float, default=0.20,
                   help="minimum max-minus-min population frequency")
    p.add_argument("--min-distance", type=int, default=100)
    args = p.parse_args()
    if args.chain and not args.reverse_chain:
        p.error("--reverse-chain is required when --chain is used")
    lo = LiftOver(args.chain) if args.chain else None
    rev = LiftOver(args.reverse_chain) if args.reverse_chain else None
    fasta = pysam.FastaFile(args.reference)
    kept = []
    with open(args.input) as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        missing = set(BASE_COLUMNS) - set(rows.fieldnames or ())
        if missing:
            raise SystemExit(f"Panel missing columns: {', '.join(sorted(missing))}")
        populations = tuple(x for x in (rows.fieldnames or ()) if x not in BASE_COLUMNS)
        if len(populations) < 2:
            raise SystemExit("Panel must have at least two population-frequency columns")
        for row in rows:
            try:
                c, pos = chrom(row["chrom"]), int(row["pos"])
                ref, alt = row["ref"].upper(), row["alt"].upper()
                freqs = [float(row[x]) for x in populations]
            except (KeyError, TypeError, ValueError):
                continue
            if c not in AUTOSOMES or len(ref) != 1 or len(alt) != 1 or ref == alt or any(not 0 <= x <= 1 for x in freqs):
                continue
            if max(freqs) < args.min_af or min(freqs) > args.max_af or max(freqs) - min(freqs) < args.min_spread:
                continue
            if lo:
                hits = lo.convert_coordinate(c, pos - 1)
                if len(hits) != 1:
                    continue
                tc, tp, strand, _ = hits[0]
                if tc != c:
                    continue
                if rev:
                    back = rev.convert_coordinate(tc, int(tp))
                    if len(back) != 1 or chrom(back[0][0]) != c or int(back[0][1]) != pos - 1:
                        continue
                if strand == "-":
                    ref, alt = ref.translate(COMP), alt.translate(COMP)
                c, pos = tc, int(tp) + 1
            target_ref = fasta.fetch(c, pos - 1, pos).upper()
            if target_ref == alt:
                ref, alt = alt, ref
                freqs = [1.0 - value for value in freqs]
            elif target_ref != ref:
                continue
            kept.append((c, pos, ref, alt, freqs))
    kept.sort(key=lambda x: (int(x[0][3:]), x[1]))
    selected, last = [], {}
    seen = set()
    for row in kept:
        key = (row[0], row[1])
        if key in seen:
            continue
        if row[1] - last.get(row[0], -10**12) < args.min_distance:
            continue
        selected.append(row); last[row[0]] = row[1]; seen.add(key)
    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["chrom", "pos", "ref", "alt", *populations])
        for c, pos, ref, alt, freqs in selected:
            writer.writerow([c, pos, ref, alt, *[f"{x:.8g}" for x in freqs]])
    print(f"Retained {len(selected)} of {len(kept)} validated markers")


if __name__ == "__main__":
    main()
