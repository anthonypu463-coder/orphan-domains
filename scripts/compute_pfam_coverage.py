#!/usr/bin/env python3
from __future__ import annotations
import csv
from pathlib import Path
from typing import Dict, List, Tuple

MANIFEST = Path("data/afdb/index/reference_manifest.tsv")           # has 'uniprot' and 'residues'
DOMS     = Path("data/outputs/hmmscan/pfam_domains_nonoverlap.tsv") # from 2C
OUT      = Path("data/outputs/hmmscan/pfam_coverage_nonoverlap.tsv")

def to_int(s: str) -> int:
    try: return int(s)
    except Exception: return 0

def load_lengths(p: Path) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    with p.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        # Expect columns: uniprot, residues, ...
        for row in rdr:
            acc = row.get("uniprot", "").strip()
            if not acc: 
                continue
            L = to_int(row.get("residues", "0"))
            lengths[acc] = L
    if not lengths:
        raise SystemExit("No lengths parsed from manifest")
    return lengths

def load_intervals(p: Path) -> Dict[str, List[Tuple[int,int]]]:
    spans: Dict[str, List[Tuple[int,int]]] = {}
    with p.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            acc = row.get("uniprot", "").strip()
            if not acc:
                continue
            a, b = to_int(row.get("env_from", "")), to_int(row.get("env_to", ""))
            if a <= 0 or b <= 0 or a > b:
                continue
            spans.setdefault(acc, []).append((a, b))
    return spans

def covered_length(intervals: List[Tuple[int,int]]) -> int:
    if not intervals:
        return 0
    # No adjacency merging: only true overlaps are coalesced.
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    total = 0
    cur_l, cur_r = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_r:                 # overlap
            cur_r = max(cur_r, e)
        else:                          # gap
            total += (cur_r - cur_l + 1)
            cur_l, cur_r = s, e
    total += (cur_r - cur_l + 1)
    return total

def main() -> None:
    lengths = load_lengths(MANIFEST)
    spans = load_intervals(DOMS)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as fh:
        fh.write("\t".join(["uniprot","residues","pfam_domains","covered_residues","coverage_frac"]) + "\n")
        for acc, L in lengths.items():
            iv = spans.get(acc, [])
            cov = covered_length(iv)
            frac = (cov / L) if L > 0 else 0.0
            # Count domains as rows in nonoverlap file for that protein
            fh.write(f"{acc}\t{L}\t{len(iv)}\t{cov}\t{frac:.4f}\n")
    print(f"Wrote {OUT}")

if __name__ == "__main__":
    main()
