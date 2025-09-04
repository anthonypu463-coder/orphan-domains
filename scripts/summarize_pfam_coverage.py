#!/usr/bin/env python3
from collections import defaultdict
from pathlib import Path
import csv

HITS = Path("data/outputs/hmmscan/pfam_hits.tsv")
MANI = Path("data/afdb/index/reference_manifest.tsv")
OUT  = Path("data/outputs/hmmscan/pfam_coverage.tsv")

def load_lengths(manifest: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with manifest.open() as fh:
        next(fh)  # header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[2].isdigit():
                lengths[parts[0]] = int(parts[2])
    return lengths

def main():
    if not HITS.exists():
        raise SystemExit(f"Missing {HITS}")
    if not MANI.exists():
        raise SystemExit(f"Missing {MANI}")
    lengths = load_lengths(MANI)
    covered = defaultdict(int)
    with HITS.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            acc = row["uniprot"]
            try:
                covered[acc] += int(row["dom_len"])
            except Exception:
                pass
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as fh:
        fh.write("uniprot\tresidues\tpfam_covered\tcoverage_frac\n")
        for acc, L in lengths.items():
            cov = covered.get(acc, 0)
            frac = (cov / L) if L > 0 else 0.0
            fh.write(f"{acc}\t{L}\t{cov}\t{frac:.4f}\n")
    print(f"Wrote {OUT}")

if __name__ == "__main__":
    main()
