#!/usr/bin/env python3
from __future__ import annotations
import csv, json, math
from pathlib import Path
from typing import Optional

IN_TSV   = Path("data/outputs/architecture/pfam_architecture_enhanced.tsv")
OUT_DIR  = Path("data/outputs/orphans")
OUT_TSV  = OUT_DIR / "orphan_candidates.tsv"
OUT_LIST = OUT_DIR / "orphan_accessions.txt"
OUT_SUM  = OUT_DIR / "summary.json"
THRESH   = 0.20  # locked threshold for orphan-rich selection

def to_int(s: str) -> int:
    try: return int(s)
    except Exception: return 0

def to_float(s: str) -> float:
    try: return float(s)
    except Exception: return 0.0

def largest_gap_len(gap_lengths: str) -> int:
    if not gap_lengths: return 0
    best = 0
    for tok in gap_lengths.split(";"):
        tok = tok.strip()
        if not tok: continue
        try:
            best = max(best, int(tok))
        except Exception:
            pass
    return best

def main() -> None:
    if not IN_TSV.exists():
        raise SystemExit(f"Missing {IN_TSV}; run 3D first.")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    total = 0
    with IN_TSV.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            total += 1
            cov = to_float(r.get("coverage_frac","0"))
            if cov < THRESH:
                L  = to_int(r.get("residues","0"))
                cov_res = to_int(r.get("covered_residues","0"))
                n_dom = to_int(r.get("n_domains","0"))
                gaps = r.get("gaps","")
                gap_lengths = r.get("gap_lengths","")
                rows.append({
                    "uniprot": r["uniprot"],
                    "residues": L,
                    "covered_residues": cov_res,
                    "coverage_frac": cov,
                    "unannotated_residues": max(L - cov_res, 0),
                    "unannotated_frac": max(1.0 - cov, 0.0),
                    "n_domains": n_dom,
                    "largest_gap": largest_gap_len(gap_lengths),
                    "gaps": gaps,
                    "gap_lengths": gap_lengths,
                })

    # Sort by lowest coverage (most orphan-like first), then by largest gap
    rows.sort(key=lambda x: (x["coverage_frac"], -x["largest_gap"], -x["unannotated_residues"]))

    # Write TSV
    with OUT_TSV.open("w") as fh:
        fh.write("\t".join([
            "uniprot","residues","covered_residues","coverage_frac",
            "unannotated_residues","unannotated_frac","n_domains","largest_gap","gaps","gap_lengths"
        ]) + "\n")
        for r in rows:
            fh.write("\t".join([
                r["uniprot"],
                str(r["residues"]),
                str(r["covered_residues"]),
                f"{r['coverage_frac']:.4f}",
                str(r["unannotated_residues"]),
                f"{r['unannotated_frac']:.4f}",
                str(r["n_domains"]),
                str(r["largest_gap"]),
                r["gaps"],
                r["gap_lengths"],
            ]) + "\n")

    # Write accession list
    with OUT_LIST.open("w") as fh:
        for r in rows:
            fh.write(r["uniprot"] + "\n")

    # Summary
    n = len(rows)
    mean_cov = (sum(r["coverage_frac"] for r in rows) / n) if n else 0.0
    OUT_SUM.write_text(json.dumps({
        "threshold": THRESH,
        "total_proteins": total,
        "orphan_candidates": n,
        "mean_coverage_in_candidates": round(mean_cov, 6),
    }, indent=2))

    print(f"Wrote {OUT_TSV}")
    print(f"Wrote {OUT_LIST}")
    print(f"Wrote {OUT_SUM}")

if __name__ == "__main__":
    main()
