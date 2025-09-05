#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
from typing import List, Tuple
import numpy as np

# Inputs: per-protein features from 6D
FEAT_ROOT = Path("data/outputs/structures/features")

# Outputs
OUT_ROOT = Path("data/outputs/segments")

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["orphanZ70_reps","strict26k_reps","orphanZ70_all"], default="orphanZ70_reps")
    ap.add_argument("--threshold", type=float, default=50.0, help="pLDDT cutoff for low-confidence")
    ap.add_argument("--min-run", type=int, default=3, help="min consecutive residues below threshold to keep")
    ap.add_argument("--limit", type=int, default=0, help="limit proteins (0=all)")
    return ap.parse_args()

def find_runs_below(plddt: np.ndarray, thr: float, min_run: int) -> List[Tuple[int,int]]:
    spans: List[Tuple[int,int]] = []
    N = int(plddt.size)
    i = 0
    while i < N:
        if plddt[i] < thr:
            j = i + 1
            while j < N and plddt[j] < thr:
                j += 1
            if (j - i) >= min_run:
                # 1-based residue indices
                spans.append((i + 1, j))  # inclusive
            i = j
        else:
            i += 1
    return spans

def main() -> None:
    args = parse_args()
    feat_dir = FEAT_ROOT / args.set
    idx = feat_dir / "features_index.tsv"
    if not idx.exists():
        raise SystemExit(f"Missing features index: {idx}. Run 6D first.")
    out_dir = OUT_ROOT / args.set
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / "lowconf_spans.tsv"
    out_json = out_dir / "lowconf_summary.json"

    rows = []
    with idx.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr: rows.append(r)
    if args.limit and args.limit < len(rows):
        rows = rows[:args.limit]

    n_proc = 0
    n_with_spans = 0
    total_spans = 0

    with out_tsv.open("w") as out:
        out.write("\t".join([
            "uniprot","n_res","n_low","low_frac","n_spans","spans","cut_sites","plddt_min","plddt_mean","plddt_max"
        ]) + "\n")
        for r in rows:
            acc = r["uniprot"]
            npz_path = Path(r["npz_path"])
            if not npz_path.exists():
                continue
            data = np.load(npz_path, allow_pickle=True)
            plddt = data["plddt"].astype(float)
            N = int(plddt.size)
            if N == 0:
                out.write(f"{acc}\t0\t0\t0.0000\t0\t\t\t0.00\t0.00\t0.00\n")
                n_proc += 1
                continue

            spans = find_runs_below(plddt, args.threshold, args.min_run)
            n_low = int((plddt < args.threshold).sum())
            low_frac = (n_low / N) if N > 0 else 0.0
            cut_sites = [int((s + e) // 2) for (s, e) in spans]
            out.write("\t".join([
                acc, str(N), str(n_low), f"{low_frac:.4f}", str(len(spans)),
                ";".join(f"{s}-{e}" for (s,e) in spans),
                ";".join(map(str, cut_sites)),
                f"{plddt.min():.2f}", f"{plddt.mean():.2f}", f"{plddt.max():.2f}"
            ]) + "\n")

            n_proc += 1
            if spans:
                n_with_spans += 1
                total_spans += len(spans)

    out_json.write_text(json.dumps({
        "set": args.set,
        "threshold": args.threshold,
        "min_run": args.min_run,
        "processed": n_proc,
        "proteins_with_lowconf_spans": n_with_spans,
        "total_lowconf_spans": total_spans
    }, indent=2))
    print(f"Wrote {out_tsv}")
    print(f"Wrote {out_json}")

if __name__ == "__main__":
    main()
