#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
from typing import List, Tuple, Optional, Dict
import numpy as np

# Inputs
SEG_ROOT = Path("data/outputs/segments")                     # has <SET>/lowconf_spans.tsv from 7A
FEAT_ROOT = Path("data/outputs/structures/features")         # has <SET>/by_acc/*.npz and features_index.tsv

# Outputs
OUT_ROOT = Path("data/outputs/segments")

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["orphanZ70_reps","strict26k_reps","orphanZ70_all"], default="orphanZ70_reps")
    ap.add_argument("--min-fragment", type=int, default=50, help="keep fragments with length >= this")
    ap.add_argument("--limit", type=int, default=0, help="limit proteins (0=all)")
    return ap.parse_args()

def merge_adjacent_spans(spans: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    if not spans: return []
    spans = sorted(spans)
    out = [list(spans[0])]
    for s,e in spans[1:]:
        ps,pe = out[-1]
        if s <= pe + 1:          # adjacent or overlapping
            out[-1][1] = max(pe, e)
        else:
            out.append([s,e])
    return [(s,e) for s,e in out]

def cut_positions(spans: List[Tuple[int,int]]) -> List[int]:
    # midpoint (floor) per merged span
    return [ (s + e) // 2 for (s,e) in spans ]

def compute_fragment_stats(plddt: Optional[np.ndarray], s: int, e: int) -> Dict[str, float]:
    if plddt is None or plddt.size == 0:
        return {"plddt_min": 0.0, "plddt_mean": 0.0, "plddt_max": 0.0, "n_low_lt50": 0}
    # indices are 1-based in spans
    frag = plddt[s-1:e]
    return {
        "plddt_min": float(frag.min()),
        "plddt_mean": float(frag.mean()),
        "plddt_max": float(frag.max()),
        "n_low_lt50": int((frag < 50.0).sum()),
    }

def main() -> None:
    args = parse_args()
    set_dir = OUT_ROOT / args.set
    spans_tsv = set_dir / "lowconf_spans.tsv"
    if not spans_tsv.exists():
        raise SystemExit(f"Missing lowconf spans: {spans_tsv}. Run 7A first.")

    feat_dir = FEAT_ROOT / args.set
    npz_dir = feat_dir / "by_acc"
    out_dir = OUT_ROOT / args.set
    out_dir.mkdir(parents=True, exist_ok=True)
    out_frag = out_dir / "segments_by_plddt.tsv"
    out_json = out_dir / "segments_summary.json"

    # Load rows
    rows = []
    with spans_tsv.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr: rows.append(r)
    if args.limit and args.limit < len(rows):
        rows = rows[:args.limit]

    kept_frag = 0
    drop_short = 0
    proteins = 0
    with out_frag.open("w") as out:
        out.write("\t".join([
            "uniprot","n_res","fragment_id","start","end","length",
            "plddt_min","plddt_mean","plddt_max","n_low_lt50","kept"
        ]) + "\n")

        for r in rows:
            acc = r["uniprot"]
            n_res = int(r["n_res"] or 0)
            # parse spans like "85-95;210-220"
            spans_str = r["spans"].strip()
            spans: List[Tuple[int,int]] = []
            if spans_str:
                for tok in spans_str.split(";"):
                    if not tok: continue
                    a,b = tok.split("-")
                    s,e = int(a), int(b)
                    if 1 <= s <= e <= n_res:
                        spans.append((s,e))
            spans = merge_adjacent_spans(spans)

            # compute cuts and fragments
            cuts = cut_positions(spans)
            cuts = sorted(set(c for c in cuts if 1 <= c < n_res))

            # attempt to load pLDDT for stats
            plddt = None
            npz_path = npz_dir / f"{acc}.npz"
            if npz_path.exists():
                try:
                    data = np.load(npz_path, allow_pickle=True)
                    plddt = data["plddt"].astype(float)
                except Exception:
                    plddt = None

            proteins += 1
            # generate fragments
            frag_id = 0
            start = 1
            all_cuts = cuts + [n_res]  # sentinel to flush last fragment
            for c in all_cuts:
                end = c if c == n_res else c
                if end < start:   # guard
                    continue
                # if this is sentinel, end==n_res; otherwise, next start is c+1
                length = end - start + 1
                stats = compute_fragment_stats(plddt, start, end)
                keep = int(length >= args.min_fragment)
                if keep: kept_frag += 1
                else: drop_short += 1
                frag_id += 1
                out.write("\t".join([
                    acc, str(n_res), f"F{frag_id:03d}",
                    str(start), str(end), str(length),
                    f"{stats['plddt_min']:.2f}", f"{stats['plddt_mean']:.2f}", f"{stats['plddt_max']:.2f}",
                    str(stats["n_low_lt50"]), str(keep)
                ]) + "\n")
                start = c + 1
            # If no cuts at all, the loop above creates exactly one fragment [1,n_res]

    out_json.write_text(json.dumps({
        "set": args.set,
        "min_fragment": args.min_fragment,
        "proteins": proteins,
        "fragments_kept": kept_frag,
        "fragments_dropped_short": drop_short,
        "fragments_total": kept_frag + drop_short
    }, indent=2))
    print(f"Wrote {out_frag}")
    print(f"Wrote {out_json}")

if __name__ == "__main__":
    main()
