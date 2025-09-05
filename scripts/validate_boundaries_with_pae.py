#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, os
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np

SEG_ROOT   = Path("data/outputs/segments")             # segments_final.tsv (7D)
STRUCT_DIR = Path("data/structures")                   # staged by 6A: by_acc/<ACC>/pae.json

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["orphanZ70_reps","strict26k_reps","orphanZ70_all"], default="orphanZ70_reps")
    ap.add_argument("--win", type=int, default=10, help="interface window size")
    ap.add_argument("--limit", type=int, default=0, help="limit number of proteins (0=all)")
    ap.add_argument("--high", type=float, default=15.0, help="PAE>this supports a split")
    ap.add_argument("--low", type=float, default=5.0, help="PAE<this suggests merge")
    return ap.parse_args()

def load_segments(set_name: str) -> Dict[str, List[Tuple[int,int,str]]]:
    """Return acc -> list of (start,end,frag_id) ordered by start."""
    seg_f = SEG_ROOT / set_name / "segments_final.tsv"
    if not seg_f.exists():
        raise SystemExit(f"Missing {seg_f} (run 7D first).")
    by: Dict[str, List[Tuple[int,int,str]]] = {}
    with seg_f.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            acc = r["uniprot"]
            s, e = int(r["start"]), int(r["end"])
            fid  = r["fragment_id"]
            by.setdefault(acc, []).append((s,e,fid))
    for acc in by:
        by[acc].sort(key=lambda x: (x[0], x[1]))
    return by

def staged_pae_path(set_name: str, acc: str) -> Path:
    # data/structures/<SET>/by_acc/<ACC>/pae.json
    return STRUCT_DIR / set_name / "by_acc" / acc / "pae.json"

def load_pae_matrix(p: Path) -> Optional[np.ndarray]:
    if not p.exists(): return None
    try:
        data = json.loads(p.read_text())
    except Exception:
        return None
    # Try known keys; otherwise search
    mat = None
    for k in ("predicted_aligned_error","pae","pae_matrix"):
        v = data.get(k) if isinstance(data, dict) else None
        if isinstance(v, list) and v and isinstance(v[0], list):
            mat = v; break
    if mat is None:
        # brute search for first square numeric matrix
        def find_square(o, depth=0):
            if depth>3: return None
            if isinstance(o, list) and o and isinstance(o[0], list):
                n=len(o)
                if all(isinstance(r, list) and len(r)==n for r in o):
                    return o
            if isinstance(o, dict):
                for vv in o.values():
                    x=find_square(vv, depth+1)
                    if x is not None: return x
            if isinstance(o, list):
                for vv in o:
                    x=find_square(vv, depth+1)
                    if x is not None: return x
            return None
        mat = find_square(data)
    if mat is None: return None
    try:
        arr = np.array(mat, dtype=float)
        if arr.ndim==2 and arr.shape[0]==arr.shape[1]:
            return arr
    except Exception:
        return None
    return None

def clamp(a:int, lo:int, hi:int)->int:
    return max(lo, min(hi, a))

def interface_score(pae: np.ndarray, a_end: int, b_start: int, win: int) -> Tuple[float,float,int,float]:
    """Average (and median) PAE between A's last win residues and B's first win residues (1-based inputs)."""
    n = pae.shape[0]
    a_lo = clamp(a_end - win + 1, 1, n)
    a_hi = clamp(a_end, 1, n)
    b_lo = clamp(b_start, 1, n)
    b_hi = clamp(b_start + win - 1, 1, n)
    # convert to 0-based slices
    A = np.arange(a_lo-1, a_hi, dtype=int)
    B = np.arange(b_lo-1, b_hi, dtype=int)
    if A.size==0 or B.size==0:
        return (float("nan"), float("nan"), 0, float("nan"))
    sub = pae[np.ix_(A,B)]
    avg = float(np.nanmean(sub))
    med = float(np.nanmedian(sub))
    n_pairs = int(sub.size)
    frac_gt15 = float(np.mean(sub > 15.0))
    return (avg, med, n_pairs, frac_gt15)

def classify(avg: float, low: float, high: float) -> str:
    if not np.isfinite(avg): return "no_pae"
    if avg < low:  return "suggest_merge"
    if avg > high: return "support_split"
    return "ambiguous"

def main() -> None:
    args = parse_args()
    segs = load_segments(args.set)
    if args.limit and args.limit < len(segs):
        keep = dict(sorted(segs.items())[:args.limit])
        segs = keep

    out_dir = SEG_ROOT / args.set
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / "pae_interface_checks.tsv"
    out_js  = out_dir / "pae_interface_summary.json"

    proteins = 0
    checks = 0
    c_support = c_merge = c_amb = 0
    missing_pae = 0

    with out_tsv.open("w") as out:
        out.write("\t".join([
            "uniprot","fragA_id","fragA_start","fragA_end",
            "fragB_id","fragB_start","fragB_end",
            "win","avg_pae","median_pae","n_pairs","frac_gt15",
            "classification","note"
        ]) + "\n")

        for acc, lst in segs.items():
            proteins += 1
            if len(lst) < 2:
                # single fragment → no interface to check
                continue
            pae_path = staged_pae_path(args.set, acc)
            pae = load_pae_matrix(pae_path)
            if pae is None:
                missing_pae += 1
                for i in range(len(lst)-1):
                    a = lst[i]; b = lst[i+1]
                    out.write("\t".join([
                        acc, a[2], str(a[0]), str(a[1]),
                        b[2], str(b[0]), str(b[1]),
                        str(args.win), "", "", "0", "",
                        "no_pae", "missing_pae_json"
                    ]) + "\n")
                    checks += 1
                continue

            N = pae.shape[0]
            # For safety, skip if indices exceed PAE size (rare); note and continue
            for i in range(len(lst)-1):
                a = lst[i]; b = lst[i+1]
                if a[1] < 1 or b[0] > N:
                    out.write("\t".join([
                        acc, a[2], str(a[0]), str(a[1]),
                        b[2], str(b[0]), str(b[1]),
                        str(args.win), "", "", "0", "",
                        "no_pae", "indices_out_of_bounds"
                    ]) + "\n")
                    checks += 1
                    continue
                avg, med, n_pairs, frac = interface_score(pae, a_end=a[1], b_start=b[0], win=args.win)
                cls = classify(avg, low=args.low, high=args.high)
                if cls == "support_split": c_support += 1
                elif cls == "suggest_merge": c_merge += 1
                elif cls == "ambiguous": c_amb += 1
                checks += 1
                out.write("\t".join([
                    acc, a[2], str(a[0]), str(a[1]),
                    b[2], str(b[0]), str(b[1]),
                    str(args.win),
                    f"{avg:.2f}" if np.isfinite(avg) else "",
                    f"{med:.2f}" if np.isfinite(med) else "",
                    str(n_pairs), f"{frac:.3f}" if np.isfinite(frac) else "",
                    cls, ""
                ]) + "\n")

    out_js.write_text(json.dumps({
        "set": args.set,
        "proteins_considered": proteins,
        "interfaces_checked": checks,
        "class_counts": {
            "support_split": c_support,
            "suggest_merge": c_merge,
            "ambiguous": c_amb,
            "no_pae": missing_pae  # count proteins with missing pae; individual rows also labeled
        },
        "window": args.win,
        "thresholds": {"low_merge": args.low, "high_split": args.high},
        "output_tsv": str(out_tsv)
    }, indent=2))
    print(f"Wrote {out_tsv}")
    print(f"Wrote {out_js}")

if __name__ == "__main__":
    main()
