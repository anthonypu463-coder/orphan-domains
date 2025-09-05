#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np

SEG_ROOT   = Path("data/outputs/segments")
STRUCT_DIR = Path("data/structures")
FEAT_ROOT  = Path("data/outputs/structures/features")

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["orphanZ70_reps","strict26k_reps","orphanZ70_all"], default="orphanZ70_reps")
    ap.add_argument("--win", type=int, default=10, help="interface window for PAE")
    ap.add_argument("--low", type=float, default=5.0, help="PAE < low → merge")
    ap.add_argument("--high", type=float, default=15.0, help="PAE > high → supports split")
    ap.add_argument("--internal_high", type=float, default=18.0, help="PAE > internal_high → consider internal split")
    ap.add_argument("--min_fragment", type=int, default=50, help="fragments must be ≥ this after ops")
    ap.add_argument("--min_split_len", type=int, default=120, help="only consider internal cuts in fragments ≥ this")
    ap.add_argument("--min_sep", type=int, default=30, help="min distance between cut and fragment ends")
    ap.add_argument("--limit", type=int, default=0, help="limit proteins (0=all)")
    return ap.parse_args()

def load_segments(set_name: str) -> Dict[str, List[Dict[str,int]]]:
    seg_f = SEG_ROOT / set_name / "segments_final.tsv"
    if not seg_f.exists():
        raise SystemExit(f"Missing {seg_f}; run 7D first.")
    by: Dict[str, List[Dict[str,int]]] = {}
    with seg_f.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            d = {
                "uniprot": r["uniprot"],
                "fid": r["fragment_id"],
                "start": int(r["start"]),
                "end": int(r["end"]),
                "length": int(r["length"]),
            }
            by.setdefault(d["uniprot"], []).append(d)
    for acc in by:
        by[acc].sort(key=lambda x: (x["start"], x["end"]))
    return by

def load_pae_checks(set_name: str) -> Dict[Tuple[str,str,str], str]:
    """Map (acc, fidA, fidB) -> classification."""
    tsv = SEG_ROOT / set_name / "pae_interface_checks.tsv"
    if not tsv.exists():
        return {}
    m: Dict[Tuple[str,str,str], str] = {}
    with tsv.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            m[(r["uniprot"], r["fragA_id"], r["fragB_id"])] = r["classification"]
    return m

def staged_pae_path(set_name: str, acc: str) -> Path:
    return STRUCT_DIR / set_name / "by_acc" / acc / "pae.json"

def load_pae_matrix(p: Path) -> Optional[np.ndarray]:
    if not p.exists(): return None
    try:
        data = json.loads(p.read_text())
    except Exception:
        return None
    mat = None
    for k in ("predicted_aligned_error","pae","pae_matrix"):
        v = data.get(k) if isinstance(data, dict) else None
        if isinstance(v, list) and v and isinstance(v[0], list):
            mat = v; break
    if mat is None:
        return None
    try:
        arr = np.array(mat, dtype=float)
        if arr.ndim==2 and arr.shape[0]==arr.shape[1]:
            return arr
    except Exception:
        return None
    return None

def clamp(i:int, lo:int, hi:int)->int:
    return max(lo, min(hi, i))

def interface_avg(pae: np.ndarray, a_end: int, b_start: int, win: int) -> float:
    n = pae.shape[0]
    a_lo = clamp(a_end - win + 1, 1, n); a_hi = clamp(a_end, 1, n)
    b_lo = clamp(b_start, 1, n);         b_hi = clamp(b_start + win - 1, 1, n)
    A = np.arange(a_lo-1, a_hi, dtype=int); B = np.arange(b_lo-1, b_hi, dtype=int)
    if A.size==0 or B.size==0: return np.nan
    return float(np.nanmean(pae[np.ix_(A,B)]))

def internal_best_cut(pae: np.ndarray, s:int, e:int, win:int, internal_high:float, min_sep:int) -> Optional[int]:
    """Return a single strong internal cut position if any (conservative)."""
    n = pae.shape[0]
    lo = s + max(win, min_sep); hi = e - max(win, min_sep)
    if hi - lo + 1 <= 0: return None
    best_pos, best_avg = None, -1.0
    for pos in range(lo, hi+1):
        a_end = pos - 1; b_start = pos
        avg = interface_avg(pae, a_end, b_start, win)
        if np.isnan(avg): continue
        if avg > best_avg:
            best_avg, best_pos = avg, pos
    if best_pos is not None and best_avg >= internal_high:
        return best_pos
    return None

def load_plddt(acc: str, set_name: str) -> Optional[np.ndarray]:
    npz = FEAT_ROOT / set_name / "by_acc" / f"{acc}.npz"
    if not npz.exists(): return None
    try:
        arr = np.load(npz, allow_pickle=True)
        return arr["plddt"].astype(float)
    except Exception:
        return None

def frag_stats(plddt: Optional[np.ndarray], s:int, e:int) -> Tuple[float,float,float,int]:
    if plddt is None or plddt.size==0:
        return 0.0,0.0,0.0,0
    seg = plddt[s-1:e]
    return float(seg.min()), float(seg.mean()), float(seg.max()), int((seg<50.0).sum())

def main() -> None:
    args = parse_args()
    segs = load_segments(args.set)
    checks = load_pae_checks(args.set)

    # limit proteins deterministically
    accs = sorted(segs.keys())
    if args.limit and args.limit < len(accs):
        accs = accs[:args.limit]

    # outputs
    out_dir = SEG_ROOT / args.set
    out_ref = out_dir / "segments_refined.tsv"
    out_log = out_dir / "pae_actions.tsv"
    out_sum = out_dir / "segments_refine_summary.json"

    merges = splits = 0
    proteins = 0

    with out_ref.open("w") as fout, out_log.open("w") as flog:
        fout.write("\t".join(["uniprot","fragment_id","start","end","length","plddt_min","plddt_mean","plddt_max","n_low_lt50","source"])+"\n")
        flog.write("\t".join(["uniprot","action","details"])+"\n")

        for acc in accs:
            lst = segs[acc]  # list of dict with keys: fid, start, end, length
            if len(lst)==0:
                continue
            proteins += 1

            # 1) MERGE adjacent pairs if suggest_merge
            merged: List[Dict[str,int]] = []
            i = 0
            while i < len(lst):
                cur = lst[i]
                if i < len(lst)-1:
                    nxt = lst[i+1]
                    cls = checks.get((acc, cur["fid"], nxt["fid"]))
                    if cls == "suggest_merge":
                        # merge cur+nxt
                        new = {"uniprot": acc, "fid": f"{cur['fid']}+{nxt['fid']}",
                               "start": cur["start"], "end": nxt["end"]}
                        new["length"] = new["end"] - new["start"] + 1
                        merged.append(new)
                        flog.write("\t".join([acc,"merge",f"{cur['fid']}+{nxt['fid']} (PAE<low)"])+"\n")
                        merges += 1
                        i += 2
                        continue
                merged.append(cur)
                i += 1

            # 2) OPTIONAL INTERNAL SPLIT — conservative (at most one per fragment)
            pae = load_pae_matrix(staged_pae_path(args.set, acc))
            refined: List[Dict[str,int]] = []
            for frag in merged:
                s,e = frag["start"], frag["end"]
                L = e - s + 1
                if pae is None or L < args.min_split_len:
                    refined.append(frag)
                    continue
                cut = internal_best_cut(pae, s, e, args.win, args.internal_high, args.min_sep)
                if cut is None:
                    refined.append(frag)
                    continue
                # enforce ≥ min_fragment on both sides
                left_len  = (cut-1) - s + 1
                right_len = e - cut + 1
                if left_len < args.min_fragment or right_len < args.min_fragment:
                    refined.append(frag)  # don't split if would create small piece
                    continue
                # commit split
                left = {"uniprot": acc, "fid": f"{frag['fid']}_L", "start": s, "end": cut-1}
                right= {"uniprot": acc, "fid": f"{frag['fid']}_R", "start": cut, "end": e}
                left["length"]  = left["end"]  - left["start"]  + 1
                right["length"] = right["end"] - right["start"] + 1
                refined.extend([left,right])
                flog.write("\t".join([acc,"split",f"{frag['fid']} at {cut} (PAE>internal_high)"])+"\n")
                splits += 1

            # 3) Recompute stats, renumber F001.., write
            plddt = load_plddt(acc, args.set)
            refined.sort(key=lambda x: (x["start"], x["end"]))
            for j, fr in enumerate(refined, start=1):
                pmin,pmean,pmax,nlow = frag_stats(plddt, fr["start"], fr["end"])
                fout.write("\t".join([
                    acc, f"F{j:03d}", str(fr["start"]), str(fr["end"]), str(fr["length"]),
                    f"{pmin:.2f}", f"{pmean:.2f}", f"{pmax:.2f}", str(nlow),
                    "pae_merge" if "+" in fr["fid"] else ("pae_split" if fr["fid"].endswith(("_L","_R")) else "kept")
                ]) + "\n")

    out_sum.write_text(json.dumps({
        "set": args.set,
        "proteins": proteins,
        "merges": merges,
        "splits": splits,
        "win": args.win,
        "thresholds": {"low_merge": args.low, "high_split": args.high, "internal_high": args.internal_high},
        "min_fragment": args.min_fragment,
        "min_split_len": args.min_split_len,
        "min_sep": args.min_sep,
        "refined_table": str(out_ref),
        "actions_log": str(out_log)
    }, indent=2))
    print(f"Wrote {out_ref}")
    print(f"Wrote {out_log}")
    print(f"Wrote {out_sum}")

if __name__ == "__main__":
    main()
