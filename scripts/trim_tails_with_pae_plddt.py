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
    ap.add_argument("--pae_thr", type=float, default=15.0, help="tail-vs-core PAE avg > this → weakly attached")
    ap.add_argument("--plddt_thr", type=float, default=50.0, help="tail mean pLDDT < this → disordered")
    ap.add_argument("--tail_max", type=int, default=25, help="max residues to trim from each end")
    ap.add_argument("--core_skip", type=int, default=10, help="skip these residues around the cut when defining core")
    ap.add_argument("--min_fragment", type=int, default=50, help="final fragments must be ≥ this")
    ap.add_argument("--limit", type=int, default=0, help="limit proteins (0=all)")
    return ap.parse_args()

def load_refined(set_name: str) -> Dict[str, List[Dict[str,int]]]:
    ref = SEG_ROOT / set_name / "segments_refined.tsv"
    if not ref.exists():
        raise SystemExit(f"Missing {ref}; run 8B first.")
    by: Dict[str, List[Dict[str,int]]] = {}
    with ref.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            d = {"uniprot": r["uniprot"], "fid": r["fragment_id"],
                 "start": int(r["start"]), "end": int(r["end"]),
                 "length": int(r["length"])}
            by.setdefault(d["uniprot"], []).append(d)
    for acc in by:
        by[acc].sort(key=lambda x: (x["start"], x["end"]))
    return by

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
    if mat is None: return None
    try:
        arr = np.array(mat, dtype=float)
        if arr.ndim==2 and arr.shape[0]==arr.shape[1]: return arr
    except Exception:
        return None
    return None

def load_plddt(acc: str, set_name: str) -> Optional[np.ndarray]:
    npz = FEAT_ROOT / set_name / "by_acc" / f"{acc}.npz"
    if not npz.exists(): return None
    try:
        arr = np.load(npz, allow_pickle=True)
        return arr["plddt"].astype(float)
    except Exception:
        return None

def clamp(i:int, lo:int, hi:int)->int:
    return max(lo, min(hi, i))

def pae_avg_tail_vs_core(pae: np.ndarray, tail_idx: np.ndarray, core_idx: np.ndarray) -> float:
    if tail_idx.size==0 or core_idx.size==0: return np.nan
    sub = pae[np.ix_(tail_idx, core_idx)]
    return float(np.nanmean(sub)) if sub.size else np.nan

def try_trim_end(pae: Optional[np.ndarray], plddt: Optional[np.ndarray], s:int, e:int,
                 end: str, tail_max:int, core_skip:int, pae_thr:float, plddt_thr:float, min_fragment:int) -> Tuple[int,int,int]:
    """
    Returns (new_s, new_e, trimmed_residues).
    'end' is 'N' or 'C'.
    """
    if pae is None or plddt is None or plddt.size==0: return s, e, 0
    n = pae.shape[0]
    new_s, new_e = s, e
    trimmed = 0
    # scan one-by-one, conservative
    for k in range(1, tail_max+1):
        if end == "N":
            t_start, t_end = s, s + k - 1
            rem_len = e - (s + k - 1) + 1
            if rem_len < min_fragment: break
            tail = np.arange(t_start-1, t_end, dtype=int)
            core_lo = clamp(t_end + 1 + core_skip, 1, n)
            core_hi = clamp(e - core_skip, 1, n)
            core = np.arange(core_lo-1, core_hi, dtype=int) if core_hi >= core_lo else np.array([], dtype=int)
            tail_pl = float(plddt[t_start-1 : t_end].mean()) if t_end >= t_start else 100.0
            pae_avg = pae_avg_tail_vs_core(pae, tail, core)
            if np.isfinite(pae_avg) and tail_pl < plddt_thr and pae_avg > pae_thr:
                new_s = s + k
                trimmed = k
            else:
                break
        else:
            t_start, t_end = e - k + 1, e
            rem_len = (e - k) - s + 1
            if rem_len < min_fragment: break
            tail = np.arange(t_start-1, t_end, dtype=int)
            core_lo = clamp(s + core_skip, 1, n)
            core_hi = clamp(t_start - 1 - core_skip, 1, n)
            core = np.arange(core_lo-1, core_hi, dtype=int) if core_hi >= core_lo else np.array([], dtype=int)
            tail_pl = float(plddt[t_start-1 : t_end].mean()) if t_end >= t_start else 100.0
            pae_avg = pae_avg_tail_vs_core(pae, tail, core)
            if np.isfinite(pae_avg) and tail_pl < plddt_thr and pae_avg > pae_thr:
                new_e = e - k
                trimmed = k
            else:
                break
    return new_s, new_e, trimmed

def frag_stats(plddt: Optional[np.ndarray], s:int, e:int) -> Tuple[float,float,float,int]:
    if plddt is None or plddt.size==0:
        return 0.0,0.0,0.0,0
    seg = plddt[s-1:e]
    return float(seg.min()), float(seg.mean()), float(seg.max()), int((seg<50.0).sum())

def main() -> None:
    args = parse_args()
    set_name = args.set
    segs = load_refined(set_name)
    accs = sorted(segs.keys())
    if args.limit and args.limit < len(accs):
        accs = accs[:args.limit]

    out_dir = SEG_ROOT / set_name
    out_trim = out_dir / "segments_trimmed.tsv"
    out_log  = out_dir / "pae_tail_trims.tsv"
    out_sum  = out_dir / "segments_trim_summary.json"

    n_prot = 0
    n_frags = 0
    n_trim_N = n_trim_C = 0

    with out_trim.open("w") as fout, out_log.open("w") as flog:
        fout.write("\t".join(["uniprot","fragment_id","start","end","length","plddt_min","plddt_mean","plddt_max","n_low_lt50","source"])+"\n")
        flog.write("\t".join(["uniprot","fragment_id","action","trimmed_res","pae_thr","plddt_thr"])+"\n")

        for acc in accs:
            plddt = load_plddt(acc, set_name)
            pae = load_pae_matrix(staged_pae_path(set_name, acc))
            if plddt is None or pae is None or pae.shape[0] != plddt.size:
                # Keep unchanged if we lack confidence data
                for fr in segs[acc]:
                    pmin,pmean,pmax,nlow = frag_stats(plddt, fr["start"], fr["end"])
                    fout.write("\t".join([acc, fr["fid"], str(fr["start"]), str(fr["end"]), str(fr["length"]),
                                          f"{pmin:.2f}", f"{pmean:.2f}", f"{pmax:.2f}", str(nlow), "kept_no_pae_or_plddt"])+"\n")
                    n_frags += 1
                n_prot += 1
                continue

            n_prot += 1
            for fr in segs[acc]:
                s, e = fr["start"], fr["end"]
                # N-tail
                new_s, new_e, tN = try_trim_end(pae, plddt, s, e, end="N",
                                                tail_max=args.tail_max, core_skip=args.core_skip,
                                                pae_thr=args.pae_thr, plddt_thr=args.plddt_thr,
                                                min_fragment=args.min_fragment)
                if tN > 0:
                    s = new_s; n_trim_N += 1
                    flog.write("\t".join([acc, fr["fid"], "trim_N", str(tN), str(args.pae_thr), str(args.plddt_thr)])+"\n")
                # C-tail
                new_s, new_e, tC = try_trim_end(pae, plddt, s, e, end="C",
                                                tail_max=args.tail_max, core_skip=args.core_skip,
                                                pae_thr=args.pae_thr, plddt_thr=args.plddt_thr,
                                                min_fragment=args.min_fragment)
                if tC > 0:
                    e = new_e; n_trim_C += 1
                    flog.write("\t".join([acc, fr["fid"], "trim_C", str(tC), str(args.pae_thr), str(args.plddt_thr)])+"\n")

                L = e - s + 1
                if L < args.min_fragment:
                    # Safety: never emit <50. If we trimmed too much, revert to original.
                    s, e = fr["start"], fr["end"]; L = e - s + 1
                pmin,pmean,pmax,nlow = frag_stats(plddt, s, e)
                fout.write("\t".join([acc, fr["fid"], str(s), str(e), str(L),
                                      f"{pmin:.2f}", f"{pmean:.2f}", f"{pmax:.2f}", str(nlow),
                                      "pae_tail_trim" if (tN>0 or tC>0) else "kept"])+"\n")
                n_frags += 1

    out_sum.write_text(json.dumps({
        "set": set_name,
        "proteins": n_prot,
        "fragments": n_frags,
        "tails_trimmed_N": n_trim_N,
        "tails_trimmed_C": n_trim_C,
        "params": {
            "pae_thr": args.pae_thr, "plddt_thr": args.plddt_thr,
            "tail_max": args.tail_max, "core_skip": args.core_skip,
            "min_fragment": args.min_fragment
        },
        "refined_input": str((SEG_ROOT / set_name / "segments_refined.tsv")),
        "trimmed_table": str(out_trim),
        "actions_log": str(out_log)
    }, indent=2))
    print(f"Wrote {out_trim}")
    print(f"Wrote {out_log}")
    print(f"Wrote {out_sum}")

if __name__ == "__main__":
    main()
