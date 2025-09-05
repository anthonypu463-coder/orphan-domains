#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
from typing import List, Tuple, Dict, Optional
import numpy as np

SEG_ROOT   = Path("data/outputs/segments")                   # from 7B
FEAT_ROOT  = Path("data/outputs/structures/features")        # NPZ per protein

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["orphanZ70_reps","strict26k_reps","orphanZ70_all"], default="orphanZ70_reps")
    ap.add_argument("--min-fragment", type=int, default=50, help="minimum fragment length to keep")
    ap.add_argument("--limit", type=int, default=0, help="limit proteins (0=all)")
    return ap.parse_args()

def load_plddt(acc: str, set_name: str) -> Optional[np.ndarray]:
    npz = FEAT_ROOT / set_name / "by_acc" / f"{acc}.npz"
    if not npz.exists(): return None
    try:
        arr = np.load(npz, allow_pickle=True)
        return arr["plddt"].astype(float)
    except Exception:
        return None

def recompute_stats(plddt: Optional[np.ndarray], s: int, e: int) -> Dict[str, float]:
    if plddt is None or plddt.size == 0:
        return {"plddt_min": 0.0, "plddt_mean": 0.0, "plddt_max": 0.0, "n_low_lt50": 0}
    frag = plddt[s-1:e]
    return {
        "plddt_min": float(frag.min()),
        "plddt_mean": float(frag.mean()),
        "plddt_max": float(frag.max()),
        "n_low_lt50": int((frag < 50.0).sum()),
    }

def choose_merge_dir(left_len: int, right_len: int, left_mean: float, right_mean: float) -> str:
    if left_len > right_len: return "left"
    if right_len > left_len: return "right"
    # tie: higher pLDDT mean
    if left_mean > right_mean: return "left"
    if right_mean > left_mean: return "right"
    return "left"  # final tie-break

def main() -> None:
    args = parse_args()
    set_dir = SEG_ROOT / args.set
    seg_tsv = set_dir / "segments_by_plddt.tsv"
    if not seg_tsv.exists():
        raise SystemExit(f"Missing segments: {seg_tsv}. Run 7B first.")

    # Load segments grouped by protein
    by_acc: Dict[str, List[dict]] = {}
    with seg_tsv.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            acc = r["uniprot"]
            r["start"] = int(r["start"]); r["end"] = int(r["end"]); r["length"] = int(r["length"])
            r["plddt_mean"] = float(r["plddt_mean"])
            by_acc.setdefault(acc, []).append(r)
    # ensure ordered
    for acc in by_acc:
        by_acc[acc].sort(key=lambda x: (x["start"], x["end"]))

    if args.limit and args.limit < len(by_acc):
        # keep deterministic subset
        keep = dict(sorted(by_acc.items())[:args.limit])
        by_acc = keep

    out_dir = set_dir
    out_filtered = out_dir / "segments_filtered.tsv"
    out_merges   = out_dir / "segments_merges.tsv"
    out_summary  = out_dir / "segments_filtered_summary.json"

    kept = 0
    dropped = 0
    merges = 0
    prots = 0

    with out_filtered.open("w") as f_out, out_merges.open("w") as m_out:
        f_out.write("\t".join([
            "uniprot","fragment_id_new","start","end","length",
            "plddt_min","plddt_mean","plddt_max","n_low_lt50","source_frags","action"
        ]) + "\n")
        m_out.write("\t".join([
            "uniprot","action","target_frag","neighbor","merge_dir","result_start","result_end","result_length"
        ]) + "\n")

        for acc, frags in by_acc.items():
            prots += 1
            plddt = load_plddt(acc, args.set)
            # Build working list of (start,end,src_ids,mean,len)
            work: List[dict] = []
            for i, r in enumerate(frags, start=1):
                work.append({
                    "start": r["start"], "end": r["end"], "src": [r["fragment_id"]],
                    "mean": r["plddt_mean"], "length": r["length"]
                })

            # Iteratively resolve short fragments
            changed = True
            while changed:
                changed = False
                # Drop terminal shorts
                if work and work[0]["length"] < args.min_fragment:
                    # drop first
                    dropped += 1
                    work.pop(0)
                    changed = True
                    continue
                if work and work[-1]["length"] < args.min_fragment:
                    dropped += 1
                    work.pop()
                    changed = True
                    continue
                # Merge interior shorts
                for i in range(1, len(work)-1):
                    frag = work[i]
                    if frag["length"] >= args.min_fragment:
                        continue
                    left, right = work[i-1], work[i+1]
                    dir_ = choose_merge_dir(left["length"], right["length"], left["mean"], right["mean"])
                    if dir_ == "left":
                        # merge into left: extend left.end
                        new = {
                            "start": left["start"],
                            "end": frag["end"],
                            "src": left["src"] + frag["src"],
                            "mean": (left["mean"]*left["length"] + frag["mean"]*frag["length"]) / (left["length"]+frag["length"]) if (left["length"]+frag["length"])>0 else 0.0,
                            "length": left["length"] + frag["length"]
                        }
                        work[i-1] = new
                        # remove current
                        work.pop(i)
                        merges += 1
                        m_out.write("\t".join([acc,"merge","+".join(frag['src']), "+".join(left['src']), "left",
                                               str(new["start"]), str(new["end"]), str(new["length"])]) + "\n")
                        changed = True
                        break
                    else:
                        # merge into right: extend right.start
                        new = {
                            "start": frag["start"],
                            "end": right["end"],
                            "src": frag["src"] + right["src"],
                            "mean": (right["mean"]*right["length"] + frag["mean"]*frag["length"]) / (right["length"]+frag["length"]) if (right["length"]+frag["length"])>0 else 0.0,
                            "length": right["length"] + frag["length"]
                        }
                        work[i+1] = new
                        work.pop(i)
                        merges += 1
                        m_out.write("\t".join([acc,"merge","+".join(frag['src']), "+".join(right['src']), "right",
                                               str(new["start"]), str(new["end"]), str(new["length"])]) + "\n")
                        changed = True
                        break

            # Emit final kept fragments with recomputed stats
            for j, frag in enumerate(work, start=1):
                stats = recompute_stats(plddt, frag["start"], frag["end"])
                f_out.write("\t".join([
                    acc, f"F{j:03d}", str(frag["start"]), str(frag["end"]), str(frag["length"]),
                    f"{stats['plddt_min']:.2f}", f"{stats['plddt_mean']:.2f}", f"{stats['plddt_max']:.2f}",
                    str(stats["n_low_lt50"]), "+".join(frag["src"]), "keep"
                ]) + "\n")
                kept += 1

    out_summary.write_text(json.dumps({
        "set": args.set,
        "min_fragment": args.min_fragment,
        "proteins": prots,
        "fragments_kept_final": kept,
        "fragments_merged": merges,
        "fragments_dropped_terminal": dropped
    }, indent=2))
    print(f"Wrote {out_filtered}")
    print(f"Wrote {out_merges}")
    print(f"Wrote {out_summary}")

if __name__ == "__main__":
    main()
