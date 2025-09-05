#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
from typing import Dict, List, Optional
import numpy as np

SEG_ROOT   = Path("data/outputs/segments")
FEAT_ROOT  = Path("data/outputs/structures/features")
FEAT_INDEX = "features_index.tsv"

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["orphanZ70_reps","strict26k_reps","orphanZ70_all"], default="orphanZ70_reps")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--high_thr", type=float, default=70.0)
    ap.add_argument("--low_thr", type=float, default=50.0)
    return ap.parse_args()

def load_plddt(acc: str, set_name: str) -> Optional[np.ndarray]:
    npz = FEAT_ROOT / set_name / "by_acc" / f"{acc}.npz"
    if not npz.exists(): return None
    try:
        arr = np.load(npz, allow_pickle=True)
        return arr["plddt"].astype(float)
    except Exception:
        return None

def load_features_index(set_name: str) -> Dict[str, int]:
    idx = FEAT_ROOT / set_name / FEAT_INDEX
    m: Dict[str, int] = {}
    if not idx.exists(): return m
    with idx.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            try: m[r["uniprot"]] = int(r["n_res"])
            except Exception: pass
    return m

def main() -> None:
    args = parse_args()
    set_dir  = SEG_ROOT / args.set
    seg_filt = set_dir / "segments_filtered.tsv"
    if not seg_filt.exists():
        raise SystemExit(f"Missing {seg_filt}. Run 7C first.")

    # Group filtered segments by protein
    by_acc: Dict[str, List[dict]] = {}
    order: List[str] = []
    with seg_filt.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            acc = r["uniprot"]
            if acc not in by_acc:
                order.append(acc); by_acc[acc] = []
            # ensure numeric fields are parsed
            r["start"] = int(r["start"]); r["end"] = int(r["end"]); r["length"] = int(r["length"])
            r["plddt_min"] = float(r["plddt_min"]); r["plddt_mean"] = float(r["plddt_mean"]); r["plddt_max"] = float(r["plddt_max"])
            r["n_low_lt50"] = int(r["n_low_lt50"])
            by_acc[acc].append(r)
    for acc in by_acc:
        by_acc[acc].sort(key=lambda x: (x["start"], x["end"]))

    if args.limit and args.limit < len(order):
        keep = set(order[:args.limit])
        by_acc = {k: v for k, v in by_acc.items() if k in keep}
        order = [a for a in order if a in keep]

    out_final = set_dir / "segments_final.tsv"
    out_log   = set_dir / "segments_edge_cases.tsv"
    out_sum   = set_dir / "segments_edge_summary.json"

    kept_rows = forced_single = dropped_all_low = normal_kept = 0

    with out_final.open("w") as fout, out_log.open("w") as flog:
        fout.write("\t".join([
            "uniprot","fragment_id","start","end","length",
            "plddt_min","plddt_mean","plddt_max","n_low_lt50","source"
        ]) + "\n")
        flog.write("\t".join(["uniprot","status","min_plddt","mean_plddt","max_plddt","note"]) + "\n")

        for acc in order:
            plddt = load_plddt(acc, args.set)
            if plddt is None or plddt.size == 0:
                # emit filtered as-is
                for j, r in enumerate(by_acc[acc], start=1):
                    fout.write("\t".join([
                        acc, f"F{j:03d}",
                        str(r["start"]), str(r["end"]), str(r["length"]),
                        f"{r['plddt_min']:.2f}", f"{r['plddt_mean']:.2f}", f"{r['plddt_max']:.2f}",
                        str(r["n_low_lt50"]), "filtered"
                    ]) + "\n")
                    kept_rows += 1
                normal_kept += 1
                continue

            pmin, pmax, pmean = float(plddt.min()), float(plddt.max()), float(plddt.mean())
            N = int(plddt.size)

            if pmin >= args.high_thr:
                forced_single += 1
                n_low = int((plddt < 50.0).sum())
                fout.write("\t".join([
                    acc, "F001", "1", str(N), str(N),
                    f"{pmin:.2f}", f"{pmean:.2f}", f"{pmax:.2f}", str(n_low), "forced_single_high70"
                ]) + "\n")
                flog.write("\t".join([acc, "single_domain_highconf", f"{pmin:.2f}", f"{pmean:.2f}", f"{pmax:.2f}", "override_to_full_length"]) + "\n")
                kept_rows += 1
                continue

            if pmax < args.low_thr:
                dropped_all_low += 1
                flog.write("\t".join([acc, "no_stable_domain", f"{pmin:.2f}", f"{pmean:.2f}", f"{pmax:.2f}", "dropped"]) + "\n")
                continue

            normal_kept += 1
            for j, r in enumerate(by_acc[acc], start=1):
                fout.write("\t".join([
                    acc, f"F{j:03d}",
                    str(r["start"]), str(r["end"]), str(r["length"]),
                    f"{r['plddt_min']:.2f}", f"{r['plddt_mean']:.2f}", f"{r['plddt_max']:.2f}",
                    str(r["n_low_lt50"]), "filtered"
                ]) + "\n")
                kept_rows += 1

    out_sum.write_text(json.dumps({
        "set": args.set,
        "proteins_total": len(order),
        "final_segments": kept_rows,
        "forced_single_full": forced_single,
        "dropped_no_stable": dropped_all_low,
        "normal_kept": normal_kept,
        "high_thr": args.high_thr,
        "low_thr": args.low_thr
    }, indent=2))
    print(f"Wrote {out_final}")
    print(f"Wrote {out_log}")
    print(f"Wrote {out_sum}")

if __name__ == "__main__":
    main()
