#!/usr/bin/env python3
from __future__ import annotations
import csv, json, statistics
from pathlib import Path
from typing import Dict, List

# Inputs
META = Path("data/outputs/orphans/orphan_metadata.tsv")             # from 4B
MANI = Path("data/afdb/index/reference_manifest.tsv")               # has plddt_mean and residues

# Outputs
OUT_DIR = Path("data/outputs/orphans")
QC_TSV  = OUT_DIR / "orphan_candidates_qc.tsv"
VIOL_TSV= OUT_DIR / "coverage_violations.tsv"
STRUCT_TXT = OUT_DIR / "orphan_structural_subset.txt"
SUM_JSON= OUT_DIR / "qc_summary.json"

CUTOFF = 0.20
PLDDT_MIN = 50.0  # flag low-confidence proteins

def to_float(x: str) -> float:
    try: return float(x)
    except Exception: return 0.0

def to_int(x: str) -> int:
    try: return int(x)
    except Exception: return 0

def load_manifest(path: Path) -> Dict[str, Dict[str, float]]:
    m: Dict[str, Dict[str, float]] = {}
    with path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            acc = r["uniprot"].strip()
            if not acc: continue
            m[acc] = {
                "residues": to_int(r.get("residues","0")),
                "plddt_mean": to_float(r.get("plddt_mean","0")),
            }
    return m

def main() -> None:
    if not META.exists(): raise SystemExit(f"Missing {META}")
    if not MANI.exists(): raise SystemExit(f"Missing {MANI}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    mani = load_manifest(MANI)

    # Load candidate metadata
    rows: List[Dict[str, str]] = []
    with META.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            acc = r["uniprot"].strip()
            if not acc: continue
            rows.append(r)

    # QC + join
    qc_rows: List[Dict[str, str]] = []
    viol: List[Dict[str, str]] = []
    plddts: List[float] = []
    lengths: List[int] = []
    low_conf: List[str] = []
    structural_subset: List[str] = []

    for r in rows:
        acc = r["uniprot"]
        cov = to_float(r.get("coverage_frac","0"))
        L_meta = to_int(r.get("residues","0"))
        m = mani.get(acc, {"residues": L_meta, "plddt_mean": 0.0})
        pmean = m["plddt_mean"]
        L = m["residues"] or L_meta

        # Coverage check
        cov_ok = cov < CUTOFF
        if not cov_ok:
            viol.append({"uniprot": acc, "coverage_frac": f"{cov:.4f}"})

        # Confidence flags
        low = pmean < PLDDT_MIN
        if low: low_conf.append(acc)
        else: structural_subset.append(acc)

        plddts.append(pmean)
        lengths.append(L)

        qc_rows.append({
            "uniprot": acc,
            "residues": str(L),
            "coverage_frac": f"{cov:.4f}",
            "n_domains": r.get("n_domains",""),
            "plddt_mean": f"{pmean:.2f}",
            "low_confidence": "1" if low else "0"
        })

    # Write QC table
    with QC_TSV.open("w") as fh:
        fh.write("\t".join(["uniprot","residues","coverage_frac","n_domains","plddt_mean","low_confidence"]) + "\n")
        for q in qc_rows:
            fh.write("\t".join([q["uniprot"], q["residues"], q["coverage_frac"], q["n_domains"], q["plddt_mean"], q["low_confidence"]]) + "\n")

    # Coverage violations (should be empty)
    with VIOL_TSV.open("w") as fh:
        fh.write("\t".join(["uniprot","coverage_frac"]) + "\n")
        for v in viol:
            fh.write(f"{v['uniprot']}\t{v['coverage_frac']}\n")

    # Structural subset (pLDDT >= 50)
    with STRUCT_TXT.open("w") as fh:
        for acc in structural_subset:
            fh.write(acc + "\n")

    # Summaries
    def safe_stats(vals: List[float]) -> Dict[str, float]:
        if not vals: return {"mean": 0.0, "median": 0.0, "min": 0.0, "max": 0.0, "p10": 0.0, "p90": 0.0}
        vals_sorted = sorted(vals)
        n = len(vals_sorted)
        return {
            "mean": round(statistics.fmean(vals_sorted), 3),
            "median": round(statistics.median(vals_sorted), 3),
            "min": round(vals_sorted[0], 3),
            "max": round(vals_sorted[-1], 3),
            "p10": round(vals_sorted[int(0.10*(n-1))], 3),
            "p90": round(vals_sorted[int(0.90*(n-1))], 3),
        }

    SUM_JSON.write_text(json.dumps({
        "cutoff": CUTOFF,
        "candidates": len(rows),
        "coverage_violations": len(viol),
        "low_confidence_count": len(low_conf),
        "structural_subset_count": len(structural_subset),
        "plddt_mean_stats": safe_stats(plddts),
        "length_stats": safe_stats([float(x) for x in lengths]),
    }, indent=2))

    print(f"Wrote {QC_TSV}")
    print(f"Wrote {VIOL_TSV} (expected near-empty)")
    print(f"Wrote {STRUCT_TXT}")
    print(f"Wrote {SUM_JSON}")

if __name__ == "__main__":
    main()
