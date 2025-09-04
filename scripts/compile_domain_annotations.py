#!/usr/bin/env python3
from __future__ import annotations
import csv, json
from pathlib import Path
from typing import Dict, List, Tuple, Optional

MANIFEST = Path("data/afdb/index/reference_manifest.tsv")               # protein lengths
COVER    = Path("data/outputs/hmmscan/pfam_coverage_nonoverlap.tsv")    # 2D output
DOMS     = Path("data/outputs/hmmscan/pfam_domains_nonoverlap.tsv")     # 2C output
OUT_TSV  = Path("data/outputs/hmmscan/domain_annotations.tsv")
OUT_JSON = Path("data/outputs/hmmscan/domain_annotations.jsonl")

def to_int(s: str) -> int:
    try: return int(s)
    except Exception: return 0

def to_float(s: str) -> Optional[float]:
    try: return float(s)
    except Exception: return None

def load_lengths(p: Path) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    with p.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            acc = row.get("uniprot", "").strip()
            if acc:
                lengths[acc] = to_int(row.get("residues","0"))
    if not lengths:
        raise SystemExit("No lengths found in manifest")
    return lengths

def load_coverage(p: Path) -> Dict[str, Tuple[int, int, float]]:
    cov: Dict[str, Tuple[int,int,float]] = {}
    with p.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            acc = row["uniprot"]
            cov[acc] = (
                to_int(row.get("pfam_domains","0")),
                to_int(row.get("covered_residues","0")),
                float(row.get("coverage_frac","0") or 0.0),
            )
    return cov

def load_domains(p: Path) -> Dict[str, List[dict]]:
    by: Dict[str, List[dict]] = {}
    with p.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            acc = row["uniprot"]
            s, e = to_int(row["env_from"]), to_int(row["env_to"])
            if s <= 0 or e <= 0 or s > e:
                continue
            by.setdefault(acc, []).append({
                "pfam_id": row["pfam_id"],
                "pfam_acc": row["pfam_acc"],
                "start": s,
                "end": e,
                "score": to_float(row.get("dom_score","")),
                "full_score": to_float(row.get("full_score","")),
                "description": row.get("description",""),
            })
    # sort within protein
    for lst in by.values():
        lst.sort(key=lambda d: (d["start"], d["end"]))
    return by

def format_domains_compact(domains: List[dict]) -> str:
    # Example token: EVA_Class_A(PF19429.4):38-110|60.40
    toks = []
    for d in domains:
        sc = f"{d['score']:.2f}" if d["score"] is not None else ""
        toks.append(f"{d['pfam_id']}({d['pfam_acc']}):{d['start']}-{d['end']}|{sc}")
    return "; ".join(toks)

def main() -> None:
    lengths = load_lengths(MANIFEST)
    coverage = load_coverage(COVER)
    domains = load_domains(DOMS)

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    tsv = OUT_TSV.open("w")
    tsv.write("\t".join(["uniprot","residues","coverage_frac","n_domains","domains"]) + "\n")

    with OUT_JSON.open("w") as js:
        for acc, L in lengths.items():
            n_dom, covered, frac = coverage.get(acc, (0, 0, 0.0))
            dom_list = domains.get(acc, [])
            # TSV row
            tsv.write(f"{acc}\t{L}\t{frac:.4f}\t{len(dom_list)}\t{format_domains_compact(dom_list)}\n")
            # JSONL row
            js.write(json.dumps({
                "uniprot": acc,
                "residues": L,
                "coverage_frac": round(frac, 6),
                "covered_residues": covered,
                "n_domains": len(dom_list),
                "domains": dom_list,
            }, ensure_ascii=False) + "\n")

    tsv.close()
    print(f"Wrote {OUT_TSV}")
    print(f"Wrote {OUT_JSON}")

if __name__ == "__main__":
    main()
