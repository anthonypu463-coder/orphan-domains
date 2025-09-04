#!/usr/bin/env python3
from __future__ import annotations
import csv, json
from pathlib import Path
from typing import Dict, List, Tuple, Optional

MANIFEST   = Path("data/afdb/index/reference_manifest.tsv")            # uniprot, residues
DOM_FINAL  = Path("data/outputs/architecture/pfam_domains_final.tsv")  # uniprot, pfam_id/acc, start, end, scores
OUT_TSV    = Path("data/outputs/architecture/pfam_architecture_enhanced.tsv")
OUT_JSONL  = Path("data/outputs/architecture/pfam_architecture_enhanced.jsonl")

def to_int(s: str) -> int:
    try: return int(s)
    except Exception: return 0

def to_float(s: str) -> Optional[float]:
    try: return float(s)
    except Exception: return None

def norm_acc(acc: str) -> str:
    acc = (acc or "").strip()
    if not acc: return ""
    i = acc.find(".")
    return acc[:i] if i > 0 else acc

def load_lengths(p: Path) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    with p.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            acc = r.get("uniprot","").strip()
            if acc:
                lengths[acc] = to_int(r.get("residues","0"))
    if not lengths:
        raise SystemExit(f"No proteins in {p}")
    return lengths

def load_domains(p: Path) -> Dict[str, List[dict]]:
    per: Dict[str, List[dict]] = {}
    with p.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            acc = r["uniprot"].strip()
            s, e = to_int(r["start"]), to_int(r["end"])
            if not acc or s <= 0 or e < s: 
                continue
            per.setdefault(acc, []).append({
                "pfam_id": r.get("pfam_id",""),
                "pfam_acc": r.get("pfam_acc",""),
                "pfam_acc_root": norm_acc(r.get("pfam_acc","")),
                "start": s,
                "end": e,
                "length": e - s + 1,
                "dom_score": to_float(r.get("dom_score","")),
                "full_score": to_float(r.get("full_score","")),
                "note": r.get("note",""),
            })
    for lst in per.values():
        lst.sort(key=lambda d: (d["start"], d["end"]))
    return per

def gaps_from_intervals(L: int, iv: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    gaps: List[Tuple[int,int]] = []
    if L <= 0: 
        return gaps
    pos = 1
    for s,e in iv:
        if s > pos:
            gaps.append((pos, s-1))
        pos = max(pos, e+1)
    if pos <= L:
        gaps.append((pos, L))
    return gaps

def main() -> None:
    lengths = load_lengths(MANIFEST)
    domains = load_domains(DOM_FINAL)

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TSV.open("w") as tsv, OUT_JSONL.open("w") as jsl:
        tsv.write("\t".join([
            "uniprot","residues","covered_residues","coverage_frac",
            "n_domains","domains_compact","gaps","gap_lengths"
        ]) + "\n")

        for acc, L in lengths.items():
            doms = domains.get(acc, [])
            iv   = [(d["start"], d["end"]) for d in doms]
            covered = sum(e - s + 1 for s, e in iv)
            frac = (covered / L) if L > 0 else 0.0
            gaps = gaps_from_intervals(L, iv)
            gap_lengths = [e - s + 1 for s, e in gaps]

            # Compact domain string e.g., PF19429:37-111|60.40;PF00001:...
            compact_parts: List[str] = []
            for d in doms:
                sc = f"{d['dom_score']:.2f}" if d["dom_score"] is not None else ""
                token = d["pfam_acc_root"] or d["pfam_id"]
                compact_parts.append(f"{token}:{d['start']}-{d['end']}|{sc}")
            compact = "; ".join(compact_parts) if compact_parts else ""

            # TSV row
            tsv.write("\t".join([
                acc, str(L), str(covered), f"{frac:.4f}",
                str(len(doms)), compact,
                ";".join(f"{s}-{e}" for s,e in gaps) if gaps else "",
                ";".join(str(x) for x in gap_lengths) if gap_lengths else ""
            ]) + "\n")

            # JSONL row (richer)
            jsl.write(json.dumps({
                "uniprot": acc,
                "residues": L,
                "covered_residues": covered,
                "coverage_frac": round(frac, 6),
                "n_domains": len(doms),
                "domains": doms,            # list of dicts with start/end/scores
                "gaps": gaps,               # list of [start, end]
                "gap_lengths": gap_lengths,
            }, ensure_ascii=False) + "\n")

    print(f"Wrote {OUT_TSV}")
    print(f"Wrote {OUT_JSONL}")

if __name__ == "__main__":
    main()
