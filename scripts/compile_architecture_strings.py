#!/usr/bin/env python3
from __future__ import annotations
import csv, json
from pathlib import Path
from typing import Dict, List, Tuple

MANIFEST = Path("data/afdb/index/reference_manifest.tsv")          # uniprot, residues
DOM_FINAL= Path("data/outputs/architecture/pfam_domains_final.tsv")# uniprot, pfam_acc, start, end, ...
OUT_TSV  = Path("data/outputs/architecture/pfam_architecture_strings.tsv")
OUT_JSON = Path("data/outputs/architecture/pfam_architecture_strings.jsonl")

def to_int(s: str) -> int:
    try: return int(s)
    except Exception: return 0

def norm_acc(acc: str) -> str:
    acc = (acc or "").strip()
    if not acc: return ""
    i = acc.find(".")
    return acc[:i] if i > 0 else acc

def load_lengths(path: Path) -> Dict[str, int]:
    m: Dict[str, int] = {}
    with path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            acc = r.get("uniprot","").strip()
            if acc:
                m[acc] = to_int(r.get("residues","0"))
    if not m:
        raise SystemExit(f"No proteins in {path}")
    return m

def load_domains(path: Path) -> Dict[str, List[Tuple[int,int,str]]]:
    by: Dict[str, List[Tuple[int,int,str]]] = {}
    if not path.exists():
        return by
    with path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            acc = r["uniprot"]
            s, e = to_int(r["start"]), to_int(r["end"])
            if s <= 0 or e < s: 
                continue
            token = norm_acc(r.get("pfam_acc","")) or r.get("pfam_id","")
            by.setdefault(acc, []).append((s, e, token))
    # sort by coordinates
    for lst in by.values():
        lst.sort(key=lambda x: (x[0], x[1]))
    return by

def merge_overlaps(iv: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    if not iv: return []
    iv = sorted(iv)
    out = [iv[0]]
    for s,e in iv[1:]:
        ps, pe = out[-1]
        if s <= pe:
            out[-1] = (ps, max(pe, e))
        else:
            out.append((s,e))
    return out

def gaps_from_intervals(L: int, iv: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    gaps: List[Tuple[int,int]] = []
    if L <= 0:
        return gaps
    pos = 1
    for s, e in iv:
        if s > pos:
            gaps.append((pos, s-1))
        pos = max(pos, e+1)
    if pos <= L:
        gaps.append((pos, L))
    return gaps

def fmt_arch(tokens: List[str]) -> str:
    return "–".join(tokens) if tokens else "None"

def fmt_spans(spans: List[Tuple[int,int]]) -> str:
    return ";".join(f"{s}-{e}" for s,e in spans) if spans else "None"

def fmt_lengths(spans: List[Tuple[int,int]]) -> str:
    return ";".join(str(e - s + 1) for s,e in spans) if spans else "0"

def main() -> None:
    lengths = load_lengths(MANIFEST)
    domains = load_domains(DOM_FINAL)

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    tsv = OUT_TSV.open("w")
    tsv.write("\t".join([
        "uniprot","residues","coverage_frac","n_domains",
        "domain_arch","n_unannotated","unannotated_spans","unannotated_lengths"
    ]) + "\n")

    with OUT_JSON.open("w") as js:
        for acc, L in lengths.items():
            lst = domains.get(acc, [])
            iv = [(s,e) for s,e,_ in lst]
            iv_merged = merge_overlaps(iv)  # final should be non-overlap; safe to merge
            covered = sum(e - s + 1 for s,e in iv_merged)
            frac = (covered / L) if L > 0 else 0.0
            tokens = [tok for _,_,tok in lst]
            arch = fmt_arch(tokens)
            gaps = gaps_from_intervals(L, iv_merged)
            row = [
                acc, str(L), f"{frac:.4f}", str(len(lst)),
                arch, str(len(gaps)), fmt_spans(gaps), fmt_lengths(gaps)
            ]
            tsv.write("\t".join(row) + "\n")
            js.write(json.dumps({
                "uniprot": acc,
                "residues": L,
                "coverage_frac": round(frac, 6),
                "n_domains": len(lst),
                "domain_arch": tokens,
                "unannotated_spans": gaps,
                "unannotated_lengths": [e - s + 1 for s,e in gaps],
            }) + "\n")

    tsv.close()
    print(f"Wrote {OUT_TSV}")
    print(f"Wrote {OUT_JSON}")

if __name__ == "__main__":
    main()
