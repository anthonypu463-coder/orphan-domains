#!/usr/bin/env python3
from __future__ import annotations
import csv, json
from pathlib import Path
from typing import Dict, List

ENHANCED   = Path("data/outputs/architecture/pfam_architecture_enhanced.tsv")  # from 3D
MANIFEST   = Path("data/afdb/index/reference_manifest.tsv")                    # plddt_mean
FASTA_ALL  = Path("data/afdb/index/swissprot_sequences.fasta")                 # from 2A
OUT_DIR    = Path("data/outputs/orphans")

META_OUT   = OUT_DIR / "orphan_zero_struct70_metadata.tsv"
FASTA_OUT  = OUT_DIR / "orphan_zero_struct70.fasta"
LIST_OUT   = OUT_DIR / "orphan_zero_struct70_accessions.txt"
MISS_OUT   = OUT_DIR / "orphan_zero_struct70_missing_sequences.txt"
SUM_OUT    = OUT_DIR / "orphan_zero_struct70_summary.json"

PLDDT_MIN  = 70.0
EPS        = 1e-9  # strict zero coverage

def to_int(s: str) -> int:
    try: return int(s)
    except Exception: return 0

def to_float(s: str) -> float:
    try: return float(s)
    except Exception: return 0.0

def load_manifest(path: Path) -> Dict[str, float]:
    m: Dict[str, float] = {}
    with path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            acc = r.get("uniprot","").strip()
            if not acc: continue
            m[acc] = to_float(r.get("plddt_mean","0"))
    return m

def load_enhanced(path: Path) -> List[dict]:
    rows: List[dict] = []
    with path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            rows.append(r)
    return rows

def collect_seqs(fasta: Path, wanted: set[str]) -> Dict[str,str]:
    seqs: Dict[str,str] = {}
    with fasta.open() as fh:
        acc=None; chunks=[]
        for ln in fh:
            if ln.startswith(">"):
                if acc in wanted and chunks:
                    seqs[acc] = "".join(chunks)
                acc = ln[1:].strip().split()[0]
                chunks=[]
            else:
                if acc in wanted:
                    chunks.append(ln.strip())
        if acc in wanted and chunks:
            seqs[acc] = "".join(chunks)
    return seqs

def main() -> None:
    assert ENHANCED.exists(), f"Missing {ENHANCED}"
    assert MANIFEST.exists(), f"Missing {MANIFEST}"
    assert FASTA_ALL.exists(), f"Missing {FASTA_ALL}"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    mani_plddt = load_manifest(MANIFEST)
    enh = load_enhanced(ENHANCED)

    # Filter: strictly zero coverage and mean pLDDT >= 70
    picked: List[dict] = []
    for r in enh:
        cov = to_float(r.get("coverage_frac","0"))
        if cov > EPS: 
            continue
        acc = r["uniprot"].strip()
        pmean = mani_plddt.get(acc, 0.0)
        if pmean >= PLDDT_MIN:
            picked.append({
                "uniprot": acc,
                "residues": r.get("residues",""),
                "covered_residues": r.get("covered_residues",""),
                "coverage_frac": f"{cov:.4f}",
                "n_domains": r.get("n_domains",""),
                "plddt_mean": f"{pmean:.2f}",
                "domains_compact": r.get("domains_compact",""),
                "gaps": r.get("gaps",""),
                "gap_lengths": r.get("gap_lengths",""),
            })

    # Order: longest first, then highest pLDDT
    picked.sort(key=lambda x: (-to_int(x["residues"]), -to_float(x["plddt_mean"])))

    # Write metadata
    with META_OUT.open("w") as fh:
        cols = ["uniprot","residues","covered_residues","coverage_frac",
                "n_domains","plddt_mean","domains_compact","gaps","gap_lengths"]
        fh.write("\t".join(cols) + "\n")
        for r in picked:
            fh.write("\t".join(r[c] for c in cols) + "\n")

    # Write accession list
    with LIST_OUT.open("w") as fh:
        for r in picked:
            fh.write(r["uniprot"] + "\n")

    # FASTA
    subset = set(r["uniprot"] for r in picked)
    seqs = collect_seqs(FASTA_ALL, subset)
    with FASTA_OUT.open("w") as fh:
        for r in picked:
            a = r["uniprot"]; s = seqs.get(a)
            if not s: 
                continue
            fh.write(f">{a}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i+60] + "\n")

    # Missing sequences report
    missing = [a for a in subset if a not in seqs]
    with MISS_OUT.open("w") as fh:
        for a in missing:
            fh.write(a + "\n")

    # Summary
    count = len(picked)
    max_cov = max((to_float(r["coverage_frac"]) for r in picked), default=0.0)
    min_plddt = min((to_float(r["plddt_mean"]) for r in picked), default=0.0)
    SUM_OUT.write_text(json.dumps({
        "criteria": {"coverage_frac": 0.0, "plddt_mean_min": PLDDT_MIN},
        "selected": count,
        "missing_sequences": len(missing),
        "max_coverage_frac_in_set": round(max_cov, 6),
        "min_plddt_in_set": round(min_plddt, 3)
    }, indent=2))

    print(f"Wrote {META_OUT}")
    print(f"Wrote {FASTA_OUT}")
    print(f"Wrote {LIST_OUT} and {MISS_OUT}")
    print(f"Wrote {SUM_OUT}")

if __name__ == "__main__":
    main()
