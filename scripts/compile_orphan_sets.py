#!/usr/bin/env python3
from __future__ import annotations
import csv, json
from pathlib import Path
from typing import Dict, List, Tuple

# Inputs (from 3D + 4A)
ENHANCED = Path("data/outputs/architecture/pfam_architecture_enhanced.tsv")
CANDIDATES = Path("data/outputs/orphans/orphan_candidates.tsv")
FASTA_ALL = Path("data/afdb/index/swissprot_sequences.fasta")

# Outputs
OUT_DIR = Path("data/outputs/orphans")
FASTA_OUT = OUT_DIR / "orphan_sequences.fasta"
META_TSV  = OUT_DIR / "orphan_metadata.tsv"
MISSING_TXT = OUT_DIR / "missing_sequences.txt"
SUMMARY_JSON = OUT_DIR / "orphan_sets_summary.json"

def to_int(s: str) -> int:
    try: return int(s)
    except Exception: return 0

def to_float(s: str) -> float:
    try: return float(s)
    except Exception: return 0.0

def largest_gap_len(gap_lengths: str) -> int:
    best = 0
    for tok in (gap_lengths or "").split(";"):
        tok = tok.strip()
        if not tok: 
            continue
        try:
            best = max(best, int(tok))
        except Exception:
            pass
    return best

def load_candidates(path: Path) -> List[str]:
    # Keep order (sorted by lowest coverage in 4A)
    order: List[str] = []
    with path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            acc = r.get("uniprot","").strip()
            if acc:
                order.append(acc)
    return order

def load_enhanced(path: Path) -> Dict[str, dict]:
    m: Dict[str, dict] = {}
    with path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            acc = r["uniprot"].strip()
            if not acc: 
                continue
            # Keep only fields we need
            L  = to_int(r.get("residues","0"))
            cov_res = to_int(r.get("covered_residues","0"))
            cov = to_float(r.get("coverage_frac","0"))
            n_dom = to_int(r.get("n_domains","0"))
            gaps = r.get("gaps","")
            gap_lengths = r.get("gap_lengths","")
            m[acc] = {
                "uniprot": acc,
                "residues": L,
                "covered_residues": cov_res,
                "coverage_frac": cov,
                "unannotated_residues": max(L - cov_res, 0),
                "unannotated_frac": max(1.0 - cov, 0.0),
                "n_domains": n_dom,
                "largest_gap": largest_gap_len(gap_lengths),
                "domains_compact": r.get("domains_compact",""),
                "gaps": gaps,
                "gap_lengths": gap_lengths,
            }
    return m

def collect_sequences(subset: set[str], fasta_path: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    # Robust FASTA reader without external libs
    with fasta_path.open() as fh:
        acc = None
        chunks: List[str] = []
        for ln in fh:
            if ln.startswith(">"):
                if acc in subset and chunks:
                    seqs[acc] = "".join(chunks)
                # header like >A0A023FBW4
                acc = ln[1:].strip().split()[0]
                chunks = []
            else:
                if acc in subset:
                    chunks.append(ln.strip())
        # flush last
        if acc in subset and chunks:
            seqs[acc] = "".join(chunks)
    return seqs

def main() -> None:
    assert ENHANCED.exists(), f"Missing {ENHANCED}"
    assert CANDIDATES.exists(), f"Missing {CANDIDATES}"
    assert FASTA_ALL.exists(), f"Missing {FASTA_ALL}"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    order = load_candidates(CANDIDATES)
    meta_all = load_enhanced(ENHANCED)

    # Build metadata TSV for candidates, preserving order
    with META_TSV.open("w") as fh:
        fh.write("\t".join([
            "uniprot","residues","covered_residues","coverage_frac",
            "unannotated_residues","unannotated_frac","n_domains",
            "largest_gap","domains_compact","gaps","gap_lengths"
        ]) + "\n")
        for acc in order:
            m = meta_all.get(acc)
            if not m:
                # Candidate without enhanced row (should not happen)
                continue
            fh.write("\t".join([
                m["uniprot"], str(m["residues"]), str(m["covered_residues"]),
                f"{m['coverage_frac']:.4f}",
                str(m["unannotated_residues"]), f"{m['unannotated_frac']:.4f}",
                str(m["n_domains"]), str(m["largest_gap"]),
                m["domains_compact"], m["gaps"], m["gap_lengths"]
            ]) + "\n")

    # Collect sequences for candidates (one pass)
    subset = set(order)
    seqs = collect_sequences(subset, FASTA_ALL)

    # Report missing and write FASTA in candidate order
    missing = [acc for acc in order if acc not in seqs]
    with MISSING_TXT.open("w") as fh:
        for acc in missing:
            fh.write(acc + "\n")

    with FASTA_OUT.open("w") as fh:
        for acc in order:
            if acc in seqs:
                fh.write(f">{acc}\n")
                s = seqs[acc]
                for i in range(0, len(s), 60):
                    fh.write(s[i:i+60] + "\n")

    # Summary
    SUMMARY_JSON.write_text(json.dumps({
        "candidates": len(order),
        "with_sequence": len(seqs),
        "missing_sequence": len(missing),
        "fraction_missing": round((len(missing) / len(order)) if order else 0.0, 6)
    }, indent=2))
    print(f"Wrote {META_TSV}")
    print(f"Wrote {FASTA_OUT}")
    print(f"Wrote {MISSING_TXT}")
    print(f"Wrote {SUMMARY_JSON}")

if __name__ == "__main__":
    main()
