#!/usr/bin/env python3
from __future__ import annotations
import csv, json, shutil
from pathlib import Path
from typing import Dict, List

# Inputs from prior steps
META_ALL  = Path("data/outputs/orphans/orphan_metadata.tsv")        # 4B
FASTA_ALL = Path("data/outputs/orphans/orphan_sequences.fasta")     # 4B
STRUCT_OK = Path("data/outputs/orphans/orphan_structural_subset.txt")# 4C (pLDDT>=50)

# Outputs
OUT_DIR   = Path("data/outputs/orphans")
CAN_META  = OUT_DIR / "orphan_20pct_metadata.tsv"
CAN_FASTA = OUT_DIR / "orphan_20pct.fasta"
STR_META  = OUT_DIR / "orphan_strict_metadata.tsv"
STR_FASTA = OUT_DIR / "orphan_strict.fasta"
STRS_META = OUT_DIR / "orphan_strict_structural_metadata.tsv"
STRS_FASTA= OUT_DIR / "orphan_strict_structural.fasta"
REPORT    = OUT_DIR / "orphan_datasets_summary.json"

TARGET_STRICT = 26000  # ~26k target

def load_meta(path: Path) -> List[dict]:
    rows: List[dict] = []
    with path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            rows.append(r)
    return rows

def load_structural_ok(path: Path) -> set[str]:
    if not path.exists():
        return set()
    return set(ln.strip() for ln in path.read_text().splitlines() if ln.strip())

def collect_seqs(fasta: Path) -> Dict[str,str]:
    seqs: Dict[str,str] = {}
    with fasta.open() as fh:
        acc=None; chunks=[]
        for ln in fh:
            if ln.startswith(">"):
                if acc is not None:
                    seqs[acc] = "".join(chunks)
                acc = ln[1:].strip().split()[0]
                chunks=[]
            else:
                if acc is not None:
                    chunks.append(ln.strip())
        if acc is not None:
            seqs[acc] = "".join(chunks)
    return seqs

def write_meta(rows: List[dict], path: Path) -> None:
    if not rows:
        path.write_text("")
        return
    cols = list(rows[0].keys())
    with path.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r.get(c,"") for c in cols) + "\n")

def write_fasta(accs: List[str], seqs: Dict[str,str], path: Path) -> int:
    wrote = 0
    with path.open("w") as fh:
        for a in accs:
            s = seqs.get(a)
            if not s: continue
            fh.write(f">{a}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i+60] + "\n")
            wrote += 1
    return wrote

def main() -> None:
    assert META_ALL.exists() and FASTA_ALL.exists(), "Missing 4B outputs"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    meta_rows = load_meta(META_ALL)            # already sorted by lowest coverage
    seqs      = collect_seqs(FASTA_ALL)
    struct_ok = load_structural_ok(STRUCT_OK)

    # Canonical orphan-rich deliverable (<20%): copy as-is
    shutil.copy2(META_ALL, CAN_META)
    shutil.copy2(FASTA_ALL, CAN_FASTA)

    # Strict ~26k by lowest coverage
    strict_rows = meta_rows[:min(TARGET_STRICT, len(meta_rows))]
    strict_accs = [r["uniprot"] for r in strict_rows]
    write_meta(strict_rows, STR_META)
    n_strict = write_fasta(strict_accs, seqs, STR_FASTA)

    # Structural strict (intersection with pLDDT>=50, preserve order, cap at target)
    strict_struct_rows: List[dict] = []
    strict_struct_accs: List[str] = []
    for r in strict_rows:
        a = r["uniprot"]
        if a in struct_ok:
            strict_struct_rows.append(r)
            strict_struct_accs.append(a)
        if len(strict_struct_rows) >= TARGET_STRICT:
            break
    write_meta(strict_struct_rows, STRS_META)
    n_strict_struct = write_fasta(strict_struct_accs, seqs, STRS_FASTA)

    # Summaries
    def parse_float(s: str) -> float:
        try: return float(s)
        except Exception: return 0.0

    max_cov_all   = max((parse_float(r["coverage_frac"]) for r in meta_rows), default=0.0)
    max_cov_str   = max((parse_float(r["coverage_frac"]) for r in strict_rows), default=0.0)
    max_cov_str_s = max((parse_float(r["coverage_frac"]) for r in strict_struct_rows), default=0.0)

    REPORT.write_text(json.dumps({
        "canonical": {
            "metadata_rows": len(meta_rows),
            "fasta_seqs": len(seqs),
            "max_coverage_frac": round(max_cov_all, 6)
        },
        "strict": {
            "target": TARGET_STRICT,
            "metadata_rows": len(strict_rows),
            "fasta_seqs": n_strict,
            "max_coverage_frac": round(max_cov_str, 6)
        },
        "strict_structural": {
            "target": TARGET_STRICT,
            "metadata_rows": len(strict_struct_rows),
            "fasta_seqs": n_strict_struct,
            "max_coverage_frac": round(max_cov_str_s, 6)
        }
    }, indent=2))
    print(f"Wrote {CAN_META} and {CAN_FASTA}")
    print(f"Wrote {STR_META} ({len(strict_rows)}) and {STR_FASTA} ({n_strict})")
    print(f"Wrote {STRS_META} ({len(strict_struct_rows)}) and {STRS_FASTA} ({n_strict_struct})")
    print(f"Wrote {REPORT}")

if __name__ == "__main__":
    main()
