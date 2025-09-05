#!/usr/bin/env python3
from __future__ import annotations
import csv, json
from pathlib import Path
from typing import Dict, List, Tuple

# Inputs
MAP_TSV   = Path("data/outputs/cluster26k/strict26k_clusters.tsv")     # rep \t member
FASTA_IN  = Path("data/outputs/orphans/orphan_strict.fasta")           # sequences for 26k
MANIFEST  = Path("data/afdb/index/reference_manifest.tsv")             # has residues, plddt_mean

# Outputs
OUT_DIR   = Path("data/outputs/cluster26k")
REPS_TSV  = OUT_DIR / "cluster26k_representatives.tsv"
MEMB_TSV  = OUT_DIR / "cluster26k_members.tsv"
REPS_FASTA= OUT_DIR / "strict26k_reps_refined.fasta"
SUMMARY   = OUT_DIR / "rep_selection_summary.json"

def load_plddt_len(manifest: Path) -> Tuple[Dict[str, float], Dict[str, int]]:
    plddt: Dict[str, float] = {}
    length: Dict[str, int] = {}
    with manifest.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            acc = r["uniprot"].strip()
            if not acc: continue
            try:
                plddt[acc]  = float(r.get("plddt_mean","") or 0.0)
                length[acc] = int(r.get("residues","") or 0)
            except Exception:
                continue
    return plddt, length

def load_mapping(map_tsv: Path) -> Dict[str, List[str]]:
    # rep -> unique members; include the rep in its own cluster
    clusters: Dict[str, List[str]] = {}
    with map_tsv.open() as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln: continue
            rep, mem = ln.split("\t")[:2]
            if rep not in clusters:
                clusters[rep] = [rep]
            if mem not in clusters[rep]:
                clusters[rep].append(mem)
    return clusters

def load_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    with path.open() as fh:
        acc=None; buf=[]
        for ln in fh:
            if ln.startswith(">"):
                if acc is not None:
                    seqs[acc] = "".join(buf)
                acc = ln[1:].strip().split()[0]
                buf=[]
            else:
                if acc is not None:
                    buf.append(ln.strip())
        if acc is not None:
            seqs[acc] = "".join(buf)
    return seqs

def main() -> None:
    assert MAP_TSV.exists(), f"Missing {MAP_TSV}"
    assert FASTA_IN.exists(), f"Missing {FASTA_IN}"
    assert MANIFEST.exists(), f"Missing {MANIFEST}"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    plddt, length = load_plddt_len(MANIFEST)
    seqs = load_fasta(FASTA_IN)
    clusters = load_mapping(MAP_TSV)

    # Representative selection: highest pLDDT, tie by length, then lexicographic
    reps_rows: List[Dict[str, str]] = []
    memb_rows: List[Tuple[str, str]] = []
    n_missing_seq = 0
    n_missing_plddt = 0

    def score(acc: str) -> Tuple[int, float, int, str]:
        p = plddt.get(acc, None)
        L = length.get(acc, len(seqs.get(acc,"")))
        has_p = 1 if p is not None else 0
        return (has_p, p if p is not None else 0.0, L if L is not None else 0, acc)

    # Stable cluster IDs
    cluster_ids = []
    for i, rep in enumerate(sorted(clusters)):  # deterministic order
        cluster_ids.append((f"CLS{(i+1):06d}", rep))

    for cid, rep0 in cluster_ids:
        members = sorted(set(clusters.get(rep0, [rep0])))
        for m in members:
            memb_rows.append((cid, m))

        # choose best representative
        best = max(members, key=score)
        if best not in seqs:
            n_missing_seq += 1  # should be 0
            continue
        if best not in plddt:
            n_missing_plddt += 1

        reps_rows.append({
            "cluster_id": cid,
            "representative": best,
            "plddt_mean": f"{plddt.get(best, 0.0):.2f}",
            "length": str(length.get(best, len(seqs[best]))),
            "n_members": str(len(members)),
            "mmseqs_rep": rep0,
            "selected_by": "plddt_then_length",
            "members": ";".join(members),
        })

    # Write mappings
    with REPS_TSV.open("w") as fh:
        fh.write("\t".join(["cluster_id","representative","plddt_mean","length","n_members","mmseqs_rep","selected_by","members"]) + "\n")
        for r in reps_rows:
            fh.write("\t".join([r[k] for k in ["cluster_id","representative","plddt_mean","length","n_members","mmseqs_rep","selected_by","members"]]) + "\n")

    with MEMB_TSV.open("w") as fh:
        fh.write("\t".join(["cluster_id","member"]) + "\n")
        for cid, m in memb_rows:
            fh.write(f"{cid}\t{m}\n")

    # Write representative sequences FASTA
    with REPS_FASTA.open("w") as fh:
        for r in reps_rows:
            acc = r["representative"]
            s = seqs.get(acc, "")
            if not s: continue
            fh.write(f">{acc}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i+60] + "\n")

    # Summary
    summary = {
        "clusters": len(reps_rows),
        "total_members": len(memb_rows),
        "missing_rep_sequences": n_missing_seq,
        "missing_rep_plddt": n_missing_plddt,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2))
    print(f"Wrote {REPS_TSV}")
    print(f"Wrote {MEMB_TSV}")
    print(f"Wrote {REPS_FASTA}")
    print(f"Wrote {SUMMARY}")

if __name__ == "__main__":
    main()
