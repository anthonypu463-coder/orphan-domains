#!/usr/bin/env bash
# Cluster the primary set (0% Pfam & mean pLDDT >=70) with MMseqs2 linclust.
# Thresholds: min-seq-id=0.4, cov-mode=5, -c 0.8 (coverage of the shorter sequence).
set -euo pipefail
ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"

IN_FASTA="${1:-$ROOT/data/outputs/orphans/orphan_zero_struct70.fasta}"
DB_DIR="$ROOT/data/outputs/cluster/db"
TMP_DIR="$ROOT/data/outputs/cluster/tmp"
SEQDB="$DB_DIR/orphanZ70_DB"
CLUS="$DB_DIR/orphanZ70_clu"
REPDB="$DB_DIR/orphanZ70_rep"
OUT_TSV="$ROOT/data/outputs/cluster/orphanZ70_clusters.tsv"
OUT_REP_FASTA="$ROOT/data/outputs/cluster/orphanZ70_reps.fasta"
OUT_PARAMS="$ROOT/data/outputs/cluster/cluster_params.json"
OUT_STATS="$ROOT/data/outputs/cluster/cluster_stats.json"

THREADS="${THREADS:-16}"

mkdir -p "$DB_DIR" "$TMP_DIR" "$ROOT/data/outputs/cluster"

command -v mmseqs >/dev/null 2>&1 || { echo "mmseqs not found in PATH"; exit 1; }
[ -s "$IN_FASTA" ] || { echo "Missing input FASTA: $IN_FASTA"; exit 1; }

mmseqs createdb "$IN_FASTA" "$SEQDB"

mmseqs linclust "$SEQDB" "$CLUS" "$TMP_DIR" \
  --min-seq-id 0.4 -c 0.8 --cov-mode 5 --threads "$THREADS"

mmseqs createtsv "$SEQDB" "$SEQDB" "$CLUS" "$OUT_TSV"

mmseqs result2repseq "$SEQDB" "$CLUS" "$REPDB"
mmseqs convert2fasta "$REPDB" "$OUT_REP_FASTA"

python - <<'PY'
import json, os, pathlib
root = pathlib.Path(__file__).resolve().parents[1]
(out := root / "cluster" / "cluster_params.json").write_text(json.dumps({
  "input_fasta": str((root.parent / "orphans" / "orphan_zero_struct70.fasta").resolve()),
  "algo": "mmseqs linclust",
  "min_seq_id": 0.4,
  "cov": 0.8,
  "cov_mode": 5,
  "threads": int(os.environ.get("THREADS","16"))
}, indent=2))
print(f"Wrote {out}")
PY

python - <<'PY'
from pathlib import Path
import json
root = Path(__file__).resolve().parents[1]
tsv = root / "cluster" / "orphanZ70_clusters.tsv"
fasta = root.parent / "orphans" / "orphan_zero_struct70.fasta"
repfa = root / "cluster" / "orphanZ70_reps.fasta"
report = root / "cluster" / "cluster_stats.json"

rep_to_members = {}
with tsv.open() as fh:
    for ln in fh:
        ln = ln.strip()
        if not ln: continue
        rep, mem = ln.split("\t")[:2]
        rep_to_members.setdefault(rep, set()).add(mem)

def fasta_ids(path):
    ids=[]
    with path.open() as fh:
        for ln in fh:
            if ln.startswith(">"):
                ids.append(ln[1:].strip().split()[0])
    return ids

all_ids = fasta_ids(fasta)
N_in = len(all_ids)
members = set().union(*rep_to_members.values()) if rep_to_members else set()
N_mapped = len(members)
N_rep = len(rep_to_members)
sizes = sorted((len(s) for s in rep_to_members.values()), reverse=True)
singletons = sum(1 for k in sizes if k==1)
top = []
for rep, s in sorted(rep_to_members.items(), key=lambda kv: len(kv[1]), reverse=True)[:5]:
    top.append({"rep": rep, "size": len(s)})

report.write_text(json.dumps({
    "input_sequences": N_in,
    "clusters": N_rep,
    "mapped_members": N_mapped,
    "singletons": singletons,
    "largest_cluster_size": sizes[0] if sizes else 0,
    "median_cluster_size": sizes[len(sizes)//2] if sizes else 0,
    "top_clusters": top
}, indent=2))
print("Wrote", report)
PY

echo "Done. Mapping: $OUT_TSV  | Reps: $OUT_REP_FASTA"
