#!/usr/bin/env bash
set -euo pipefail
ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
IN_FASTA="${1:-$ROOT/data/outputs/orphans/orphan_strict.fasta}"
DB_DIR="$ROOT/data/outputs/cluster26k/db"
TMP_DIR="$ROOT/data/outputs/cluster26k/tmp"
SEQDB="$DB_DIR/strict26k_DB"
CLUS="$DB_DIR/strict26k_clu"
REPDB="$DB_DIR/strict26k_rep"
OUT_MAP="$ROOT/data/outputs/cluster26k/strict26k_clusters.tsv"
OUT_REP="$ROOT/data/outputs/cluster26k/strict26k_reps.fasta"
THREADS="${THREADS:-16}"

mkdir -p "$DB_DIR" "$TMP_DIR" "$(dirname "$OUT_MAP")"
command -v mmseqs >/dev/null || { echo "mmseqs not found"; exit 1; }
[ -s "$IN_FASTA" ] || { echo "Missing input FASTA: $IN_FASTA"; exit 1; }

mmseqs createdb "$IN_FASTA" "$SEQDB"
mmseqs linclust "$SEQDB" "$CLUS" "$TMP_DIR" \
  --min-seq-id 0.4 -c 0.8 --cov-mode 5 --cluster-mode 1 --threads "$THREADS"
mmseqs createtsv "$SEQDB" "$SEQDB" "$CLUS" "$OUT_MAP"
mmseqs result2repseq "$SEQDB" "$CLUS" "$REPDB"
mmseqs convert2fasta "$REPDB" "$OUT_REP"
echo "Done: $OUT_MAP | $OUT_REP"
