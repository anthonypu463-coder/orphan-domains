#!/usr/bin/env bash
set -euo pipefail
ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
FASTA="${1:-$ROOT/data/afdb/index/swissprot_sequences.fasta}"
OUTDIR="${2:-$ROOT/data/outputs/hmmscan/shards}"
SHARDS="${SHARDS:-16}"

[ -s "$FASTA" ] || { echo "Missing FASTA: $FASTA"; exit 1; }
mkdir -p "$OUTDIR"
rm -f "$OUTDIR"/shard_*.faa

# Round-robin split by records to balance lengths approximately
awk -v S="$SHARDS" -v OUT="$OUTDIR" '
  BEGIN{rec=0}
  /^>/{rec++; idx=(rec-1)%S; fn=sprintf("%s/shard_%05d.faa", OUT, idx); print > fn; next}
  {print >> fn}
' "$FASTA"

ls -1 "$OUTDIR"/shard_*.faa | wc -l
