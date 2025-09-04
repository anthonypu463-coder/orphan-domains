#!/usr/bin/env bash
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
HMMDB="${HMMDB:-$ROOT/data/pfam/Pfam-A.hmm}"
SHARD_DIR="$ROOT/data/outputs/hmmscan/shards"
OUT_DIR="$ROOT/data/outputs/hmmscan/domtblout"
LOG_DIR="$ROOT/data/outputs/hmmscan/logs"

JOBS="${JOBS:-8}"         # concurrent hmmscan processes
THREADS_PER="${THREADS_PER:-1}"  # --cpu for each hmmscan

command -v hmmscan >/dev/null 2>&1 || { echo "hmmscan not found"; exit 1; }
for suf in h3m h3i h3f h3p; do
  [ -s "${HMMDB}.${suf}" ] || { echo "Missing HMMER index: ${HMMDB}.${suf}"; exit 1; }
done

mkdir -p "$OUT_DIR" "$LOG_DIR"

# Build a runlist with one command per shard
RUNLIST="$(mktemp)"
trap 'rm -f "$RUNLIST"' EXIT

for shard in "$SHARD_DIR"/shard_*.faa; do
  base="$(basename "$shard" .faa)"
  dom="$OUT_DIR/${base}.domtblout"
  log="$LOG_DIR/${base}.log"
  echo "hmmscan --cut_ga --noali --cpu ${THREADS_PER} --domtblout '$dom' '$HMMDB' '$shard' > /dev/null 2> '$log'" >> "$RUNLIST"
done

# Execute with bounded parallelism
xargs -I CMD -P "$JOBS" bash -lc "CMD" < "$RUNLIST"

echo "hmmscan completed: $(ls -1 "$OUT_DIR"/*.domtblout | wc -l) parts"
