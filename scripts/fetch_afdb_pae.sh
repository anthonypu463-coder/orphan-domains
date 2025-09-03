#!/usr/bin/env bash
set -euo pipefail

ACC_LIST="data/afdb/index/swissprot_accessions.txt"
DEST="data/afdb/pae"
PARALLEL="${PARALLEL:-8}"

test -s "$ACC_LIST" || { echo "Missing $ACC_LIST"; exit 1; }
mkdir -p "$DEST"

fetch_one() {
  acc="$1"
  out="${DEST}/${acc}.json"
  tmp="${out}.tmp"
  url="https://alphafold.ebi.ac.uk/files/AF-${acc}-F1-predicted_aligned_error_v4.json"
  if [ -s "$out" ]; then
    exit 0
  fi
  # Use wget for portability; retry a few times; move into place atomically on success
  if wget -q --tries=3 --timeout=30 -O "$tmp" "$url"; then
    mv "$tmp" "$out"
    echo "OK $acc"
  else
    rm -f "$tmp"
    echo "FAIL $acc" >&2
  fi
  # Gentle pacing to avoid hammering the server
  sleep 0.05
}

export -f fetch_one
export DEST
# xargs drives limited parallelism
xargs -a "$ACC_LIST" -n1 -P"$PARALLEL" bash -c 'fetch_one "$@"' _

# Summarize
ok=$(ls -1 "$DEST"/*.json 2>/dev/null | wc -l || echo 0)
echo "PAE JSON files present: $ok"
