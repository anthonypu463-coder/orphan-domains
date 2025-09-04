#!/usr/bin/env bash
set -euo pipefail
ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
ACC_LIST="${ACC_LIST:-$ROOT/data/afdb/index/swissprot_accessions.txt}"
PAE_DIR="${PAE_DIR:-$ROOT/data/afdb/pae}"
INTERVAL="${INTERVAL:-5}"

[ -s "$ACC_LIST" ] || { echo "Missing $ACC_LIST"; exit 1; }
mkdir -p "$PAE_DIR"
total=$(wc -l < "$ACC_LIST" | tr -d '[:space:]')

while :; do
  clear || true
  now=$(date -Iseconds)
  done_count=$(find "$PAE_DIR" -maxdepth 1 -type f -name '*.json' 2>/dev/null | wc -l | tr -d '[:space:]')
  tmp_count=$(find "$PAE_DIR" -maxdepth 1 -type f -name '*.tmp' 2>/dev/null | wc -l | tr -d '[:space:]')
  remaining=$(( total - done_count ))
  pct=$(awk -v d="$done_count" -v t="$total" 'BEGIN{ if(t>0) printf "%.2f", (d*100)/t; else print "0.00"}')
  echo "AFDB PAE download progress — $now"
  printf "Completed: %d / %d (%s%%)\nPending:   %d\nTemp files: %d\n" "$done_count" "$total" "$pct" "$remaining" "$tmp_count"
  pgrep -f "scripts/fetch_afdb_pae.sh" >/dev/null 2>&1 && echo "Fetcher: RUNNING" || echo "Fetcher: not detected"
  active_wget=$(pgrep -fa "wget.*alphafold\.ebi\.ac\.uk/files/AF-.*predicted_aligned_error_v4\.json" 2>/dev/null | wc -l | tr -d '[:space:]')
  echo "Active wget: $active_wget"
  if (( remaining > 0 )); then
    echo "Examples of remaining accessions:"
    if (( done_count > 0 )); then
      comm -23 <(sort "$ACC_LIST") <(find "$PAE_DIR" -maxdepth 1 -type f -name '*.json' -printf '%f\n' | sed 's/\.json$//' | sort) | head -n 5
    else
      head -n 5 "$ACC_LIST"
    fi
  fi
  sleep "$INTERVAL"
done
