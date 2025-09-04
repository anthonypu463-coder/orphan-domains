#!/usr/bin/env bash
# Progress monitor for:
#  1) scripts/organize_by_acc.sh  → data/afdb/by_acc/<ACC>/…
#  2) scripts/build_reference_manifest.py → data/afdb/index/reference_manifest.tsv
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
ACC_LIST="${ACC_LIST:-$ROOT/data/afdb/index/swissprot_accessions.txt}"
BYACC="${BYACC:-$ROOT/data/afdb/by_acc}"
MANIFEST_OUT="${MANIFEST_OUT:-$ROOT/data/afdb/index/reference_manifest.tsv}"
INTERVAL="${INTERVAL:-10}"
BAR=30
FAST="${FAST:-0}"   # FAST=1 to skip symlink counts for lighter scans

[[ -s "$ACC_LIST" ]] || { echo "Missing accession list: $ACC_LIST"; exit 1; }
mkdir -p "$BYACC"
total=$(wc -l < "$ACC_LIST" | tr -d '[:space:]')

prev_dirs=0 prev_rows=0 prev_ts=$(date +%s)

while :; do
  now=$(date +%s)

  # by-acc progress
  done_dirs=$(find "$BYACC" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l | tr -d '[:space:]')
  if [[ "$FAST" -eq 1 ]]; then
    linked_models="-" ; linked_pae="-"
  else
    linked_models=$(find "$BYACC" -mindepth 2 -maxdepth 2 -type l -name 'model_v4.cif.gz' 2>/dev/null | wc -l | tr -d '[:space:]')
    linked_pae=$(find "$BYACC" -mindepth 2 -maxdepth 2 -type l -name 'pae.json' 2>/dev/null | wc -l | tr -d '[:space:]')
  fi

  # manifest progress
  if [[ -f "$MANIFEST_OUT" ]]; then
    rows=$(($(wc -l < "$MANIFEST_OUT" | tr -d '[:space:]') - 1))
    (( rows < 0 )) && rows=0
  else
    rows=0
  fi

  pct_dirs=$(awk -v a="$done_dirs" -v b="$total" 'BEGIN{ if(b>0) printf "%.2f", (a*100)/b; else print "0.00"}')
  pct_rows=$(awk -v a="$rows" -v b="$total" 'BEGIN{ if(b>0) printf "%.2f", (a*100)/b; else print "0.00"}')

  dt=$(( now - prev_ts )); (( dt<=0 )) && dt=1
  prev_dirs=$done_dirs; prev_rows=$rows; prev_ts=$now

  fill_dirs=$(( (done_dirs * BAR) / (total>0?total:1) ))
  fill_rows=$(( (rows * BAR) / (total>0?total:1) ))
  bar_dirs="$(printf '%*s' "$fill_dirs" '' | tr ' ' '#')$(printf '%*s' "$((BAR - fill_dirs))" '' | tr ' ' '-')"
  bar_rows="$(printf '%*s' "$fill_rows" '' | tr ' ' '#')$(printf '%*s' "$((BAR - fill_rows))" '' | tr ' ' '-')"

  clear || true
  echo "Monitor — by-acc organization & reference manifest    $(date -Iseconds)"
  echo
  pgrep -fa "organize_by_acc\.sh" >/dev/null 2>&1 && echo "organize_by_acc.sh: RUNNING" || echo "organize_by_acc.sh: not detected"
  pgrep -fa "build_reference_manifest\.py" >/dev/null 2>&1 && echo "build_reference_manifest.py: RUNNING" || echo "build_reference_manifest.py: not detected"
  echo
  printf "by_acc dirs:   [%s] %s%%  (%d / %d)\n" "$bar_dirs" "$pct_dirs" "$done_dirs" "$total"
  printf "symlinks:      models=%s  pae=%s\n" "$linked_models" "$linked_pae"
  echo
  printf "manifest rows: [%s] %s%%  (%d / %d)\n" "$bar_rows" "$pct_rows" "$rows" "$total"
  [[ -f "$MANIFEST_OUT" ]] && tail -n 2 "$MANIFEST_OUT" | sed 's/^/last rows: /'
  echo
  echo "Update every ${INTERVAL}s (FAST=$FAST). Ctrl-C to exit."
  sleep "$INTERVAL"

  if [[ "$done_dirs" -ge "$total" && "$rows" -ge "$total" ]]; then
    if ! pgrep -fa "organize_by_acc\.sh" >/dev/null 2>&1 && ! pgrep -fa "build_reference_manifest\.py" >/dev/null 2>&1; then
      echo "All done."
      exit 0
    fi
  fi
done
