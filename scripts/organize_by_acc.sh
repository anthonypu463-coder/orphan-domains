#!/usr/bin/env bash
# Create by-acc symlinks in parallel. Idempotent.
set -euo pipefail
ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
ACC_LIST="$ROOT/data/afdb/index/swissprot_accessions.txt"
MANIFEST="$ROOT/data/afdb/index/swissprot_manifest.tsv"
BYACC="$ROOT/data/afdb/by_acc"
PAE_DIR="$ROOT/data/afdb/pae"
WORK="$ROOT/data/afdb/index/byacc_work.tsv"
PARALLEL="${PARALLEL:-8}"

[ -s "$ACC_LIST" ] || { echo "Missing $ACC_LIST"; exit 1; }
[ -s "$MANIFEST" ] || { echo "Missing $MANIFEST"; exit 1; }
mkdir -p "$BYACC"

# Join acc -> cif path into a worklist: "<ACC>\t<cif_path>"
awk 'NR==FNR && FNR>1 {m[$1]=$2; next} {if(length($0) && ($0 in m)) print $0 "\t" m[$0]}' "$MANIFEST" "$ACC_LIST" > "$WORK"

link_one() {
  acc="$1"; cif="$2"
  base="$(basename "$cif")"
  d="$BYACC/$acc"
  mkdir -p "$d"
  # Relative symlinks to avoid expensive path resolution
  ln -sf "../cif_gz/$base" "$d/model_v4.cif.gz"
  if [ -s "$PAE_DIR/$acc.json" ]; then
    ln -sf "../pae/$acc.json" "$d/pae.json"
  fi
}
export -f link_one
export BYACC PAE_DIR

xargs -a "$WORK" -n2 -P"$PARALLEL" bash -c 'link_one "$@"' _
echo "Organized by accession at $BYACC"
