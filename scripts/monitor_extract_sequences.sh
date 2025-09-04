#!/usr/bin/env bash
set -euo pipefail
ROOT="$(git -C ~/orphan-domains rev-parse --show-toplevel 2>/dev/null || echo ~/orphan-domains)"
ACC_LIST="${ACC_LIST:-$ROOT/data/afdb/index/swissprot_accessions.txt}"
FASTA="${FASTA:-$ROOT/data/afdb/index/swissprot_sequences.fasta}"
INTERVAL="${INTERVAL:-5}"

[ -s "$ACC_LIST" ] || { echo "Missing $ACC_LIST"; exit 1; }
total=$(wc -l < "$ACC_LIST" | tr -d '[:space:]')

prev_done=0
prev_t=$(date +%s)
ema_rate=0

count_headers() { grep -a -c '^>' "$FASTA" 2>/dev/null || echo 0; }

while :; do
  clear || true
  now=$(date +%s)

  if [ -s "$FASTA" ]; then
    done=$(count_headers)
    size=$(stat -c %s "$FASTA" 2>/dev/null || echo 0)
    mtime=$(date -r "$FASTA" -Iseconds 2>/dev/null || stat -c %y "$FASTA")
  else
    done=0; size=0; mtime="(not created)"
  fi

  dt=$((now - prev_t)); drows=$((done - prev_done))
  if (( dt > 0 )); then
    inst=$(awk -v a="$drows" -v b="$dt" 'BEGIN{ if(b>0) printf "%.6f", a/b; else print 0 }')
  else
    inst=0
  fi
  ema_rate=$(awk -v ema="$ema_rate" -v inst="$inst" 'BEGIN{ printf "%.6f", (0.7*ema + 0.3*inst) }')
  pct=$(awk -v d="$done" -v t="$total" 'BEGIN{ if(t>0) printf "%.2f", (d*100)/t; else print "0.00"}')

  if awk -v r="$ema_rate" 'BEGIN{ exit !(r>0) }'; then
    eta_secs=$(awk -v rem="$((total - done))" -v r="$ema_rate" 'BEGIN{ if(r>0){printf "%.0f", rem/r}else{print -1} }')
    h=$((eta_secs/3600)); m=$(( (eta_secs%3600)/60 )); s=$((eta_secs%60))
    eta=$(printf "%02d:%02d:%02d" "$h" "$m" "$s")
  else
    eta="estimating…"
  fi

  echo "Sequence extraction monitor — $(date -Iseconds)"
  printf "Completed: %'d / %'d (%s%%)\n" "$done" "$total" "$pct"
  printf "Δrecords (last %ss): %'d\n" "$INTERVAL" "$drows"
  printf "Rate (smoothed): %.2f seqs/sec\n" "$ema_rate"
  printf "ETA: %s\n" "$eta"
  printf "FASTA: %s (size %'d bytes; last write %s)\n" "$FASTA" "$size" "$mtime"

  if pgrep -f "$ROOT/scripts/extract_sequences.py" >/dev/null 2>&1; then
    echo "Extractor: RUNNING"
  else
    echo "Extractor: not detected"
  fi

  echo
  echo "Recent headers:"
  tail -n 200 "$FASTA" 2>/dev/null | grep -a '^>' | tail -n 3 | sed 's/^/  /' || echo "  (no records yet)"

  prev_done=$done
  prev_t=$now
  sleep "$INTERVAL"
done
