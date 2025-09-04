#!/usr/bin/env bash
set -euo pipefail
ROOT="$(git -C ~/orphan-domains rev-parse --show-toplevel 2>/dev/null || echo ~/orphan-domains)"
ACC_LIST="$ROOT/data/afdb/index/swissprot_accessions.txt"
OUT="$ROOT/data/afdb/index/reference_manifest.tsv"
INTERVAL="${INTERVAL:-5}"

[ -s "$ACC_LIST" ] || { echo "Missing $ACC_LIST"; exit 1; }

total=$(wc -l < "$ACC_LIST" | tr -d '[:space:]')
prev_done=0
prev_t=$(date +%s)
ema_rate=0  # rows/sec, simple EMA

fmt_eta() {
  local secs="$1"
  if (( secs < 0 )); then secs=0; fi
  local d=$((secs/86400)); local h=$(( (secs%86400)/3600 ))
  local m=$(( (secs%3600)/60 )); local s=$(( secs%60 ))
  if (( d>0 )); then printf "%dd %02d:%02d:%02d" "$d" "$h" "$m" "$s"
  else printf "%02d:%02d:%02d" "$h" "$m" "$s"; fi
}

while :; do
  clear || true
  now=$(date +%s)

  if [ -f "$OUT" ]; then
    done_rows=$(( $(wc -l < "$OUT") - 1 ))
    (( done_rows < 0 )) && done_rows=0
  else
    done_rows=0
  fi

  # instantaneous rate over last interval
  dt=$(( now - prev_t ))
  drows=$(( done_rows - prev_done ))
  if (( dt > 0 )); then
    inst_rate=$(awk -v a="$drows" -v b="$dt" 'BEGIN{ if(b>0) printf "%.6f", a/b; else print 0 }')
  else
    inst_rate=0
  fi

  # exponential moving average to smooth ETA (alpha=0.3)
  ema_rate=$(awk -v ema="$ema_rate" -v inst="$inst_rate" 'BEGIN{ printf "%.6f", (0.7*ema + 0.3*inst) }')

  remaining=$(( total - done_rows ))
  pct=$(awk -v d="$done_rows" -v t="$total" 'BEGIN{ if(t>0) printf "%.2f", (d*100)/t; else print "0.00" }')

  if awk -v r="$ema_rate" 'BEGIN{ exit !(r>0) }'; then
    eta_secs=$(awk -v rem="$remaining" -v r="$ema_rate" 'BEGIN{ if(r>0){printf "%.0f", rem/r}else{print -1} }')
    eta_human=$(fmt_eta "${eta_secs:-0}")
  else
    eta_human="estimating…"
  fi

  echo "AFDB reference manifest build — $(date -Iseconds)"
  printf "Completed: %'d / %'d (%s%%)\n" "$done_rows" "$total" "$pct"
  printf "Δrows (last %ss): %'d\n" "$INTERVAL" "$drows"
  printf "Rate (smoothed): %.2f rows/sec\n" "$ema_rate"
  printf "ETA: %s\n" "$eta_human"

  if pgrep -f "$ROOT/scripts/build_reference_manifest.py" >/dev/null 2>&1; then
    echo "Builder: RUNNING"
    ps -o pid,pcpu,pmem,etime,comm,args -p "$(pgrep -f -n "$ROOT/scripts/build_reference_manifest.py")" | sed -e '1,1!s/^/  /'
  else
    echo "Builder: not detected"
  fi

  if [ -f "$OUT" ]; then
    echo
    echo "Last updated: $(date -r "$OUT" -Iseconds 2>/dev/null || stat -c %y "$OUT")"
    echo "Tail (3 rows):"
    tail -n 3 "$OUT" | sed 's/^/  /'
  fi

  prev_done=$done_rows
  prev_t=$now
  sleep "$INTERVAL"
done
