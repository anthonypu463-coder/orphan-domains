#!/usr/bin/env bash
set -euo pipefail
ROOT="$(git -C ~/orphan-domains rev-parse --show-toplevel 2>/dev/null || echo ~/orphan-domains)"
SHARD_DIR="${SHARD_DIR:-$ROOT/data/outputs/hmmscan/shards}"
OUT_DIR="${OUT_DIR:-$ROOT/data/outputs/hmmscan/domtblout}"
LOG_DIR="${LOG_DIR:-$ROOT/data/outputs/hmmscan/logs}"
INTERVAL="${INTERVAL:-5}"

shopt -s nullglob
total() { ls -1 "$SHARD_DIR"/shard_*.faa 2>/dev/null | wc -l | tr -d '[:space:]'; }
done_parts() { ls -1 "$OUT_DIR"/shard_*.domtblout 2>/dev/null | wc -l | tr -d '[:space:]'; }
bases() { for f in "$1"/shard_*.${2}; do b=$(basename "$f"); echo "${b%%.*}"; done | sort -u; }

T=$(total)
(( T > 0 )) || { echo "No shards found in $SHARD_DIR"; exit 1; }

prev_done=0; prev_t=$(date +%s); ema_rate=0

fmt_eta(){ local s=$1; ((s<0)) && s=0; printf "%02d:%02d:%02d" $((s/3600)) $(((s%3600)/60)) $((s%60)); }

while :; do
  clear || true
  now=$(date +%s)
  D=$(done_parts)
  dt=$((now - prev_t)); dparts=$((D - prev_done))
  inst_rate=$(awk -v a="$dparts" -v b="$dt" 'BEGIN{ if(b>0) printf "%.6f", a/b; else print 0 }')
  ema_rate=$(awk -v ema="$ema_rate" -v inst="$inst_rate" 'BEGIN{ printf "%.6f", (0.7*ema + 0.3*inst) }')
  pct=$(awk -v d="$D" -v t="$T" 'BEGIN{ if(t>0) printf "%.2f", (d*100)/t; else print "0.00"}')
  rem=$((T - D))
  eta=$( if awk -v r="$ema_rate" 'BEGIN{exit !(r>0)}'; then awk -v rem="$rem" -v r="$ema_rate" 'BEGIN{printf "%.0f", rem/r}'; else echo -1; fi )
  echo "Pfam hmmscan monitor — $(date -Iseconds)"
  printf "Completed: %'d / %'d (%s%%)\n" "$D" "$T" "$pct"
  printf "Δshards (last %ss): %'d\n" "$INTERVAL" "$dparts"
  printf "Rate (smoothed): %.2f shards/sec\n" "$ema_rate"
  printf "ETA: %s\n" "$(fmt_eta "$eta")"

  if pgrep -f "$ROOT/scripts/run_hmmscan.sh" >/dev/null 2>&1; then
    echo "Runner: RUNNING"
  else
    echo "Runner: not detected"
  fi
  active=$(pgrep -fa "hmmscan .*Pfam-A\.hmm" 2>/dev/null | wc -l | tr -d '[:space:]')
  echo "Active hmmscan processes: $active"
  if (( active > 0 )); then
    ps -o pid,pcpu,pmem,etime,comm,args -C hmmscan | head -n 6
  fi

  echo
  echo "Recent logs:"
  ls -1t "$LOG_DIR"/shard_*.log 2>/dev/null | head -n 2 | while read -r L; do
    echo "  $(basename "$L")  (updated $(date -r "$L" +%T))"
    tail -n 3 "$L" | sed 's/^/    /'
  done

  # Show a few missing shard basenames (useful to confirm progress)
  echo
  echo "Example remaining shards:"
  comm -23 <(bases "$SHARD_DIR" faa) <(bases "$OUT_DIR" domtblout) | head -n 5 | sed 's/^/  /' || true

  prev_done=$D; prev_t=$now
  sleep "$INTERVAL"
done
