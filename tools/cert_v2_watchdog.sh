#!/bin/bash
# cert_v2_watchdog.sh <レーン番号> <レーン総数> <スレッド数> <profile> [出力の接頭辞] [juliaup チャネル]
#     — σ(β,Δ) 本番候補の認証 v2 (`tools/certify_sigma_v2.jl`) を見張り付きで回す (260819Cl)
#
# `tools/cert_watchdog.sh` (v1) と同じ構え: **JSONL の mtime** で生死を見る (行ごとに flush)。
# ⚠ v2 は 1 行に 18 窓あり、1 行 = 最大 ~1 時間 (1 スレッド・n_q=1216) かかる。
#   そのため STALL は 1 行の最長より十分長く取る (行の途中は JSONL が動かない)。
#   ⇒ 窓ごとには書かない (行の全窓が揃ってから書く = 再開の単位が行) ので、
#      STALL = 7200 s (2 時間)。健全なレーンは ≤ 1 時間に 1 回書く。
# ⚠ **完走 ≠ 健全。**終わったら必ず `--summary` を読むこと。
set -u
lane=${1:?レーン番号 (0 始まり)}
nlane=${2:?レーン総数}
nthr=${3:-1}
profile=${4:-pilot}
prefix=${5:-cert_v2_${profile}}
chan=${6:-+1.11}
cd "$(dirname "$0")/.." || exit 1
out="$PWD/../${prefix}_lane${lane}.jsonl"
log="$PWD/../${prefix}_lane${lane}_log.txt"
touch "$out"
STALL=7200
for attempt in $(seq 1 200); do
  echo "=== cert v2 lane $lane/$nlane ($profile) attempt $attempt start julia$chan -t $nthr $(date '+%F %T') ===" >> "$log"
  touch "$out"
  julia $chan --project=. -t "$nthr" --gcthreads=1 tools/certify_sigma_v2.jl \
        "$out" --profile "$profile" --lane "$lane/$nlane" >> "$log" 2>&1 &
  jpid=$!
  while kill -0 $jpid 2>/dev/null; do
    sleep 60
    now=$(date +%s)
    mt=$(stat -c %Y "$out" 2>/dev/null || echo "$now")
    stall=$((now - mt))
    if [ $stall -gt $STALL ]; then
      echo "=== watchdog: jsonl stalled ${stall}s, killing pid $jpid $(date '+%F %T') ===" >> "$log"
      kill -9 $jpid 2>/dev/null
      sleep 10
      break
    fi
  done
  wait $jpid 2>/dev/null
  rc=$?
  echo "=== cert v2 lane $lane/$nlane attempt $attempt exit=$rc $(date '+%F %T') ===" >> "$log"
  if [ $rc -eq 0 ]; then
    echo "=== cert v2 lane $lane/$nlane COMPLETE $(date '+%F %T') ===" >> "$log"
    exit 0
  fi
  sleep 10
done
