#!/bin/bash
# lane_watchdog.sh <レーン番号> [レーン総数] [スレッド数] [juliaup チャネル] [tags] [出力先]
#                                        — 本番生成レーンの見張り付き実行 (260808Cl)
#
# ReciPro 側 `handout/lane_watchdog.sh` (v3 生成で使ったもの) を Temari へ移し、
# v4 用に引数の順序と既定値を整理した版。
#
# なぜ必要か: Julia は高割り当て・多スレッドの長時間バッチで
# EXCEPTION_ACCESS_VIOLATION を起こす (2026-08-04 実測: 1.12 は gc_mark_objarray、
# 1.11 は sweep_malloced_memory)。しかも**プロセスが死に切らず wedged になり、
# ログだけ止まる**ので、死活監視は「プロセス生存」ではなく「ログの mtime 停滞」で
# 行う必要がある。診断で分かっていること:
#   - RAM 126G 中 78G 空き・ページファイルほぼ未使用 → 資源枯渇ではない
#   - --gcthreads=1 で既に nmarkthreads=1 / nsweepthreads=0 (GC 並列は最小構成)
#   - --heap-size-hint は GC 回数を変えない
#   - ionization.jl に unsafe/ccall/pointer 演算はあるが @threads は互いに素な
#     添字への書き込みのみ → 自コードのデータ競合ではなくランタイム側の問題
# よって対策は「落ちたら即座に拾い直す」。出力はチャネル単位で原子的、さらに
# E0 行単位のチェックポイント (*.partial.jsonl) があるので損失は最大 1 行。
#
# ⚠ **完走 ≠ 健全。**v3 では GC クラッシュ由来のメモリ破損が
#    チェックポイント経由で 1 行だけ生き残った前例がある。QC を必ず通すこと。
set -u
lane=${1:?レーン番号 (0 始まり)}
nlane=${2:-8}
nthr=${3:-4}
chan=${4:-+1.11}
tags=${5:-}
outdir=${6:-}
cd "$(dirname "$0")/.." || exit 1
log="../temari_v4_lane${lane}_log.txt"
opts=""
[ -n "$tags" ] && opts="$opts --tags $tags"
[ -n "$outdir" ] && opts="$opts --out $outdir"
for attempt in $(seq 1 60); do
  echo "=== lane $lane/$nlane attempt $attempt start julia$chan -t $nthr$opts $(date '+%F %T') ===" >> "$log"
  julia $chan -t "$nthr" --gcthreads=1 src/gen_production.jl \
        --lane "$lane/$nlane" $opts >> "$log" 2>&1 &
  jpid=$!
  while kill -0 $jpid 2>/dev/null; do
    sleep 60
    now=$(date +%s)
    mt=$(stat -c %Y "$log" 2>/dev/null || echo "$now")
    if [ $((now - mt)) -gt 900 ]; then   # 15 分停滞 = wedged とみなす
      echo "=== watchdog: log stalled >15min, killing pid $jpid $(date '+%F %T') ===" >> "$log"
      kill -9 $jpid 2>/dev/null
      sleep 10
      break
    fi
  done
  wait $jpid 2>/dev/null
  rc=$?
  echo "=== lane $lane/$nlane attempt $attempt exit=$rc $(date '+%F %T') ===" >> "$log"
  if [ $rc -eq 0 ]; then
    echo "=== lane $lane/$nlane COMPLETE $(date '+%F %T') ===" >> "$log"
    exit 0
  fi
  sleep 10
done
echo "=== lane $lane/$nlane gave up after 60 attempts $(date '+%F %T') ===" >> "$log"
exit 1
