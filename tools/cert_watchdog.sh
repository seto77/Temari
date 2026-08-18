#!/bin/bash
# cert_watchdog.sh <レーン番号> [レーン総数] [スレッド数] [juliaup チャネル] [出力の接頭辞]
#                          — σ(β,Δ) 全格子認証レーンの見張り付き実行 (260819Cl)
#
# `tools/lane_watchdog.sh` (本番生成用) と同じ構えだが、**見張る対象が違う**:
#
#   ⚠⚠ ログではなく **JSONL の mtime** を見る。
#   `certify_sigma.jl` は進捗を **20 行ごと**にしか印字しないので、1 行 ~30 s なら
#   ログは 10 分間動かない = 15 分規則の余裕が 1.5 倍しか無い。JSONL は**行ごとに
#   flush** されるので、こちらが正しい生存信号になる。
#
# なぜ必要か: Windows Julia は長時間・高割り当てのバッチで EXCEPTION_ACCESS_VIOLATION
# を起こす (v4 生成中に 5 回、うち 2 回は **wedged** = プロセスが残りログだけ止まる)。
# 20 時間級の実行では**必ず遭遇する前提**で組む。
#
# 損失は最大 1 行 (JSONL は行ごとに flush、再開は既存行の読み飛ばし)。
# ⚠ **完走 ≠ 健全。**終わったら必ず `--summary` を読むこと。
set -u

lane=${1:?レーン番号 (0 始まり)}
nlane=${2:-8}
nthr=${3:-4}
chan=${4:-+1.11}
prefix=${5:-cert_sigma_v1}
cd "$(dirname "$0")/.." || exit 1

out="$PWD/../${prefix}_lane${lane}.jsonl"
log="$PWD/../${prefix}_lane${lane}_log.txt"
touch "$out"

STALL=900            # JSONL が 15 分動かなければ wedged とみなす
                     # (健全なレーンは 30-60 s ごとに 1 行書く = 15-30 倍の余裕)

for attempt in $(seq 1 200); do
  echo "=== cert lane $lane/$nlane attempt $attempt start julia$chan -t $nthr $(date '+%F %T') ===" >> "$log"
  touch "$out"       # 起動 (コンパイル) 中に誤判定しないよう mtime を打ち直す
  julia $chan --project=. -t "$nthr" --gcthreads=1 tools/certify_sigma.jl \
        "$out" --lane "$lane/$nlane" >> "$log" 2>&1 &
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
  echo "=== cert lane $lane/$nlane attempt $attempt exit=$rc $(date '+%F %T') ===" >> "$log"
  if [ $rc -eq 0 ]; then
    echo "=== cert lane $lane/$nlane COMPLETE $(date '+%F %T') ===" >> "$log"
    exit 0
  fi
  sleep 10
done
echo "=== cert lane $lane/$nlane gave up after 200 attempts $(date '+%F %T') ===" >> "$log"
exit 1
