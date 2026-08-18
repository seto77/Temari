#!/bin/bash
# run_cert_fleet.sh [レーン数] [スレッド数] [接頭辞] — σ(β,Δ) 全格子認証のフリート起動
#                                                      (260819Cl)
#
# 16 物理コア / 32 論理に **8 プロセス × 4 スレッド**で敷き詰める。
# 単一プロセス 32 スレッドだと、既定の窓が GL 16 点しか無いので半分のスレッドが遊ぶ
# (オラクルの 48 点も 2 波)。1 行あたりのノード評価は 16 窓 × (16+48) = 1024 で
# 固定なので、細かいプロセスに割った方が理想 (32 ノード時間/行) に近づく。
#
# ⚠ 各レーンは**別の JSONL** へ書く。集計は
#   julia +1.11 --project=. tools/certify_sigma.jl ../cert_sigma_v1_lane*.jsonl --summary
#
# ⚠ 出力はリポの**外** (`../`) に置く — 505 MB の refs と同じで、公開リポに
#   中間生成物を入れないため。
set -u
nlane=${1:-8}
nthr=${2:-4}
prefix=${3:-cert_sigma_v1}
cd "$(dirname "$0")/.." || exit 1

for lane in $(seq 0 $((nlane - 1))); do
  nohup bash tools/cert_watchdog.sh "$lane" "$nlane" "$nthr" +1.11 "$prefix" \
        > "../${prefix}_lane${lane}_watchdog.txt" 2>&1 &
  echo "lane $lane 起動 (pid $!)"
  sleep 2
done
echo "全 $nlane レーン起動。生死は ../${prefix}_lane*.jsonl の mtime で見る"
