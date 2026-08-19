#!/bin/bash
# run_cert_v2_fleet.sh <profile> [レーン数] [スレッド数] [接頭辞] — 認証 v2 のフリート起動 (260819Cl)
#
#   bash tools/run_cert_v2_fleet.sh pilot 11 1        # sentinel 11 行 (時間の実測)
#   bash tools/run_cert_v2_fleet.sh deep 16 1         # 全チャネル × E₀ 3 点 (+ sentinel) — 数日
#   julia +1.11 --project=. tools/certify_sigma_v2.jl ../cert_v2_deep_lane*.jsonl --summary
#
# 負荷の方針は v1 と同じ (`tools/run_cert_fleet.sh` 冒頭): 論理 32 本を使い切らない、
# 起動をずらす (同時 JIT を避ける)、nice -n 10。レーン分割はチャネル単位。
# ⚠ 走行中は tools/{sigma_beta_delta,angular_split_v2,angular_sweep,beta_spike,certify_sigma_v2}.jl を
#   1 byte も触らない (CERT_FP_V2 が変わり、済みの行が捨てられる)。
# ⚠ 出力はリポの**外** (`../`)。
set -u
profile=${1:?profile (pilot|deep|sentinel)}
nlane=${2:-12}
nthr=${3:-1}
prefix=${4:-cert_v2_${profile}}
stagger=${STAGGER:-45}
cd "$(dirname "$0")/.." || exit 1
echo "認証 v2 フリート起動: profile=$profile  $nlane レーン × $nthr スレッド = $((nlane * nthr)) 本  接頭辞 $prefix"
for lane in $(seq 0 $((nlane - 1))); do
  nohup nice -n 10 bash tools/cert_v2_watchdog.sh "$lane" "$nlane" "$nthr" "$profile" "$prefix" +1.11 \
        > "../${prefix}_lane${lane}_watchdog.txt" 2>&1 &
  echo "  lane $lane 起動 (pid $!)  $(date '+%F %T')"
  [ "$lane" -lt $((nlane - 1)) ] && sleep "$stagger"
done
echo "全 $nlane レーン起動。生死は ../${prefix}_lane*.jsonl の mtime で見る"
