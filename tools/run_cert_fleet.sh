#!/bin/bash
# run_cert_fleet.sh [レーン数] [スレッド数] [接頭辞] — σ(β,Δ) 全格子認証のフリート起動
#                                                      (260819Cl)
#
# ## ⚠⚠ 既定を 8×4 から **4×3 へ落とした** (260819Cl、実測を受けて)
#
# 8 レーン × 4 スレッド (= 論理 32 本を使い切る) で起動したところ、**1 行も書けないうちに
# ホストが BSOD した** (2026-08-18 14:05 頃、bugcheck 0x139 param1=3 = LIST_ENTRY 破損)。
# このホストは 4×32 GB DDR5 を AM5 の 4 スロット全部に挿しており、直近 1 か月で
# **10 回**の BSOD 履歴を持つ (0x4E PFN_LIST_CORRUPT / 0x101 CLOCK_WATCHDOG_TIMEOUT /
# 0xBE / 0x1E c0000005 ×3 …)。**⇒ 負荷は下げる。**
# 詳細と、これが「Julia の GC クラッシュ」の説明になりうることは
# `docs/host_stability_2026-08-19.md`。
#
#   - **同時スレッドを 12 本に**する (論理 32 の 37 %、物理 16 の 75 %)
#   - **起動をずらす** — 8 プロセスが同時に JIT を回すのが急性の引き金。
#     実際、落ちたのは全レーンがコンパイル中の時間帯だった
#   - **nice -n 10** で対話操作を巻き込まない
#
# ## ★ 12 本の内訳は **12 レーン × 1 スレッド** (260819Cl、実測で決めた)
#
# 1 行の所要を実測すると (どちらも他が 12 スレッド走っている状態):
#
#     1 スレッド 60.5 s/行 /  3 スレッド 25.6 s/行
#     ⇒ Amdahl で 直列 8.1 s + 並列 52.4 s
#
# 同じ 12 スレッドをどう割るかで実効が変わる:
#
#     4 レーン × 3 スレッド : 25.6/4 = 6.40 s/行 → 26.3 時間
#     6 レーン × 2 スレッド : 34.3/6 = 5.72 s/行 → 23.5 時間
#     **12 レーン × 1 スレッド : 60.5/12 = 5.04 s/行 → 20.7 時間**
#
# **CPU の占有は同じ 12 本**のまま 21 % 速い。総 CPU 仕事量が一定なので、
# 細かく割るほど Amdahl の損が消える。
# ⚠ **レーン数の変更で既済の行は失われない** — 済み判定は
# `<接頭辞>_lane*.jsonl` を全部読む (`sibling_lane_files`)。
# ⚠ 直列 8.1 s の正体は行ごとの `prepare_channel` (束縛状態は E₀ に依らないのに
#   E₀ 行ごとに解き直している)。**これを直せばさらに 2.6 時間縮むが、
#   `certify_sigma.jl` を編集すると走行中の認証の指紋が変わって全行が無効になる**
#   ので、走り始めてからは触らない (`docs/certification_2026-08-19.md` §2.1)。
#
# ⚠ 各レーンは**別の JSONL** へ書く。集計は
#   julia +1.11 --project=. tools/certify_sigma.jl ../cert_sigma_v1_lane*.jsonl --summary
#
# ⚠ 出力はリポの**外** (`../`) に置く — 505 MB の refs と同じで、公開リポに
#   中間生成物を入れないため。
set -u
nlane=${1:-12}
nthr=${2:-1}
prefix=${3:-cert_sigma_v1}
stagger=${STAGGER:-60}
cd "$(dirname "$0")/.." || exit 1

echo "フリート起動: $nlane レーン × $nthr スレッド = $((nlane * nthr)) 本 (論理 32 / 物理 16)"
for lane in $(seq 0 $((nlane - 1))); do
  nohup nice -n 10 bash tools/cert_watchdog.sh "$lane" "$nlane" "$nthr" +1.11 "$prefix" \
        > "../${prefix}_lane${lane}_watchdog.txt" 2>&1 &
  echo "  lane $lane 起動 (pid $!)  $(date '+%F %T')"
  [ "$lane" -lt $((nlane - 1)) ] && sleep "$stagger"
done
echo "全 $nlane レーン起動。生死は ../${prefix}_lane*.jsonl の mtime で見る"
