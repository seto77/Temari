#!/bin/bash
# factors_regen_check.sh [Z...] — X14「計算の再現性」: 出荷 JSON を同じ生成器で作り直して
#                                  **byte 同一**であることを確認する (260816Cl 新設)
#
# ⚠ 揮発情報 (時刻・秒・ホスト) は runlog/ に分離してあるので、本体 JSON は byte 同一になるはず。
#   同一でなければ (a) ソース指紋が違う (別 commit / dirty) か (b) ソルバが非決定論的、のどちらか。
#   (a) は generator_source_sha256 の比較で切り分ける。
# ⚠ 既定の元素 = H (1) + 多殻の中量級 Fe (26)。重元素 (Au) は 1 本 1 時間級なので必要時だけ。
# 使い方: bash tools/factors_regen_check.sh            # 1 26
#         bash tools/factors_regen_check.sh 1 26 79
#   環境変数 PROD (既定 src/prod_factors_v1) / REPO (生成器を走らせる checkout。既定 = このリポ)
set -eu
cd "$(dirname "$0")/.."
PROD=${PROD:-src/prod_factors_v1}
REPO=${REPO:-$(pwd)}
ZS=${*:-"1 26"}
tmp=$(cygpath -m "$(mktemp -d)")   # julia (Windows) に渡すので混合形式のパスへ
echo "=== X14 再現性: Z = $ZS / 生成器 = $REPO / 出荷 = $PROD ==="
ng=0
for z in $ZS; do
  f=$(printf 'SF_Z%03d.json' "$z")
  ( cd "$REPO" && julia --startup-file=no -t 1 src/gen_factors.jl "$z" --out "$tmp" > "$tmp/$f.log" 2>&1 ) || {
    echo "[NG] Z=$z の再生成に失敗"; tail -3 "$tmp/$f.log"; ng=$((ng+1)); continue; }
  fp_ship=$(python -c "import json;print(json.load(open('$PROD/$f',encoding='utf-8'))['generator_source_sha256'])")
  fp_new=$(python -c "import json;print(json.load(open('$tmp/$f',encoding='utf-8'))['generator_source_sha256'])")
  if cmp -s "$PROD/$f" "$tmp/$f"; then
    echo "  Z=$z: byte 同一 ✅ (指紋 ${fp_ship:0:12})"
  else
    ng=$((ng+1))
    if [ "$fp_ship" != "$fp_new" ]; then
      echo "  Z=$z: 差あり — ただし**ソース指紋が違う** (出荷 ${fp_ship:0:12} / 再生成 ${fp_new:0:12})。同じ commit の checkout で再実行すること"
    else
      echo "  Z=$z: 差あり ⚠ 同一指紋なのに byte が違う → ソルバの非決定論。要調査"
      python - "$PROD/$f" "$tmp/$f" <<'PY'
import json,sys
a=json.load(open(sys.argv[1],encoding='utf-8')); b=json.load(open(sys.argv[2],encoding='utf-8'))
for k in ('f_x','f_e_A'):
    d=max(abs(x-y) for x,y in zip(a[k],b[k])); print('   max|Δ%s| = %.3e' % (k,d))
for k in ('norm_correction','n_electrons_raw'):
    print('   %s: %r vs %r' % (k,a[k],b[k]))
PY
    fi
  fi
done
rm -rf "$tmp"
[ "$ng" -eq 0 ] && echo "X14 再現性: ALL PASS" || { echo "X14 再現性: $ng 件 NG"; exit 1; }
