#!/usr/bin/env bash
# 端点切断の 2 段目 (指示書 2026-08-15 §2-4 / 事前登録 v1 §8)。
#   (a) 極端例 (原点側 = Au 79・Rn 86 / 遠方側 = H 1・Li 3・Na 11・Cs 55) を
#       **採用格子 dt/16 (stage 5) で deep** (r0/10・r0/100・rmax*1.5・rmax*2、tight τ/10、
#       元素単位ラダー) = 1 元素 5 SCF
#   (b) 事前登録 §4.1 の標本 14 元素を **粗い dt (stage 1) で deep** (縮み比を見る)
# 重複する 5 元素 (H/Na/Cs/Au/Rn) は dt/16 側を正とし、独立証拠とは数えない (codex)。
# ⚠ frozen worktree (`repo/`) の certify_endpoints.jl を使う (tool_sha を嘘にしない)。
# ⚠ 既存 JSON は **中身を検査してから** skip する (z / stage / deep 4 変種 / 全解収束 /
#   tool_sha256 が worktree の実物と一致)。食い違えば stale/ へ退避して作り直す (codex 指摘 2 回目)。
# ⚠ 3 pass 後に 20 本揃わなければ exit 1 (完走を装わない)。
# 見積り: Rn/Au の tight dt/16 は 1 本 4000〜6400 s (16 レーン競合下) × 5 ⇒ 6〜9 h。
set -u
REPO="c:/tmp/temari_factors_2026-08-16/repo"
TOOL="$REPO/tools/certify_endpoints.jl"
OUT="c:/tmp/temari_endpoints2_2026-08-16"
BASE="c:/tmp/temari_factors_2026-08-16"
LANES="${LANES:-8}"
mkdir -p "$OUT/logs" "$OUT/stale"
TOOL_SHA=$(sha256sum "$TOOL" | cut -c1-64)

# "z:stage" の並び。重いもの (stage 5) を先に
ITEMS="86:5 79:5 55:5 11:5 3:5 1:5 86:1 79:1 64:1 55:1 54:1 36:1 29:1 24:1 20:1 18:1 11:1 10:1 2:1 1:1"

# 既存 JSON が「この走の完成品」か。1 = current (skip してよい)、0 = 無い/壊れている/別物 (退避済み)
is_current() {
  local f=$1 z=$2 st=$3
  [ -f "$f" ] || return 1
  # ⚠ ヒアドキュメントを export -f した関数の中に置くと子シェルで再パースできない
  #   (2026-08-16 12:53 に "syntax error near fi" で判明 → 別ファイルの python に出した)
  if python "$BASE/ep2_is_current.py" "$f" "$z" "$st" "$TOOL_SHA"; then return 0; fi
  mv "$f" "$OUT/stale/$(basename "$f").$(date +%Y%m%dT%H%M%S).stale"
  echo "⚠ 既存 $(basename "$f") が現行の走の完成品でない → stale/ へ退避して作り直す"
  return 1
}
export -f is_current

run_one() {
  set -u
  z=${1%%:*}; st=${1##*:}
  tag=$(printf 'ep_z%03d_st%d_deep' "$z" "$st")
  is_current "$OUT/$tag.json" "$z" "$st" && return 0
  cd "$REPO" && julia --startup-file=no -t 1 "$TOOL" "$z" --stage "$st" --deep --out "$OUT" \
        > "$OUT/logs/$tag.log" 2>&1
}
export -f run_one
export OUT TOOL REPO TOOL_SHA BASE

count_current() {
  local n=0
  for it in $ITEMS; do
    z=${it%%:*}; st=${it##*:}
    tag=$(printf 'ep_z%03d_st%d_deep' "$z" "$st")
    is_current "$OUT/$tag.json" "$z" "$st" >/dev/null 2>&1 && n=$((n+1))
  done
  echo "$n"
}

for pass in 1 2 3; do
  echo "=== pass $pass  $(date '+%F %T') ==="
  echo "$ITEMS" | tr ' ' '\n' | xargs -P "$LANES" -I{} bash -c 'run_one {}'
  n=$(count_current)
  echo "=== pass $pass done: $n / 20 (中身を検査した完成品の数)  $(date '+%F %T') ==="
  [ "$n" -ge 20 ] && break
done
n=$(count_current)
if [ "$n" -ge 20 ]; then
  echo "=== ep2 fleet done: 20/20  $(date '+%F %T') ==="; exit 0
else
  echo "=== ep2 fleet INCOMPLETE: $n/20  $(date '+%F %T') ==="; exit 1
fi
