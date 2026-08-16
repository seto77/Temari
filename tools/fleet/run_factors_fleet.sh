#!/usr/bin/env bash
# dataset-factors の出荷生成フリート (リポ保存版。2026-08-16 の走は c:/tmp/temari_factors_2026-08-16/ の
# 版で発火済み — その版は既存ファイルをファイル名だけで skip していた (codex 指摘 2 回目)。
# ⚠ 本版は **skip を gen_factors.jl に任せる** (existing_is_current が版・model_id・指紋・格子 SHA を
#   検査し、食い違えば stale/ へ退避する)。3 pass 後に 86 本の**検査済み完成品**が揃わなければ exit 1。
# ⚠ frozen commit の git worktree から走らせる。-t 1 (決定論)。
set -u
REPO="${REPO:-c:/tmp/temari_factors_2026-08-16/repo}"
TOOL="$REPO/src/gen_factors.jl"
OUT="${OUT:-c:/Users/seto/source/repos/Temari/src/prod_factors_v1}"
BASE="${BASE:-c:/tmp/temari_factors_2026-08-16}"
LANES="${LANES:-8}"
mkdir -p "$OUT" "$BASE/logs"

ELEMENTS=$(seq 86 -1 1)      # 重い順

run_one() {
  set -u
  z=$1
  tag=$(printf 'SF_Z%03d' "$z")
  # ⚠ ここでファイル名だけを見て return しない。生成器が中身を検査して skip / 退避する
  cd "$REPO" && julia --startup-file=no -t 1 "$TOOL" "$z" --out "$OUT" \
        >> "$BASE/logs/$tag.log" 2>&1
}
export -f run_one
export OUT TOOL REPO BASE

# 検査済み完成品の数 (生成器と同じ規則: 版 / model_id / 指紋 / 格子 SHA / 長さ)
count_current() {
  cd "$REPO" && julia --startup-file=no -t 1 -e '
    include("src/gen_factors.jl")
    r = FactorsRecipe(); n = 0
    for z in FACTORS_Z_RANGE
        p = joinpath(ARGS[1], factors_filename(z))
        isfile(p) || continue
        d = try parse_json_file(p) catch; nothing end
        d === nothing && continue
        ok = d["dataset_version"] == factors_dataset_version(r) && d["model_id"] == factors_model_id(r) &&
             d["generator_source_sha256"] == FACTORS_SOURCE_FINGERPRINT && length(d["f_x"]) == S_N_INTERVALS + 1 &&
             length(d["f_e_A"]) == S_N_INTERVALS + 1 && d["scf"]["converged"] === true
        ok && (n += 1)
    end
    println(n)' "$OUT" 2>/dev/null | tail -1
}

for pass in 1 2 3; do
  echo "=== pass $pass  $(date '+%F %T') ==="
  echo "$ELEMENTS" | xargs -P "$LANES" -I{} bash -c 'run_one {}'
  n=$(count_current)
  echo "=== pass $pass done: $n / 86 (検査済み完成品)  $(date '+%F %T') ==="
  [ "$n" -ge 86 ] && break
done
n=$(count_current)
if [ "$n" -ge 86 ]; then echo "=== fleet done 86/86  $(date '+%F %T') ==="; exit 0
else echo "=== fleet INCOMPLETE $n/86  $(date '+%F %T') ==="; exit 1; fi
