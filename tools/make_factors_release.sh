#!/bin/bash
# make_factors_release.sh [prod_dir] [out_dir] — dataset-factors の配布物を**決定論的に**作る
#                                                 (260816Cl 新設。F 用 make_dataset_release.sh の対)
#
# ⚠ **F dataset (`temari-dataset-vX.Y.Z`) とは別 family・別 release・別 DOI** (計画書 §7)。
#   名前は `temari-factors-vX.Y.Z`、release タグは `dataset-factors-vX.Y.Z`。
# ⚠ 決定論: --sort=name / --mtime 固定 / owner 0 / gzip -n。**2 回組んで同一 SHA** になることを
#   このスクリプト自身が確認する (X14 の「梱包の決定論」)。
# ⚠ 自己完結: schema・golden・manifest・実行可能な契約 (`temari_factors_contract.py`)・
#   README・LICENSE を同梱する。**runlog/ (時刻・秒・ホスト) は同梱しない** — 表本体は再生成で
#   byte 同一だが runlog は走ごとに違うので、入れると「再生成しても同じ archive」が成り立たない
#   (codex 指摘 2 回目)。runlog は MANIFEST.md の集計と、リポ外の運用記録として残す。
# ⚠ tar の mode も正規化する (`--mode`)。並び・mtime・owner だけでは別ホストで再現しない。
# ⚠ 出荷前提: prod_dir に 86 元素 + manifest.json + MANIFEST.md + schema/factors_golden_v1.json
#   が揃っていること (無ければ止まる)。
set -eu
prod=${1:-src/prod_factors_v1}
outdir=${2:-dist}
cd "$(dirname "$0")/.."

[ -f "$prod/manifest.json" ] || { echo "manifest.json が無い: $prod (tools/make_factors_manifest.jl)" >&2; exit 1; }
[ -f "$prod/MANIFEST.md" ] || { echo "MANIFEST.md が無い: $prod" >&2; exit 1; }
[ -f schema/factors_golden_v1.json ] || { echo "schema/factors_golden_v1.json が無い (X13 未凍結)" >&2; exit 1; }
ver=$(python -c "import json; print(json.load(open('$prod/manifest.json',encoding='utf-8'))['dataset_version'])")
n=$(ls "$prod"/SF_Z???.json | wc -l)
if [ "${FACTORS_RELEASE_DRYRUN:-0}" != "1" ]; then     # 開発時の梱包試験だけ緩める (出荷では使わない)
  [ "$ver" != "0.0.0-dev" ] || { echo "dataset_version が 0.0.0-dev。出荷処方の出力ではない" >&2; exit 1; }
  [ "$n" -eq 86 ] || { echo "元素数が 86 でない ($n)" >&2; exit 1; }
fi
name="temari-factors-v${ver}"
stage="$outdir/$name"

build_once() {
  rm -rf "$stage"
  mkdir -p "$stage/schema" "$stage/tools" "$stage/licenses"
  cp "$prod"/SF_Z???.json "$stage/"
  cp "$prod/MANIFEST.md" "$prod/manifest.json" "$stage/"
  cp schema/temari_factors_v1.schema.json schema/factors_golden_v1.json "$stage/schema/"
  cp tools/temari_factors_contract.py "$stage/tools/"
  cp licenses/CC-BY-4.0.txt licenses/MIT.txt "$stage/licenses/"
  cp docs/factors_release_LICENSE.md "$stage/LICENSE.md"
  cp docs/factors_release_README.md "$stage/README.md"
  mt=$(python -c "
import re
m=open('$prod/MANIFEST.md',encoding='utf-8').read()
d=re.search(r'20\d\d-\d\d-\d\d', m)
print(d.group(0) if d else '2026-01-01')")
  ( cd "$outdir" && tar --sort=name --mtime="$mt UTC00:00" --owner=0 --group=0 \
        --numeric-owner --mode='u=rw,go=r,a+X' --format=gnu -cf - "$name" | gzip -n -9 > "$name.tar.gz" )
  echo "$mt"
}

echo "=== $name を組む (prod=$prod, $n 元素) ==="
mt=$(build_once)
sha1=$( (cd "$outdir" && sha256sum "$name.tar.gz") )
# ---- X14: もう一度組んで同一バイトになることを確認 ----------------------------------
build_once > /dev/null
sha2=$( (cd "$outdir" && sha256sum "$name.tar.gz") )
[ "$sha1" = "$sha2" ] || { echo "⚠ 梱包が決定論的でない: $sha1 vs $sha2" >&2; exit 1; }
echo "$sha1" > "$outdir/$name.tar.gz.sha256"
echo "  $sha1"
echo "  mtime を $mt に固定。2 回組んで同一 SHA (梱包の決定論 OK)"

# ---- 自己検査: 展開してから manifest・schema・契約 (+golden) を通す -------------------
echo "=== 展開して検証 ==="
tmp=$(mktemp -d)
tar --force-local -xzf "$outdir/$name.tar.gz" -C "$tmp"
julia --startup-file=no tools/make_factors_manifest.jl "$tmp/$name" --verify || {
  echo "manifest の照合に失敗" >&2; rm -rf "$tmp"; exit 1; }
PYTHONUTF8=1 python - "$tmp/$name" <<'EOF' || { echo "schema 検証に失敗" >&2; rm -rf "$tmp"; exit 1; }
import json, sys, glob, os, jsonschema
d = sys.argv[1]
sch = json.load(open(os.path.join(d, "schema", "temari_factors_v1.schema.json"), encoding="utf-8"))
v = jsonschema.Draft202012Validator(sch)
bad = 0
for f in sorted(glob.glob(os.path.join(d, "SF_Z???.json"))):
    errs = list(v.iter_errors(json.load(open(f, encoding="utf-8"))))
    if errs:
        bad += 1; print("[NG]", os.path.basename(f), errs[0].message[:120])
print("schema: %d ファイル検証, NG %d" % (len(glob.glob(os.path.join(d, 'SF_Z???.json'))), bad))
sys.exit(1 if bad else 0)
EOF
PYTHONUTF8=1 python "$tmp/$name/tools/temari_factors_contract.py" "$tmp/$name" \
    --golden "$tmp/$name/schema/factors_golden_v1.json" --negative     $( [ "${FACTORS_RELEASE_DRYRUN:-0}" = "1" ] && echo --allow-dev ) > "$tmp/contract.log" || {
  echo "契約テストに失敗"; tail -20 "$tmp/contract.log"; rm -rf "$tmp"; exit 1; }
tail -3 "$tmp/contract.log"
# ---- 展開物に Julia QC (F1–F10) も掛ける。指紋は manifest に凍結した値を期待値にする ----------
fp=$(python -c "import json; print(json.load(open('$prod/manifest.json',encoding='utf-8'))['generator_source_sha256'])")
julia --startup-file=no tools/check_factor_tables.jl "$tmp/$name" --expect-fingerprint "$fp" \
      --golden "$tmp/$name/schema/factors_golden_v1.json" \
      $( [ "${FACTORS_RELEASE_DRYRUN:-0}" = "1" ] && echo --allow-dev ) > "$tmp/qc.log" || {
  echo "Julia QC に失敗"; tail -20 "$tmp/qc.log"; rm -rf "$tmp"; exit 1; }
tail -1 "$tmp/qc.log"
echo "  ✅ 展開物で manifest 照合・schema・契約テスト (golden + 負のミュータント + scipy) ・Julia QC が通った"
rm -rf "$tmp"
rm -rf "$stage"
