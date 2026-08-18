#!/bin/bash
# make_dataset_release.sh [prod_dir] [out_dir] — データセット配布物を**決定論的に**作る
#                                                (260813Cl 追加。指示書 §7 (a))
#
# なぜ要るか: `prod*/` は .gitignore なので**データはリポジトリに入っていない**。
# コードだけが公開されていて、dataset v5.0.0 はどこにも配布されていない状態だった。
#
# ⚠⚠ **コードの release とは別立てにする** (作者指示 2026-08-13)。
#   ソフトウェア版 (`vX.Y.Z`) とデータセット版 (`dataset-vX.Y.Z`) は独立した版で、
#   CITATION.cff も「Generated datasets carry their own version, independent of
#   this software version」と宣言している。**同じ release に混ぜない。**
#
# ⚠ **決定論的に作る。**Zenodo (正本) と GitHub Release (ミラー) を
#   **同一バイト**にできなければ「どちらが本物か」を digest で言えない。そのために:
#     --sort=name        エントリ順をファイル名で固定 (readdir 順に依存させない)
#     --mtime            mtime を固定 (ビルド時刻を入れない)
#     --owner/--group 0  所有者を固定
#     gzip -n            gzip ヘッダにタイムスタンプと名前を入れない
#
# ⚠ **パッケージは自己完結**にする — データだけ配ると読み方が失われる。
#   schema・manifest・実行可能な契約 (`temari_contract.py`) を同梱する。
set -eu
prod=${1:-src/prod_v5_jl}
outdir=${2:-dist}
cd "$(dirname "$0")/.."

[ -f "$prod/manifest.json" ] || { echo "manifest.json が無い: $prod" >&2; exit 1; }
ver=$(python -c "import json,sys; print(json.load(open('$prod/manifest.json',encoding='utf-8'))['dataset_version'])")
name="temari-dataset-v${ver}"
stage="$outdir/$name"

echo "=== $name を組む (prod=$prod) ==="
rm -rf "$stage"
mkdir -p "$stage/schema" "$stage/tools" "$stage/licenses"

# ---- データ本体 + 正本 ------------------------------------------------------
cp "$prod"/F_*.json "$stage/"
cp "$prod/MANIFEST.md" "$prod/manifest.json" "$stage/"
# ---- 読み方 (これが無いとデータだけ残って規約が失われる) --------------------
cp schema/temari_dataset_v2.schema.json "$stage/schema/"
cp tools/temari_contract.py "$stage/tools/"
# ⚠ **ライセンスは 2 本立て** (作者決定 2026-08-10): データ本体は **CC-BY-4.0**、
# 同梱の loader (`temari_contract.py`) は **MIT**。どちらがどれに掛かるかを
# `LICENSE.md` に明記する — 1 枚の LICENSE を置くと、混在パッケージでは
# **どちらの条件で配られたのか読み手が判断できない**。
cp licenses/CC-BY-4.0.txt licenses/MIT.txt "$stage/licenses/"
cp docs/release/dataset_release_LICENSE.md "$stage/LICENSE.md"
cp docs/release/dataset_release_README.md "$stage/README.md"

n=$(ls "$stage"/F_*.json | wc -l)
echo "  チャネル $n / schema・manifest・contract・README 同梱"
echo "  ライセンス: データ = CC-BY-4.0 / 同梱 loader = MIT (LICENSE.md に明記)"

# ---- 決定論的な tar.gz ------------------------------------------------------
# mtime は**データセットの日付**に固定する (ビルド時刻ではない = 何度作っても同じ)
mt=$(python -c "
import json,re
m=open('$prod/MANIFEST.md',encoding='utf-8').read()
d=re.search(r'20\d\d-\d\d-\d\d', m)
print(d.group(0) if d else '2026-01-01')")
( cd "$outdir" && tar --sort=name --mtime="$mt UTC00:00" --owner=0 --group=0 \
      --numeric-owner --format=gnu -cf - "$name" | gzip -n -9 > "$name.tar.gz" )
sha=$( (cd "$outdir" && sha256sum "$name.tar.gz") )
echo "$sha" > "$outdir/$name.tar.gz.sha256"
echo "  $sha"
echo "  mtime を $mt に固定 (再実行しても同じバイトになる)"

# ---- 自己検査: 展開してから manifest と契約を通す ---------------------------
# ⚠ **配ったものをそのまま検査する。**元の prod ディレクトリを検査しても
#   「梱包で壊れていないか」の証明にならない。
echo "=== 展開して検証 ==="
tmp=$(mktemp -d)
tar -xzf "$outdir/$name.tar.gz" -C "$tmp"
julia +1.11 --startup-file=no tools/make_manifest.jl "$tmp/$name" --verify || {
  echo "manifest の照合に失敗" >&2; rm -rf "$tmp"; exit 1; }
PYTHONUTF8=1 python "$tmp/$name/tools/temari_contract.py" "$tmp/$name" > /dev/null || {
  echo "契約テストに失敗" >&2; rm -rf "$tmp"; exit 1; }
echo "  ✅ 展開物で manifest 照合と契約テストが通った"
rm -rf "$tmp"
rm -rf "$stage"

echo
echo "配布物: $outdir/$name.tar.gz  ($(du -h "$outdir/$name.tar.gz" | cut -f1))"
echo "⚠ Zenodo を正本、GitHub Release をミラーにする。DOI・版・digest を相互記載すること。"
