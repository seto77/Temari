#!/bin/bash
# sync_from_recipro.sh [--apply]
#
# ReciPro 側 (現時点の正本) から src/ へコードを同期する。既定は差分表示のみで、
# 上書きは --apply を付けたときだけ。P1 (層分割) が終わって Temari が正本に
# なったら、このスクリプトは向きを逆にするか削除する。
set -u
SRC="${TEMARI_RECIPRO_HANDOUT:-/c/Users/seto/source/repos/ReciPro/tools/IonizationGen/handout}"
DST="$(cd "$(dirname "$0")/../src" && pwd)"
FILES="ionization.jl ionization.py gen_production.jl bote_salvat.json reference_values.json"

[ -d "$SRC" ] || { echo "取り込み元が無い: $SRC"; exit 1; }
apply=0
[ "${1:-}" = "--apply" ] && apply=1

echo "取り込み元: $SRC"
echo "取り込み先: $DST"
if command -v git >/dev/null && git -C "$SRC" rev-parse --short HEAD >/dev/null 2>&1; then
  echo "元リポの HEAD: $(git -C "$SRC" rev-parse --short HEAD)"
  # 元リポに未コミット変更があると、何を取り込んだのか後から辿れない
  if [ -n "$(git -C "$SRC" status --porcelain -- $FILES)" ]; then
    echo "⚠ 取り込み元に未コミット変更があります。先に元リポでコミットしてください:"
    git -C "$SRC" status --short -- $FILES
    [ $apply -eq 1 ] && exit 1
  fi
fi

changed=0
for f in $FILES; do
  if [ ! -f "$DST/$f" ]; then
    echo "  [新規] $f"; changed=1
  elif ! cmp -s "$SRC/$f" "$DST/$f"; then
    echo "  [差分] $f ($(diff <(cat "$DST/$f") <(cat "$SRC/$f") | grep -c '^[<>]') 行)"
    changed=1
  fi
done
[ $changed -eq 0 ] && { echo "差分なし"; exit 0; }

if [ $apply -eq 1 ]; then
  for f in $FILES; do cp "$SRC/$f" "$DST/$f"; done
  echo "同期しました。src/IMPORT.md のコミット表を更新し、selftest を通してから commit すること"
else
  echo "(表示のみ。反映するには --apply)"
fi
