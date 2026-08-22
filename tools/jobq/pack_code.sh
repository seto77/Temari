#!/bin/bash
# pack_code.sh <source-tree> [--out-root ROOT] [--name NAME] [--prod-fp HEX16] [--dry-run]
#   — コードツリーを内容アドレスの tar.gz に固めて ROOT/code/ へ置く (PROTOCOL.md §1.4)
#
# なぜ git clone ではないのか (2026-08-21 実測。PROTOCOL §1.4。再導出しない):
#   - ローカル main は origin より 26 コミット先行しており、固定したい commit が GitHub に無い。
#   - clone の改行は core.autocrlf 次第で、同じ commit でも指紋が割れる。書庫は実バイトを運ぶので起きない。
#   - 走行中の F v6 フリートのツリーを下の 1 行で固めて別の場所へ展開すると
#     generator_source_fingerprint が ce058cce4fe9b31d まで一致する (= 他 PC が合流できる)。0.7 MB。
#   - 同じ内容を 2 回固めると sha256 も同じ (決定論オプションが「同じコード ⇒ 同じ id」を成立させる)。
#   - ツリー全体を固めると atom_cache で 619 MB になる ⇒ パス一覧は明示する。
#
# ⚠ 書庫に .git は入らない (パス一覧の帰結) ⇒ この書庫で走った遠隔ワーカーの F_*.json は
#   generator_commit が "unknown" になる (gen_production.jl の _git_probe が exit 128 で落ちるため)。
#   digest → commit の対応は sidecar <name>-<sha16>.json の "commit" が正本。
#   ⚠ ここに commit を書いたファイル (PACKED_COMMIT 等) を混ぜても読ませられない —
#   読む側 (_git_probe) を直すと gen_production.jl が変わり、PRODUCTION_SOURCE_FINGERPRINT が
#   ce058cce4fe9b31d から動いて走行中のフリートに合流できなくなる (実測: コメント 1 行で f8d9a89cc3c33a4e)。
#
# 識別子の定義 (この 1 行。オプションを足したり削ったりしてはいけない):
#   tar -C <tree> --sort=name --mtime='UTC 2020-01-01' --owner=0 --group=0 --numeric-owner \
#       -cf - src tools Project.toml | gzip -n -9
#
# 出力: 標準出力に 64 桁の sha256 だけ (campaign にそのまま貼れる / SHA=$(pack_code.sh …) で受けられる)。
#       人向けの報告は標準エラーへ出す。
# 終了コード: 0 = 成功 (既に同じ digest があった場合も成功) / 1 = 使い方・配布元の欠陥 / 3 = 置けなかった
set -u

PATHS="src tools Project.toml"      # ★ 識別子の一部。変えるなら世代を分けること
CODE_NAME_RE='^[a-z][a-z0-9_-]{0,31}$'

usage() { sed -n '2,26p' "$0" >&2; }
err()   { printf 'pack_code: %s\n' "$1" >&2; }

name=""; root=""; tree=""; prodfp=""; dry=0
while [ $# -gt 0 ]; do
  case "$1" in
    --out-root) [ $# -ge 2 ] || { err "--out-root に値が無い"; exit 1; }; root=$2; shift 2 ;;
    --name)     [ $# -ge 2 ] || { err "--name に値が無い"; exit 1; };     name=$2; shift 2 ;;
    --prod-fp)  [ $# -ge 2 ] || { err "--prod-fp に値が無い"; exit 1; };  prodfp=$2; shift 2 ;;
    --dry-run)  dry=1; shift ;;
    -h|--help)  usage; exit 0 ;;
    --*)        err "未知の引数 $1"; exit 1 ;;
    *)          if [ -z "$tree" ]; then tree=$1; else err "引数が多すぎる ($1)"; exit 1; fi; shift ;;
  esac
done

[ -n "$tree" ] || { err "source-tree を指定すること"; usage; exit 1; }
[ -n "$name" ] || name=${JOBQ_CODE_NAME:-temari}
[ -n "$root" ] || root=${JOBQ_ROOT:-//10.31.108.5/jobq}
root=${root%/}

[[ "$name" =~ $CODE_NAME_RE ]] || { err "--name '$name' が規則 $CODE_NAME_RE に合わない (ファイル名になる)"; exit 1; }
if [ -n "$prodfp" ] && ! [[ "$prodfp" =~ ^[0-9a-f]{16}$ ]]; then
  err "--prod-fp '$prodfp' は 16 桁の 16 進でなければならない"; exit 1
fi

# --- 配布元の検査 ------------------------------------------------------------------
[ -d "$tree" ] || { err "source-tree が無い: $tree"; exit 1; }
tree_abs=$(cd "$tree" && pwd) || { err "source-tree に cd できない: $tree"; exit 1; }
case "$tree_abs" in *'"'*|*'\'*) err "source-tree のパスに \" か \\ が入っている ($tree_abs)"; exit 1 ;; esac

miss=0
for p in $PATHS; do
  [ -e "$tree_abs/$p" ] || { err "パス一覧の $p が $tree_abs に無い"; miss=1; }
done
[ $miss -eq 0 ] || { err "パス一覧が揃っていないので何もしない (固めるのは: $PATHS)"; exit 1; }

# --- 来歴 (git。強制はしない — 識別子は digest の方) --------------------------------
commit=""; dirty=false; gitnote=""
if git -C "$tree_abs" rev-parse --git-dir >/dev/null 2>&1; then
  commit=$(git -C "$tree_abs" rev-parse HEAD 2>/dev/null || printf '')
  st=$(git -C "$tree_abs" status --porcelain -uno 2>/dev/null || printf 'ERR')
  if [ "$st" = "ERR" ]; then
    gitnote="git status が読めなかった"; dirty=true
  elif [ -n "$st" ]; then
    dirty=true
    printf '\n' >&2
    err "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    err "!! 作業ツリーが dirty です (追跡ファイルに未コミットの変更がある)。"
    err "!! 記録する commit は <sha>-dirty になります。あとから同じ commit を checkout しても"
    err "!! この書庫のバイトは再現できません — 再現できるのは digest だけです。"
    err "!! 未コミットの変更:"
    printf '%s\n' "$st" | sed 's/^/!!   /' >&2
    err "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    printf '\n' >&2
  fi
  [ -n "$commit" ] || gitnote="HEAD が読めなかった (空の repo?)"
else
  gitnote="git 管理下ではない"
fi
commit_rec=""
if [ -z "$commit" ]; then
  dirty=true
  err "警告: $tree_abs の commit が特定できない ($gitnote) — code_commit は \"\" になる (来歴のみなので票は通る)"
elif [ "$dirty" = true ]; then
  commit_rec="$commit-dirty"
else
  commit_rec="$commit"
fi

# --- 固める (まず手元に。NAS には検証済みのバイトしか置かない) -----------------------
tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/pack_code.XXXXXX") || { err "一時ディレクトリを作れない"; exit 3; }
cleanup() { rm -rf "$tmpdir"; }
trap cleanup EXIT
local_tar="$tmpdir/pack.tar.gz"

# ⚠ tar のパスは /c/… 形式 (MSYS の GNU tar は "C:" をリモートホスト扱いする。PROTOCOL §11.2)。
#   -f - なので書き出しはシェルのリダイレクトで、tar 側にホスト解釈の余地は無い。
set -o pipefail
tar -C "$tree_abs" --sort=name --mtime='UTC 2020-01-01' --owner=0 --group=0 --numeric-owner \
    -cf - $PATHS 2>"$tmpdir/tar.err" | gzip -n -9 > "$local_tar"
rc=$?
set +o pipefail
if [ $rc -ne 0 ]; then
  err "tar/gzip が失敗した (rc=$rc)"; sed 's/^/    /' "$tmpdir/tar.err" >&2; exit 3
fi
[ -s "$tmpdir/tar.err" ] && { err "tar の警告:"; sed 's/^/    /' "$tmpdir/tar.err" >&2; }

sha=$(sha256sum "$local_tar" | cut -c1-64)
[[ "$sha" =~ ^[0-9a-f]{64}$ ]] || { err "sha256 が取れなかった"; exit 3; }
sha16=${sha:0:16}
bytes=$(wc -c < "$local_tar" | tr -d ' ')
# --force-local: TMPDIR が C:\… 形式に設定されていても GNU tar が "C:" をリモートホスト扱いしない (§11.2)
nfiles=$(tar --force-local -tzf "$local_tar" | wc -l | tr -d ' ')

arc="$name-$sha16.tar.gz"
meta="$name-$sha16.json"
dest="$root/code/$arc"
dest_meta="$root/code/$meta"

{
  printf '\n'
  printf 'pack_code: source-tree = %s\n' "$tree_abs"
  printf '           paths       = %s (%s entries)\n' "$PATHS" "$nfiles"
  printf '           commit      = %s%s\n' "${commit_rec:-(none)}" "$( [ -n "$gitnote" ] && printf ' [%s]' "$gitnote")"
  printf '           bytes       = %s\n' "$bytes"
  printf '           sha256      = %s\n' "$sha"
  printf '           -> %s\n' "$dest"
  printf '           %s\n' '⚠ .git は入らない ⇒ 遠隔生成の F_*.json は generator_commit="unknown"'
  printf '             digest → commit は %s の "commit" が正本\n' "$dest_meta"
  [ -n "$prodfp" ] && printf '           prod_fp     = %s (人が貼った値。検算はしない)\n' "$prodfp"
} >&2

if [ $dry -eq 1 ]; then
  if [ -e "$dest" ]; then err "dry-run: 宛先は既にある (上書きしない)"; else err "dry-run: 宛先はまだ無い"; fi
  err "dry-run: 何も書かなかった"
  printf '%s\n' "$sha"
  exit 0
fi

# --- 置く (tmp + rename。既にある digest は上書きしない) ----------------------------
if [ ! -d "$root" ]; then
  err "ROOT $root が無い (共有が見えていない?)。ROOT 自体は作らない"; exit 3
fi
mkdir -p "$root/code" || { err "mkdir $root/code に失敗"; exit 3; }

sha_of() { sha256sum "$1" 2>/dev/null | cut -c1-64; }

put_nooverwrite() {   # $1 = 手元のファイル, $2 = 宛先, $3 = 表示名 → 0 置いた / 10 先客と同一 / 1 失敗
  local from=$1 to=$2 label=$3 tmp
  if [ -e "$to" ]; then
    if [ "$(sha_of "$to")" = "$(sha_of "$from")" ]; then return 10; fi
    err "$label: 宛先が既にあって中身が違う — 上書きしない ($to)"; return 1
  fi
  tmp="$(dirname "$to")/.tmp.$(basename "$to").$$"
  cp "$from" "$tmp" || { rm -f "$tmp"; err "$label: 一時ファイルを書けない ($tmp)"; return 1; }
  if [ "$(sha_of "$tmp")" != "$(sha_of "$from")" ]; then
    rm -f "$tmp"; err "$label: 共有へのコピーが切れている (hash 不一致)"; return 1
  fi
  # ⚠ mv -n は宛先があっても rc=0 を返す (実 NAS 実測、PROTOCOL §11.1)。判定は宛先の読み直しで行う。
  mv -n "$tmp" "$to" 2>/dev/null
  if [ -e "$tmp" ]; then rm -f "$tmp"; fi
  if [ ! -e "$to" ]; then err "$label: rename 後に宛先が無い"; return 1; fi
  if [ "$(sha_of "$to")" != "$(sha_of "$from")" ]; then
    err "$label: 宛先の hash が自分のものと違う (競合した先客が別内容 / 切れたコピー) — $to"; return 1
  fi
  return 0
}

put_nooverwrite "$local_tar" "$dest" "書庫"; rc=$?
case $rc in
  0)  err "書庫を置いた: $dest" ;;
  10) err "同じ digest の書庫が既にある (同じバイト) — 置き直さない: $dest" ;;
  *)  exit 3 ;;
esac

# --- sidecar (先客があれば触らない: packed_utc などはその書庫を説明している方が正しい) ---
if [ -e "$dest_meta" ]; then
  err "sidecar が既にある — 触らない: $dest_meta"
  [ -n "$prodfp" ] && err "  (--prod-fp を入れ直したいなら $dest_meta を人が消してから再実行する)"
else
  paths_json=$(printf '%s' "$PATHS" | tr ' ' '\n' | sed 's/.*/"&"/' | paste -sd, -)
  {
    printf '{"schema":1,"name":"%s","commit":"%s","dirty":%s,' "$name" "$commit_rec" "$dirty"
    printf '"paths":[%s],"sha256":"%s","bytes":%s,' "$paths_json" "$sha" "$bytes"
    printf '"packed_utc":"%s","packed_by":"%s","source_tree":"%s"' \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${USERNAME:-$(id -un)}@$(hostname)" "$tree_abs"
    [ -n "$prodfp" ] && printf ',"prod_fp":"%s"' "$prodfp"
    printf '}\n'
  } > "$tmpdir/meta.json"
  put_nooverwrite "$tmpdir/meta.json" "$dest_meta" "sidecar"; rc=$?
  case $rc in
    0)  err "sidecar を置いた: $dest_meta" ;;
    10) err "同じ sidecar が既にある: $dest_meta" ;;
    *)  exit 3 ;;
  esac
fi

err "campaign に貼る値: --code-sha256 $sha $( [ -n "$commit_rec" ] && printf -- '--code-commit %s' "$commit_rec")"
printf '%s\n' "$sha"
exit 0
