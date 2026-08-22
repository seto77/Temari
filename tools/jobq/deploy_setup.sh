#!/bin/bash
# deploy_setup.sh [ROOT] [--spool SPOOL] [--dry-run] — tools/jobq/ の配布物を共有へ置く (PROTOCOL.md §1.5)
#
#   ROOT  = $1 (オプション以外の最初の引数) > $JOBQ_ROOT > //10.31.108.5/jobq   ← 人が開く共有の直下
#   SPOOL = --spool > $JOBQ_SPOOL > $ROOT/spool                                  ← 機械が書くもの全部
#
#   共有の直下 (人が見る場所) へ 3 ファイル: register.cmd / unregister.cmd / README.txt (**CRLF**)
#   ROOT/setup/ へ 7 ファイル (**LF**): PIN.json bootstrap.ps1 nastest.ps1 queuectl.jl reaper.sh
#                                       worker.conf.template worker.sh
#   ROOT/code/ は空のまま作る (中身は pack_code.sh が入れる。§1.4)
#   SPOOL/ の骨組み: queue queue/.tmp running results done failed control hosts campaigns
#
#   - どれかが無い・改行が規則どおりでない → 何もせず exit 1
#     (ワーカーは SETUP_SHA256 で同期するので、欠けた組・壊れた組を配らない)
#   - 各ファイルは .tmp.<name>.<pid> に書いてから rename。宛先を読み直して hash を照合する
#   - SETUP_SHA256 は**最後に**書く。覆うのは setup/ の 7 ファイルだけ (code/ と spool/ は含めない)
#   - ROOT 自体は作らない (共有が見えていない状態で /c 直下に掘らないため)
#   テスト用: JOBQ_SETUP_SRC=<dir> で配布元を差し替えられる (既定は本スクリプトのあるディレクトリ)
set -u

SETUP_FILES="PIN.json bootstrap.ps1 nastest.ps1 queuectl.jl reaper.sh worker.conf.template worker.sh"  # LC_ALL=C 名前順
ROOT_FILES="register.cmd:register.cmd unregister.cmd:unregister.cmd share_README.txt:README.txt"       # src:dst (CRLF)
SPOOL_DIRS="queue queue/.tmp running results done failed control hosts campaigns"
OBSOLETE_DIRS="leases running/.reaping"   # 前版の名残 (2026-08-21 に廃止。消さずに知らせるだけ)
CR=$(printf '\r')

dry=0; root=""; spool=""
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) dry=1; shift ;;
    --spool)   [ $# -ge 2 ] || { printf 'deploy_setup: --spool に値が無い\n' >&2; exit 1; }; spool=$2; shift 2 ;;
    -h|--help) sed -n '2,18p' "$0"; exit 0 ;;
    --*)       printf 'deploy_setup: 未知の引数 %s\n' "$1" >&2; exit 1 ;;
    *)         if [ -z "$root" ]; then root=$1; shift; else printf 'deploy_setup: 引数が多すぎる (%s)\n' "$1" >&2; exit 1; fi ;;
  esac
done
[ -n "$root" ] || root=${JOBQ_ROOT:-//10.31.108.5/jobq}
root=${root%/}
[ -n "$spool" ] || spool=${JOBQ_SPOOL:-$root/spool}
spool=${spool%/}
src=${JOBQ_SETUP_SRC:-$(cd "$(dirname "$0")" && pwd)}

# --- 改行の検査 (LF のものは CR を 1 つも含まない / CRLF のものは全行が CRLF) --------
is_lf_only() { ! grep -qU "$CR" "$1"; }   # -U: Git for Windows の grep は text モードで行末の CR を剥がす (実測 3.0)
is_crlf() {
  local total withcr ncr nlf
  total=$(grep -cU '' "$1" 2>/dev/null || printf 0)
  withcr=$(grep -cU "$CR\$" "$1" 2>/dev/null || printf 0)
  ncr=$(tr -dc "$CR" < "$1" | wc -c); nlf=$(tr -dc '\n' < "$1" | wc -c)
  [ "$total" -gt 0 ] && [ "$total" = "$withcr" ] && [ "$ncr" = "$nlf" ]
}

# --- 配布元の検査 (欠け・改行) ---------------------------------------------------
bad=0
for f in $SETUP_FILES; do
  if [ ! -f "$src/$f" ]; then
    printf 'NG  setup/%s が無い (%s)\n' "$f" "$src" >&2; bad=1
  elif ! is_lf_only "$src/$f"; then
    printf 'NG  %s に CR が含まれる — LF に直してから配る (PROTOCOL §11.2 / .gitattributes)\n' "$f" >&2; bad=1
  fi
done
for pair in $ROOT_FILES; do
  s=${pair%%:*}; d=${pair##*:}
  if [ ! -f "$src/$s" ]; then
    printf 'NG  %s が無い (%s)\n' "$s" "$src" >&2; bad=1
  elif ! is_crlf "$src/$s"; then
    printf 'NG  %s は CRLF でなければならない (cmd.exe / Notepad が読む。宛先 %s) — 全行が CRLF か確かめる\n' "$s" "$d" >&2; bad=1
  fi
done
[ $bad -eq 0 ] || { printf 'deploy_setup: 配布元に問題があるので何もしない\n' >&2; exit 1; }

if [ ! -d "$root" ]; then
  printf 'deploy_setup: ROOT %s が無い (共有が見えていない?)。ROOT 自体は作らない\n' "$root" >&2; exit 1
fi

sha_of() { if [ -f "$1" ]; then sha256sum "$1" | cut -c1-64; else printf '%s' "-"; fi; }

note=""; [ $dry -eq 1 ] && note="  (dry-run: 何も書かない)"
printf 'deploy_setup: %s\n' "$note"
printf '  配布元 : %s\n' "$src"
printf '  ROOT   : %s   (人が開く場所: register.cmd / unregister.cmd / README.txt / setup/ / code/ / spool/)\n' "$root"
printf '  SPOOL  : %s\n' "$spool"

# --- before の hash と差分の表示 ----------------------------------------------
changed=0
for pair in $ROOT_FILES; do
  s=${pair%%:*}; d=${pair##*:}
  a=$(sha_of "$src/$s"); b=$(sha_of "$root/$d")
  if [ "$b" = "-" ]; then st=new; changed=$((changed+1))
  elif [ "$a" = "$b" ]; then st=same
  else st=changed; changed=$((changed+1)); fi
  printf '  %-8s %-24s <- %s\n' "$st" "$d" "$s"
done
for f in $SETUP_FILES; do
  a=$(sha_of "$src/$f"); b=$(sha_of "$root/setup/$f")
  if [ "$b" = "-" ]; then st=new; changed=$((changed+1))
  elif [ "$a" = "$b" ]; then st=same
  else st=changed; changed=$((changed+1)); fi
  printf '  %-8s setup/%-18s %s\n' "$st" "$f" "${a:0:16}"
done
new_sha=$(cd "$src" && sha256sum $SETUP_FILES)     # 名前順 (SETUP_FILES がそう並んでいる)、LF
old_sha=$(cat "$root/setup/SETUP_SHA256" 2>/dev/null || printf '')
if [ "$new_sha" = "$old_sha" ]; then printf '  %-8s setup/SETUP_SHA256\n' same; else printf '  %-8s setup/SETUP_SHA256\n' changed; fi
for d in setup code spool; do
  t=$root/$d; [ "$d" = spool ] && t=$spool
  [ -d "$t" ] || printf '  %-8s %s/\n' mkdir "$t"
done
for d in $SPOOL_DIRS; do
  [ -d "$spool/$d" ] || printf '  %-8s %s/%s/\n' mkdir "$spool" "$d"
done

# --- 前版の名残を知らせる (消さない) ---------------------------------------------
for d in $OBSOLETE_DIRS; do
  [ -d "$spool/$d" ] && printf '  ⚠ 前版の名残 %s/%s/ がある (2026-08-21 に廃止。中身を確かめてから人が消す)\n' "$spool" "$d"
done
for d in queue running results done failed control hosts campaigns leases; do
  if [ -d "$root/$d" ] && [ "$root/$d" != "$spool/$d" ]; then
    printf '  ⚠ 旧配置 %s/%s/ が共有の直下にある (今は %s/%s/。人が中身を確かめてから消す)\n' "$root" "$d" "$spool" "$d"
  fi
done

[ $dry -eq 1 ] && exit 0

# --- 骨組み ---------------------------------------------------------------------
mk() { mkdir -p "$1" || { printf 'deploy_setup: mkdir %s に失敗\n' "$1" >&2; exit 1; }; }
mk "$root/setup"; mk "$root/code"; mk "$spool"
for d in $SPOOL_DIRS; do mk "$spool/$d"; done

# --- 配布物 (tmp → rename、宛先を読み直して照合) --------------------------------
place() {   # $1 = 配布元, $2 = 宛先
  local from=$1 to=$2 tmp
  tmp="$(dirname "$to")/.tmp.$(basename "$to").$$"
  if ! { cp "$from" "$tmp" && mv -f "$tmp" "$to"; }; then
    rm -f "$tmp"; printf 'deploy_setup: %s の配置に失敗\n' "$to" >&2; exit 1
  fi
  if [ "$(sha_of "$to")" != "$(sha_of "$from")" ]; then
    printf 'deploy_setup: %s の宛先 hash が配布元と違う (切れたコピー?) — SETUP_SHA256 は書かない\n' "$to" >&2; exit 1
  fi
}
for pair in $ROOT_FILES; do
  s=${pair%%:*}; d=${pair##*:}
  place "$src/$s" "$root/$d"
done
for f in $SETUP_FILES; do place "$src/$f" "$root/setup/$f"; done

# --- 目印は最後 (同期中のワーカーが半端な組を掴んでも hash 不一致で次のループに直る) ---
tmp="$root/setup/.tmp.SETUP_SHA256.$$"
if ! { printf '%s\n' "$new_sha" > "$tmp" && mv -f "$tmp" "$root/setup/SETUP_SHA256"; }; then
  rm -f "$tmp"; printf 'deploy_setup: SETUP_SHA256 の配置に失敗\n' >&2; exit 1
fi
if [ "$(cat "$root/setup/SETUP_SHA256")" != "$new_sha" ]; then
  printf 'deploy_setup: SETUP_SHA256 の読み直しが一致しない\n' >&2; exit 1
fi
printf 'deploy_setup: 完了 (%d ファイル更新)。SETUP_SHA256 = %s\n' "$changed" "$(sha256sum "$root/setup/SETUP_SHA256" | cut -c1-16)"
exit 0
