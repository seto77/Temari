#!/usr/bin/env bash
# tools/require_all_pass.sh LOG — selftest のログに `ALL PASS (` で始まる行があることを要求する (260922Cl、R2 の M1 の H)。
#
# CI の selftest の step は終了コードしか見ていなかった。selftest が 1 つも走らずに exit 0 で終わる形
# (CLI の分岐の置き場を変えて `selftest` の引数が届かない、wrapper が呼び忘れる) でも合格していた
# (`docs/notes/r2_module_repro_plan_2026-09-22.md` §2 の H)。jobq の検収 (PROTOCOL.md) と同じ条件にそろえる。
#
# EXIT 0 = 行がある / 1 = 無い / 2 = 使い方・ログが無い
set -u
[ $# -eq 1 ] || { echo "usage: $0 LOG" >&2; exit 2; }
[ -f "$1" ] || { echo "require_all_pass: no log file $1" >&2; exit 2; }
if grep -qE '^ALL PASS \(' "$1"; then
  grep -E '^ALL PASS \(' "$1" | tail -1
  exit 0
fi
echo "require_all_pass: no line starting with 'ALL PASS (' in $1 — selftest did not run to the end" >&2
exit 1
