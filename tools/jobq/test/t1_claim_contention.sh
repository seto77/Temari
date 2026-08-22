#!/bin/bash
# t1_claim_contention.sh [ROOT] [N] [R] [PRIM] — 受け入れ試験 T1: 同じ票を N 個の claimer が同時に取り合う
#   (設計書 §6 T1 / PROTOCOL §4 「元ファイルへの rename が成功した者だけが所有する」)
#   ROOT = $1 (既定 = scratch の t1/) — **共有ルート**。票は `ROOT/spool/queue/`、claim は `ROOT/spool/running/`
#       (2026-08-21 の配置変更。`JOBQ_SPOOL` を渡せば spool だけを別の場所にできる)
#   N = claimer 数 (既定 16)、R = 回数 (既定 50)、
#   PRIM = claim の原始操作 (既定 mv-verify = PROTOCOL §4 の CLAIM そのもの):
#     mv         queue/<base>.json → running/<base>.<owner>.json を `mv` (rename) するだけ
#     mv-verify  mv の後 0.5 s 置いて「宛先がある・元が無い」を読み直す (§4 の `rename_settled`。連鎖 rename で
#                先を越されていたら負けと判定する。**これが実装の規則**)
#     mkdir      running/<base>.lock/ を `mkdir` (CreateDirectory の排他)
#     noclobber  running/<base>.claim を `set -o noclobber` で作る (O_EXCL create の排他)
#   各回: queue/.tmp に票を書いて queue/ へ rename → go ファイルで N 個の claimer を同時に放す →
#         それぞれ別の owner 名で取りに行く → 勝者 (成功を信じた者) を数える。
#   合格 = 毎回 勝者がちょうど 1 / (mv 系) running/ の票が元と byte 一致 / queue/ に残骸なし / 余計なファイルなし。
#
#   ⚠⚠ 2026-08-20 実測 (この PC、ローカル NTFS): **素の mv は排他ではない** — Win32 の rename は「パスで開いて
#     ハンドルで改名」なので、最初の rename の前に開いた者は全員成功する (連鎖 rename。16 並列で 1 回に 2〜16 人が
#     「勝つ」)。Cygwin の mv でも Python の os.rename でも同じ。mkdir・noclobber・mv-verify は 0 失敗。
#   ⚠ 2026-08-20 深夜、実 NAS `//10.31.108.5/jobq` では **毎ラウンド勝者ちょうど 1** だった (PROTOCOL §11.1)。
#     ローカル NTFS の連鎖 rename は共有が SMB でない経路 (テスト・将来のローカル run) で効くので、
#     実装は mv-verify を使う。**この事実は再導出しない** — 測り直すのは NAS の設定を変えたときだけ。
#   ⚠ 本番 NAS で試すときは ROOT に //10.31.108.5/jobq/t1 のような**専用サブディレクトリ**を渡す
#     (campaign 名 jobq_t1 の票しか作らないが、本番の spool/ を共有しないこと)。
#   ⚠ 全台から同時に叩く T1 (設計書) は、各 PC でこのスクリプトを同じ ROOT・別 N で回し、
#     各 PC の勝者数の合計 = R になることを目で確かめる (本スクリプトは 1 PC 内の並行性だけを見る)。
set -u
SCRATCH_DEFAULT=${JOBQ_TEST_SCRATCH:-/c/Users/seto/AppData/Local/Temp/claude/c--Users-seto-source-repos-Temari/65d3ab9f-e469-4550-9c21-5d0e4b61f0d0/scratchpad/jobq_test}
root=${1:-$SCRATCH_DEFAULT/t1}
n=${2:-16}
rounds=${3:-50}
prim=${4:-mv-verify}
spool=${JOBQ_SPOOL:-$root/spool}
case "$root" in
  //10.31.108.5/jobq|//10.31.108.5/jobq/)
    printf 'T1: 本番 ROOT 直下では回さない (//10.31.108.5/jobq/t1 のような専用サブディレクトリを渡す)\n' >&2; exit 1 ;;
esac
case "$spool" in
  //10.31.108.5/jobq/spool|//10.31.108.5/jobq/spool/*)
    printf 'T1: 本番 spool では回さない\n' >&2; exit 1 ;;
esac
case "$prim" in mv|mv-verify|mkdir|noclobber) ;; *) printf 'T1: PRIM は mv / mv-verify / mkdir / noclobber (%s)\n' "$prim" >&2; exit 1 ;; esac

q="$spool/queue"; run="$spool/running"; win="$root/.t1_win"; go="$root/.t1_go"
mkdir -p "$q/.tmp" "$run" "$win" || exit 1
rm -rf "$q"/jobq_t1_* "$q/.tmp"/* "$run"/jobq_t1_* "$win"/* "$go" 2>/dev/null

claimer() {   # $1 = round (6 桁), $2 = claimer index
  local base="jobq_t1_$1.e001" i=$2 dst
  dst="$run/$base.t1-claimer-$i-s0-b1.json"
  while [ ! -e "$go" ]; do sleep 0.02; done
  case "$prim" in
    mv)        mv "$q/$base.json" "$dst" 2>/dev/null && : > "$win/$1.$i" ;;
    mv-verify) if mv "$q/$base.json" "$dst" 2>/dev/null; then
                 sleep 0.5
                 { [ -e "$dst" ] && [ ! -e "$q/$base.json" ]; } && : > "$win/$1.$i"
               fi ;;
    mkdir)     mkdir "$run/$base.lock" 2>/dev/null && : > "$win/$1.$i" ;;
    noclobber) ( set -o noclobber; : > "$run/$base.claim" ) 2>/dev/null && : > "$win/$1.$i" ;;
  esac
}

fail=0; t0=$(date +%s); multi=0
for r in $(seq 1 "$rounds"); do
  r6=$(printf '%06d' "$r"); base="jobq_t1_$r6.e001"
  printf '{"schema":1,"campaign":"jobq_t1","jobseq":%d,"claim_epoch":1,"task":"jobq.noop","code_sha256":"","code_commit":"","args":{"seconds":0},"created_utc":"%s","issued_by":"t1"}\n' \
    "$r" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$q/.tmp/$base.json"
  mv "$q/.tmp/$base.json" "$q/$base.json" || { printf 'round %d: 投入に失敗\n' "$r"; fail=1; break; }
  orig=$(sha256sum "$q/$base.json" | cut -c1-64)
  rm -f "$go"
  pids=""
  for i in $(seq 1 "$n"); do claimer "$r6" "$i" & pids="$pids $!"; done
  sleep 0.3            # 全 claimer が go 待ちに入るまで
  : > "$go"
  wait $pids
  winners=$(ls "$win" 2>/dev/null | grep -c "^$r6\.")
  ok=1; why=""
  [ "$winners" -eq 1 ] || { ok=0; why="$why winners=$winners[$(ls "$win" | grep "^$r6\." | sed "s/^$r6\.//" | tr '\n' ' ' | sed 's/ $//')]"; }
  [ "$winners" -gt 1 ] && multi=$((multi+1))
  case "$prim" in
    mv|mv-verify)
      claimed=$(ls "$run" 2>/dev/null | grep -c "^$base\.")
      [ "$claimed" -eq 1 ] || { ok=0; why="$why running=$claimed"; }
      [ -e "$q/$base.json" ] && { ok=0; why="$why queue-leftover"; }
      if [ "$claimed" -eq 1 ]; then
        f=$(ls "$run"/"$base".*.json | head -1)
        fo=$(basename "$f" | sed "s/^$base\.t1-claimer-//; s/-s0-b1\.json$//")
        why="$why final-owner=$fo"
        [ "$(sha256sum "$f" | cut -c1-64)" = "$orig" ] || { ok=0; why="$why content-differs"; }
      fi ;;
    mkdir|noclobber)
      rm -f "$q/$base.json" ;;   # 票自体は触らない原始操作なので片付ける
  esac
  if [ $ok -eq 0 ]; then printf 'round %d: FAIL%s\n' "$r" "$why"; fail=1; fi
done
rm -f "$go"
left_q=$(ls "$q" 2>/dev/null | grep -vc '^\.tmp$')
left_tmp=$(ls "$q/.tmp" 2>/dev/null | grep -c .)
n_run=$(ls "$run" 2>/dev/null | grep -c "^jobq_t1_"); n_win=$(ls "$win" 2>/dev/null | grep -c .)
printf 'T1 claim contention (%s): ROOT=%s SPOOL=%s\n' "$prim" "$root" "$spool"
printf '  N=%d R=%d  %d s  winners=%d (期待 %d)  勝者 2 人以上の回=%d  running=%d  queue 残=%d  .tmp 残=%d\n' \
  "$n" "$rounds" "$(( $(date +%s) - t0 ))" "$n_win" "$rounds" "$multi" "$n_run" "$left_q" "$left_tmp"
[ "$n_win" -eq "$rounds" ] || fail=1
[ "$n_run" -eq "$rounds" ] || fail=1
{ [ "$left_q" -eq 0 ] && [ "$left_tmp" -eq 0 ]; } || fail=1
if [ $fail -eq 0 ]; then printf 'T1 (%s): PASS\n' "$prim"; exit 0; else printf 'T1 (%s): FAIL\n' "$prim"; exit 1; fi
