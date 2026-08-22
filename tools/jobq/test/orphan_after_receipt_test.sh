#!/bin/bash
# orphan_after_receipt_test.sh — §7 の「REISSUE の前に票自身の receipt を見直す」の負のテスト (2026-08-22)。
#
# 実演する主張:
#   1. 回収済み (orphan) の票に **その票自身の receipt** が付いていれば、reaper は **再投入しない**
#      (outcome = already_finished)。
#   2. receipt が無ければ、同じ状況で **再投入する** (epoch+1 が queue/ に出る)。
#      ⇒ 1 で再投入されなかったのは receipt を見たからであって、他の理由ではない。
#
# なぜ要るか: pass() は REAP の直前にも has_receipt を見るが、**retry_orphans は後の走査で
#   finish_orphan を呼び直す**ので、その間に本人が完走して receipt を書きうる (窓は 1 周期 = 300 s 以上)。
#   見落とすと、完了済みの計算をやり直すだけでなく、**別 CPU の再計算はバイト列が変わる** (§6.5) ため
#   publish が dup として弾き、failed/<c>/dup/ に**偽の「不一致」の証拠**が残る。
#
# 使い方:  bash tools/jobq/test/orphan_after_receipt_test.sh
#   JOBQ_TEST_SCRATCH で scratch の場所を変えられる。本番 NAS (//...) と /c/jobq は拒否する。
set -u

SCRATCH=${JOBQ_TEST_SCRATCH:-/c/Users/seto/AppData/Local/Temp/claude/jobq_test}/orphan_after_receipt
case $SCRATCH in //*|/c/jobq*) echo "refusing to run against $SCRATCH" >&2; exit 2 ;; esac
ROOT=$SCRATCH/root; SPOOL=$ROOT/spool; LOCAL=$SCRATCH/local
here=$(cd "$(dirname "$0")" && pwd); REAPER=$(cd "$here/.." && pwd)/reaper.sh

PASSED=0; FAILED=0
check(){ if [ "$2" = "$3" ]; then PASSED=$((PASSED+1)); printf 'PASS  %s (%s)\n' "$1" "$2"
         else FAILED=$((FAILED+1)); printf 'FAIL  %s (期待 %s / 実際 %s)\n' "$1" "$3" "$2"; fi; }

CAMP=orphrec; OWNER=oar-test-s0-b1
LOG=$LOCAL/logs/reaper.log
run_once() { JOBQ_ROOT=$ROOT JOBQ_SPOOL=$SPOOL JOBQ_LOCAL=$LOCAL JOBQ_CLAIM_TIMEOUT=1 \
             JOBQ_REAPER_INTERVAL=1 JOBQ_SETTLE_SECONDS=0.1 bash "$REAPER" --once > /dev/null 2>&1; }

# $1 = jobseq (6 桁), $2 = "receipt" なら done/ に receipt も置く
setup_orphan() {
  local j=$1 base="${CAMP}_$1.e001"
  mkdir -p "$SPOOL/failed/$CAMP/orphan" "$SPOOL/done/$CAMP"
  cat > "$SPOOL/failed/$CAMP/orphan/$base.$OWNER.json" <<JSON
{ "schema": 1, "campaign": "$CAMP", "jobseq": $((10#$j)), "claim_epoch": 1, "task": "jobq.noop",
  "code_sha256": "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
  "code_commit": "", "args": { "seconds": 1 },
  "created_utc": "2026-08-22T00:00:00Z", "issued_by": "orphan_after_receipt_test" }
JSON
  if [ "$2" = receipt ]; then
    cat > "$SPOOL/done/$CAMP/$base.$OWNER.json" <<JSON
{ "schema": 1, "base": "$base", "owner": "$OWNER", "task": "jobq.noop",
  "outnames": ["${CAMP}_lane$1001.jsonl"], "manifest_sha256": [], "finished_utc": "2026-08-22T00:01:00Z" }
JSON
  fi
}

rm -rf "$SCRATCH"
mkdir -p "$SPOOL/running" "$SPOOL/queue" "$SPOOL/done" "$SPOOL/failed" "$SPOOL/hosts" "$LOCAL/state" "$LOCAL/logs"
: > "$SPOOL/hosts/idle-worker-s0-b1.status.json.tmp"   # hosts/ が使えることの材料 (§7 の全滅ガードを黙らせる)
printf '{ "worker_id": "idle-worker", "slot": 0, "boot_seq": 1, "tick": 1, "base": null }\n' \
  > "$SPOOL/hosts/idle-worker-s0-b1.status.json"
rm -f "$SPOOL/hosts/idle-worker-s0-b1.status.json.tmp"

setup_orphan 000001 receipt     # 本人が完走していた orphan
setup_orphan 000002 none        # receipt の無い orphan (対照)
run_once

echo "--- 1. 票自身の receipt がある orphan ---"
check "epoch+1 は queue/ に出ない"        "$([ -f "$SPOOL/queue/${CAMP}_000001.e002.json" ] && echo yes || echo no)" "no"
check "outcome = already_finished"        "$(grep -c 'already_finished' "$LOG")"                                      "1"
check "決着の sidecar は書かれた"          "$([ -f "$SPOOL/failed/$CAMP/orphan/${CAMP}_000001.e001.$OWNER.reason.json" ] && echo yes || echo no)" "yes"
check "元の DONE receipt は無傷"           "$(grep -c '2026-08-22T00:01:00Z' "$SPOOL/done/$CAMP/${CAMP}_000001.e001.$OWNER.json")" "1"
check "reaper は FAIL receipt を書いていない" "$([ -f "$SPOOL/failed/$CAMP/${CAMP}_000001.e001.$OWNER.json" ] && echo yes || echo no)" "no"

echo "--- 2. receipt の無い orphan (対照: 同じ走査で再投入される) ---"
check "epoch+1 が queue/ に出る"          "$([ -f "$SPOOL/queue/${CAMP}_000002.e002.json" ] && echo yes || echo no)" "yes"
check "outcome = reissued"                "$(grep -c 'REISSUE orphrec_000002.e002' "$LOG")"                          "1"

echo
printf 'orphan_after_receipt_test: PASS %d / FAIL %d   (ログ: %s)\n' "$PASSED" "$FAILED" "$LOG"
[ "$FAILED" -eq 0 ]
