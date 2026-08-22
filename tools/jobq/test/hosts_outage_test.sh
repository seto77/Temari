#!/bin/bash
# hosts_outage_test.sh — §7 の「hosts/ 全滅ガード」の負のテスト (2026-08-22)。
#
# 実演する主張:
#   1. `hosts/` が**使えない** (読めも書けもしない) pass では、reaper は **strike を積まず REAP もしない**
#      (WARN を出す)。running/ だけ読める部分障害で、生きているフリートを一斉に回収しないため。
#   2. `hosts/` は**使えるが status が 1 つも無い** (運用者が掃除した・新しい共有・全台退役) 場合は、
#      **通常どおり回収する**。⇒ ガードは回収を永久に止めない。
#
# 2 が要る理由: 1 のガードだけを入れると、hosts/ を掃除した時点で死んだ claim が二度と回収されなく
#   なる。reaper は「ワーカーの生死」ではなく「共有の hosts/ が読み書きできる状態か」を判定するので、
#   1 つも読めないときは**自分で書いて読み直して**その 2 つを分ける。
#
# 使い方:  bash tools/jobq/test/hosts_outage_test.sh
#   JOBQ_TEST_SCRATCH で scratch の場所を変えられる。本番 NAS (//...) と /c/jobq は拒否する。
set -u

SCRATCH=${JOBQ_TEST_SCRATCH:-/c/Users/seto/AppData/Local/Temp/claude/jobq_test}/hosts_outage
case $SCRATCH in //*|/c/jobq*) echo "refusing to run against $SCRATCH" >&2; exit 2 ;; esac
ROOT=$SCRATCH/root; SPOOL=$ROOT/spool; LOCAL=$SCRATCH/local
here=$(cd "$(dirname "$0")" && pwd); REAPER=$(cd "$here/.." && pwd)/reaper.sh

PASSED=0; FAILED=0
ok()   { PASSED=$((PASSED + 1)); printf 'PASS  %s\n' "$1"; }
bad()  { FAILED=$((FAILED + 1)); printf 'FAIL  %s\n' "$1"; }
check(){ if [ "$2" = "$3" ]; then ok "$1 ($2)"; else bad "$1 (期待 $3 / 実際 $2)"; fi; }

CAMP=hostguard; BASE=${CAMP}_000001.e001; OWNER=hg-test-s0-b1
TICKET_NAME=$BASE.$OWNER.json

rm -rf "$SCRATCH"; mkdir -p "$SPOOL/running" "$SPOOL/queue" "$SPOOL/done" "$SPOOL/failed" "$LOCAL/state" "$LOCAL/logs"
cat > "$SPOOL/running/$TICKET_NAME" <<'JSON'
{
  "schema": 1,
  "campaign": "hostguard",
  "jobseq": 1,
  "claim_epoch": 1,
  "task": "jobq.noop",
  "code_sha256": "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
  "code_commit": "",
  "args": { "seconds": 1 },
  "created_utc": "2026-08-22T00:00:00Z",
  "issued_by": "hosts_outage_test"
}
JSON

LOG=$LOCAL/logs/reaper.log
run_once() { JOBQ_ROOT=$ROOT JOBQ_SPOOL=$SPOOL JOBQ_LOCAL=$LOCAL JOBQ_CLAIM_TIMEOUT=1 \
             JOBQ_REAPER_INTERVAL=1 JOBQ_SETTLE_SECONDS=0.1 bash "$REAPER" --once > /dev/null 2>&1; }
# orphan/ には票と .reason.json の 2 つが入る (§1.2)。票だけを数える。
orphans() { ls "$SPOOL/failed/$CAMP/orphan/" 2>/dev/null | grep -cv '[.]reason[.]json$'; }

echo "--- 1. hosts/ が使えない (ディレクトリの場所に通常ファイル = 読めも書けもしない) ---"
rm -rf "$SPOOL/hosts"; : > "$SPOOL/hosts"
run_once; sleep 2; run_once; sleep 2; run_once     # 初見 + 沈黙 2 周期 = 本来なら REAP に届く
check "ガードの WARN が出ている"          "$(grep -c 'hosts/ not readable as a whole' "$LOG")" "3"
check "STRIKE は 1 度も積まれていない"     "$(grep -c 'STRIKE' "$LOG")"                        "0"
check "REAP されていない"                  "$(grep -c 'REAP' "$LOG")"                          "0"
check "claim は running/ に残っている"     "$([ -f "$SPOOL/running/$TICKET_NAME" ] && echo yes || echo no)" "yes"
check "orphan は 1 つも無い"               "$(orphans)"                                        "0"
s=$(awk -F'\t' -v k="$BASE" '$1 == k { print $5 }' "$LOCAL/state/reaper.tsv" 2>/dev/null | tail -1)
check "状態ファイルの strikes は 0"        "${s:-0}"                                           "0"

echo "--- 2. hosts/ は書けるが、status はあって 1 つも読めない (ACL・排他ロックの署名) ---"
# ⚠ ここが seen の分岐: プローブは**成功する** (hosts/ は書ける) ので、プローブだけで判定すると
#   ガードが外れて一斉回収になる。「ファイルはあるのに 1 つも読めない」を先に障害と断定すること。
rm -f "$SPOOL/hosts"; mkdir -p "$SPOOL/hosts"; : > "$SPOOL/hosts/dead-worker-s0-b1.status.json"
run_once; sleep 2; run_once; sleep 2; run_once
check "ガードの WARN が増えた"             "$(grep -c 'hosts/ not readable as a whole' "$LOG")" "6"
check "まだ STRIKE は積まれていない"       "$(grep -c 'STRIKE' "$LOG")"                        "0"
check "まだ REAP されていない"             "$(grep -c 'REAP' "$LOG")"                          "0"
check "claim はまだ running/ にある"       "$([ -f "$SPOOL/running/$TICKET_NAME" ] && echo yes || echo no)" "yes"

echo "--- 3. hosts/ は使えて status が 1 つも無い (掃除された・全台退役) ---"
rm -rf "$SPOOL/hosts"; mkdir -p "$SPOOL/hosts"
run_once; sleep 2; run_once; sleep 2; run_once
check "ガードの WARN は増えていない"       "$(grep -c 'hosts/ not readable as a whole' "$LOG")" "6"
check "今度は STRIKE が積まれた"           "$([ "$(grep -c 'STRIKE' "$LOG")" -ge 1 ] && echo yes || echo no)" "yes"
check "今度は REAP された"                 "$([ "$(grep -c 'REAP' "$LOG")" -ge 1 ] && echo yes || echo no)"   "yes"
check "claim は running/ から消えた"       "$([ -f "$SPOOL/running/$TICKET_NAME" ] && echo yes || echo no)"   "no"
check "orphan の票が 1 つできた"           "$(orphans)"                                        "1"
check "回収の理由 sidecar がある"          "$([ -f "$SPOOL/failed/$CAMP/orphan/$BASE.$OWNER.reason.json" ] && echo yes || echo no)" "yes"
check "epoch+1 の票が queue/ にある"       "$([ -f "$SPOOL/queue/${CAMP}_000001.e002.json" ] && echo yes || echo no)" "yes"
check "プローブの残骸が hosts/ に無い"     "$(ls -a "$SPOOL/hosts" 2>/dev/null | grep -c 'reaper-probe')"      "0"

echo
printf 'hosts_outage_test: PASS %d / FAIL %d   (ログ: %s)\n' "$PASSED" "$FAILED" "$LOG"
[ "$FAILED" -eq 0 ]
