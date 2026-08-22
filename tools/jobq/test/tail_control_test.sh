#!/bin/bash
# tail_control_test.sh — Deep 終盤の slow-host standby 規則を、claim 前の境界で検査する。
# 壊れた policy / 自分の速度未登録は必ず fail-open。ちょうど閾値は止めず、strictly slower だけを止める。
set -u

SCRATCH=${JOBQ_TEST_SCRATCH:-/c/Users/seto/AppData/Local/Temp/codex-jobq-tail}
case $SCRATCH in //*|/c/jobq|/c/jobq/*) echo "refusing non-scratch path: $SCRATCH" >&2; exit 2 ;; esac
ROOT=$SCRATCH/root; SPOOL=$ROOT/spool; LOCAL=$SCRATCH/local
here=$(cd "$(dirname "$0")" && pwd); worker=$(cd "$here/.." && pwd)/worker.sh
wid=tail-slow-1234; camp=temari_sigma_deep
passed=0; failed=0
ok() { passed=$((passed+1)); printf 'PASS  %s\n' "$1"; }
bad() { failed=$((failed+1)); printf 'FAIL  %s\n' "$1"; }
check() { local d=$1; shift; if "$@" >/dev/null 2>&1; then ok "$d"; else bad "$d"; fi; }
run() { JOBQ_TAIL_RULE_TEST=1 JOBQ_ROOT=$ROOT JOBQ_SPOOL=$SPOOL JOBQ_LOCAL=$LOCAL bash "$worker" 0 2>/dev/null; }
rescue_run() { JOBQ_TAIL_RULE_TEST=rescue JOBQ_TAIL_RULE_TEST_SLEEP=1 JOBQ_ROOT=$ROOT JOBQ_SPOOL=$SPOOL JOBQ_LOCAL=$LOCAL bash "$worker" 0 2>/dev/null; }
policy() { printf '%s' "$1" > "$SPOOL/control/tail.json"; }
ticket() { : > "$SPOOL/queue/$1"; }
check_out() { local d=$1 re=$2 out; out=$(run); if [[ "$out" =~ $re ]]; then ok "$d"; else bad "$d ($out)"; fi; }

rm -rf "$SCRATCH"; mkdir -p "$SPOOL/control" "$SPOOL/queue" "$LOCAL"
printf 'JOBQ_ROOT=%s\nJOBQ_SPOOL=%s\nJOBQ_LOCAL=%s\nWORKER_ID=%s\nSLOTS=1\nTHREADS=1\n' "$ROOT" "$SPOOL" "$LOCAL" "$wid" > "$LOCAL/worker.conf"
ticket "${camp}_000001.e001.json"

check_out "sidecar 無しは fail-open" '^may_work=yes reason=$'
policy '{"schema":2,"campaign":"temari_sigma_deep","fleet_slots":2,"k_per_slot":1,"median_slowdown":1,"threshold_multiplier":1,"rescue_seconds":60,"slowdown_by_worker_id":{"tail-slow-1234":3}}'
check_out "未知 schema は fail-open" '^may_work=yes'
policy '{"schema":1,"campaign":"temari_sigma_deep","fleet_slots":2,"k_per_slot":1,"median_slowdown":1,"threshold_multiplier":1,"rescue_seconds":60,"slowdown_by_worker_id":{"other-1234":3}}'
check_out "自分の速度未登録は fail-open" '^may_work=yes'
policy '{"schema":1,"campaign":"temari_sigma_deep","fleet_slots":2,"k_per_slot":1,"median_slowdown":2,"threshold_multiplier":1.5,"rescue_seconds":60,"slowdown_by_worker_id":{"tail-slow-1234":3}}'
check_out "等速 (3.0 == cutoff) は claim 可" '^may_work=yes'
policy '{"schema":1,"campaign":"temari_sigma_deep","fleet_slots":2,"k_per_slot":1,"median_slowdown":2,"threshold_multiplier":1.5,"rescue_seconds":60,"slowdown_by_worker_id":{"tail-slow-1234":3.01}}'
check_out "残票 < K×slots かつ strictly slow は standby" '^may_work=no reason=control/tail: remaining=1 < 2; slowdown=3.01 > 3;'
# hook だけでなく main ループでも、claim より前に票を残すことを確認する。
JOBQ_MAX_IDLE_LOOPS=1 JOBQ_ROOT=$ROOT JOBQ_SPOOL=$SPOOL JOBQ_LOCAL=$LOCAL bash "$worker" 0 >/dev/null 2>&1
check "standby の main loop は Deep 票を claim しない" test -f "$SPOOL/queue/${camp}_000001.e001.json"
ticket "${camp}_000002.e001.json"
check_out "残票がちょうど K×slots なら claim 可" '^may_work=yes'
rm -f "$SPOOL/queue"/*
ticket "aaa_other_000001.e001.json"; ticket "${camp}_000001.e001.json"
check_out "先頭が別 campaign なら Deep policy は他 campaign を止めない" '^may_work=yes'
rm -f "$SPOOL/queue"/*; ticket "${camp}_000001.e001.json"
policy '{"schema":1,"campaign":"temari_sigma_deep","fleet_slots":2,"k_per_slot":1,"median_slowdown":2,"threshold_multiplier":1.5,"rescue_seconds":1,"slowdown_by_worker_id":{"tail-slow-1234":3.01}}'
out=$(rescue_run)
if [[ "$out" == 'may_work=yes reason=control/tail: rescue after '* ]]; then ok "rescue 秒後は fail-open で claim 可"; else bad "rescue 秒後は fail-open で claim 可 ($out)"; fi

printf '\ntail_control_test: PASS %d / FAIL %d\n' "$passed" "$failed"
[ "$failed" -eq 0 ]
